#!/usr/bin/env python3
import argparse
import gzip
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


PHARMCAT_DIR = Path("/opt/pharmcat")

PHARMCAT_JAR = PHARMCAT_DIR / "pharmcat-3.2.0-all.jar"
PHARMCAT_PREPROCESSOR = PHARMCAT_DIR / "pharmcat_vcf_preprocessor"
POSITIONS_SITES_VCF = PHARMCAT_DIR / "pharmcat_positions.sites.vcf.gz"
POSITIONS_REF_VCF = PHARMCAT_DIR / "pharmcat_positions.vcf.bgz"

THREADS = "5"
JAVA_MEM = "100g"
MIN_BCFTOOLS_VERSION = (1, 20)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Single-sample PharmCAT pipeline with flat output layout."
    )
    parser.add_argument("--cram", required=True, help="Input CRAM or BAM path")
    parser.add_argument("--reference", required=True, help="Reference FASTA path")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--sample-id", required=True, help="Sample ID prefix for all outputs")
    return parser.parse_args()


class PharmcatPipeline:
    def __init__(self, cram, reference_fasta, output_dir, sample_id):
        self.alignment = Path(cram).expanduser().resolve()
        self.reference_fasta = Path(reference_fasta).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.sample_id = sample_id.strip()

        if not self.sample_id:
            raise ValueError("sample_id cannot be empty")

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = self.output_dir / f"{self.sample_id}.pipeline.log"

    def log(self, message):
        line = str(message).rstrip()
        print(line, flush=True)
        with open(self.log_file, "a") as fh:
            fh.write(line + "\n")

    def run_cmd(self, cmd):
        cmd = [str(x) for x in cmd]
        self.log(f"[CMD] {' '.join(shlex.quote(x) for x in cmd)}")

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        if proc.stdout is not None:
            for line in proc.stdout:
                self.log(line.rstrip())

        ret = proc.wait()
        if ret != 0:
            raise RuntimeError(f"Command failed with exit code {ret}: {' '.join(cmd)}")

    def ensure_file(self, path, label):
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")

    def require_vcf_index(self, vcf_path, label):
        vcf_path = Path(vcf_path)
        tbi = Path(str(vcf_path) + ".tbi")
        csi = Path(str(vcf_path) + ".csi")
        if not tbi.exists() and not csi.exists():
            raise FileNotFoundError(
                f"{label} index not found. Expected one of: {tbi} or {csi}"
            )

    def require_tools(self):
        for tool in ["gatk", "samtools", "bgzip", "tabix", "java", "bcftools"]:
            if shutil.which(tool) is None:
                raise RuntimeError(f"Required tool not found in PATH: {tool}")

        self.ensure_file(PHARMCAT_JAR, "PharmCAT JAR")
        self.ensure_file(PHARMCAT_PREPROCESSOR, "PharmCAT preprocessor")
        self.ensure_file(POSITIONS_SITES_VCF, "Positions sites VCF")
        self.ensure_file(POSITIONS_REF_VCF, "Positions ref VCF")

        self.require_vcf_index(POSITIONS_SITES_VCF, "Positions sites VCF")
        self.require_vcf_index(POSITIONS_REF_VCF, "Positions ref VCF")

    def require_bcftools_version(self):
        res = subprocess.run(
            ["bcftools", "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        first_line = res.stdout.splitlines()[0].strip()
        parts = first_line.split()
        if len(parts) < 2:
            raise RuntimeError(f"Could not parse bcftools version: {first_line}")

        version_str = parts[1]
        nums = version_str.split(".")
        major = int(nums[0])
        minor = int(nums[1]) if len(nums) > 1 else 0

        if (major, minor) < MIN_BCFTOOLS_VERSION:
            raise RuntimeError(
                f"bcftools >= {MIN_BCFTOOLS_VERSION[0]}.{MIN_BCFTOOLS_VERSION[1]} required, found {version_str}"
            )

    def get_ref_dict_path(self):
        return self.reference_fasta.with_suffix(".dict")

    def ensure_reference_indexes(self):
        fai = Path(str(self.reference_fasta) + ".fai")
        ref_dict = self.get_ref_dict_path()

        if not fai.exists():
            self.log(f"[INFO] Creating FASTA index: {fai}")
            self.run_cmd(["samtools", "faidx", str(self.reference_fasta)])

        if not ref_dict.exists():
            self.log(f"[INFO] Creating sequence dictionary: {ref_dict}")
            self.run_cmd(
                [
                    "gatk",
                    "CreateSequenceDictionary",
                    "-R",
                    str(self.reference_fasta),
                    "-O",
                    str(ref_dict),
                ]
            )

    def ensure_alignment_index(self):
        if self.alignment.suffix == ".cram":
            idx1 = Path(str(self.alignment) + ".crai")
            idx2 = self.alignment.with_suffix(".crai")
            if not idx1.exists() and not idx2.exists():
                self.log(f"[INFO] Creating CRAI index: {self.alignment}")
                self.run_cmd(["samtools", "index", str(self.alignment)])

        elif self.alignment.suffix == ".bam":
            idx1 = Path(str(self.alignment) + ".bai")
            idx2 = self.alignment.with_suffix(".bai")
            if not idx1.exists() and not idx2.exists():
                self.log(f"[INFO] Creating BAI index: {self.alignment}")
                self.run_cmd(["samtools", "index", "-@", THREADS, str(self.alignment)])
        else:
            raise ValueError(f"Unsupported alignment type: {self.alignment}")

    def write_input_check(self):
        out = self.output_dir / f"{self.sample_id}.inputs.txt"
        with open(out, "w") as fh:
            fh.write(f"sample_id\t{self.sample_id}\n")
            fh.write(f"alignment\t{self.alignment}\n")
            fh.write(f"reference_fasta\t{self.reference_fasta}\n")
            fh.write(f"pharmcat_jar\t{PHARMCAT_JAR}\n")
            fh.write(f"pharmcat_preprocessor\t{PHARMCAT_PREPROCESSOR}\n")
            fh.write(f"positions_sites_vcf\t{POSITIONS_SITES_VCF}\n")
            fh.write(f"positions_ref_vcf\t{POSITIONS_REF_VCF}\n")
            fh.write(f"threads\t{THREADS}\n")
            fh.write(f"java_mem\t{JAVA_MEM}\n")

    def run_forcecall(self):
        raw_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.raw.vcf.gz"

        self.run_cmd([
            "gatk",
            "--java-options", f"-Xmx{JAVA_MEM}",
            "HaplotypeCaller",
            "--native-pair-hmm-threads", THREADS,
            "--alleles", str(POSITIONS_SITES_VCF),
            "-R", str(self.reference_fasta),
            "-I", str(self.alignment),
            "-O", str(raw_vcf),
            "-L", str(POSITIONS_SITES_VCF),
            "-ip", "20",
            "--max-mnp-distance", "1",
            "--output-mode", "EMIT_ALL_ACTIVE_SITES",
        ])

        self.log(f"[INFO] Indexing raw VCF: {raw_vcf}")
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(raw_vcf)])

        return raw_vcf

    def run_preprocessor(self, raw_vcf):
        self.run_cmd([
            str(PHARMCAT_PREPROCESSOR),
            "-vcf", str(raw_vcf),
            "-o", str(self.output_dir),
            "-bf", self.sample_id,
            "-refVcf", str(POSITIONS_REF_VCF),
            "-refFna", str(self.reference_fasta),
        ])
        return self.find_preprocessed_vcf()

    def find_preprocessed_vcf(self):
        exact_paths = [
            self.output_dir / f"{self.sample_id}.preprocessed.vcf.bgz",
            self.output_dir / f"{self.sample_id}.preprocessed.vcf.gz",
        ]
        for path in exact_paths:
            if path.exists():
                return path

        raise FileNotFoundError(
            f"Preprocessed VCF not found for sample_id '{self.sample_id}' in {self.output_dir}"
        )

    def open_text_auto(self, path):
        path = Path(path)
        if path.suffix in {".gz", ".bgz"}:
            return gzip.open(path, "rt")
        return open(path, "r")

    def fix_inf_qual(self, preprocessed_vcf):
        fixed_vcf = self.output_dir / f"{self.sample_id}.preprocessed.fixed.vcf.gz"
        self.log(f"[INFO] Fixing QUAL=inf: {preprocessed_vcf} -> {fixed_vcf}")

        with self.open_text_auto(preprocessed_vcf) as in_fh, open(fixed_vcf, "wb") as out_fh:
            proc = subprocess.Popen(
                ["bgzip", "-c"],
                stdin=subprocess.PIPE,
                stdout=out_fh,
                stderr=subprocess.PIPE,
                text=True,
            )

            if proc.stdin is None or proc.stderr is None:
                raise RuntimeError("Failed to open bgzip subprocess streams")

            for line in in_fh:
                if line.startswith("#"):
                    proc.stdin.write(line)
                    continue

                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 6 and cols[5].strip().lower() in {"inf", "infinity"}:
                    cols[5] = "."

                proc.stdin.write("\t".join(cols) + "\n")

            proc.stdin.close()
            stderr = proc.stderr.read()
            ret = proc.wait()

            if ret != 0:
                raise RuntimeError(f"bgzip failed while fixing QUAL: {stderr.strip()}")

        self.log(f"[INFO] Indexing fixed VCF: {fixed_vcf}")
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(fixed_vcf)])

        return fixed_vcf

    def flatten_pharmcat_outputs(self, tmp_dir):
        tmp_dir = Path(tmp_dir)

        for path in sorted(tmp_dir.rglob("*")):
            if not path.is_file():
                continue

            rel = path.relative_to(tmp_dir)
            flat_name = "__".join(rel.parts)

            if not (
                flat_name.startswith(f"{self.sample_id}.")
                or flat_name.startswith(f"{self.sample_id}_")
            ):
                flat_name = f"{self.sample_id}.{flat_name}"

            dest = self.output_dir / flat_name

            if dest.exists():
                if dest.is_file():
                    dest.unlink()
                else:
                    shutil.rmtree(dest, ignore_errors=True)

            shutil.move(str(path), str(dest))

        shutil.rmtree(tmp_dir, ignore_errors=True)

    def run_pharmcat(self, fixed_vcf):
        tmp_dir = self.output_dir / f".{self.sample_id}.pharmcat_tmp"
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir, ignore_errors=True)
        tmp_dir.mkdir(parents=True, exist_ok=False)

        self.run_cmd([
            "java",
            "-jar",
            str(PHARMCAT_JAR),
            "-vcf", str(fixed_vcf),
            "-o", str(tmp_dir),
        ])

        self.flatten_pharmcat_outputs(tmp_dir)

    def run(self):
        self.log("=" * 80)
        self.log(f"[INFO] sample_id : {self.sample_id}")
        self.log(f"[INFO] alignment : {self.alignment}")
        self.log(f"[INFO] reference : {self.reference_fasta}")
        self.log(f"[INFO] outdir    : {self.output_dir}")

        self.ensure_file(self.reference_fasta, "Reference FASTA")
        self.ensure_file(self.alignment, "Alignment")
        self.require_tools()
        self.require_bcftools_version()

        self.write_input_check()
        self.ensure_reference_indexes()
        self.ensure_alignment_index()

        raw_vcf = self.run_forcecall()
        pre_vcf = self.run_preprocessor(raw_vcf)
        fixed_vcf = self.fix_inf_qual(pre_vcf)
        self.run_pharmcat(fixed_vcf)

        self.log(f"[INFO] Finished: {self.sample_id}")


def main():
    args = parse_args()
    try:
        PharmcatPipeline(
            cram=args.cram,
            reference_fasta=args.reference,
            output_dir=args.outdir,
            sample_id=args.sample_id,
        ).run()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()