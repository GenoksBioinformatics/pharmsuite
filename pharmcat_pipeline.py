#!/usr/bin/env python3
import argparse
import gzip
import shlex
import shutil
import subprocess
from pathlib import Path


PHARMCAT_DIR = Path("/opt/pharmcat")

PHARMCAT_JAR = PHARMCAT_DIR / "pharmcat-3.2.0-all.jar"
PHARMCAT_PREPROCESSOR = PHARMCAT_DIR / "pharmcat_vcf_preprocessor"
POSITIONS_SITES_VCF = PHARMCAT_DIR / "pharmcat_positions.sites.vcf.gz"
POSITIONS_REF_VCF = PHARMCAT_DIR / "pharmcat_positions.vcf.bgz"

THREADS = "5"
JAVA_MEM = "100g"


def parse_args():
    p = argparse.ArgumentParser(
        description="Single-sample PharmCAT pipeline with flat output layout."
    )
    p.add_argument("--cram", required=True, help="Input CRAM or BAM path")
    p.add_argument("--reference", required=True, help="Reference FASTA path")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--sample-id", required=True, help="Sample ID prefix for all outputs")
    return p.parse_args()


class PharmcatPipelineSingle:
    def __init__(self, cram, reference_fasta, output_dir, sample_id):
        self.alignment = Path(cram).expanduser().resolve()
        self.reference_fasta = Path(reference_fasta).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.sample_id = sample_id.strip()

        if not self.sample_id:
            raise ValueError("sample_id cannot be empty")

        self.log_file = self.output_dir / f"{self.sample_id}.pipeline.log"

    def log(self, message, log_file=None):
        line = str(message).rstrip()
        print(line, flush=True)
        if log_file:
            with open(log_file, "a") as f:
                f.write(line + "\n")

    def run_cmd(self, cmd, log_file=None):
        cmd = [str(x) for x in cmd]
        self.log(f"[CMD] {' '.join(shlex.quote(x) for x in cmd)}", log_file)

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        for line in proc.stdout:
            self.log(line.rstrip(), log_file)

        ret = proc.wait()
        if ret != 0:
            raise RuntimeError(f"Command failed with exit code {ret}: {' '.join(cmd)}")

    def require_tools(self):
        for tool in ["gatk", "samtools", "bgzip", "tabix", "java", "bcftools"]:
            if shutil.which(tool) is None:
                raise RuntimeError(f"Required tool not found in PATH: {tool}")

        if not PHARMCAT_PREPROCESSOR.is_file():
            raise RuntimeError(f"PharmCAT preprocessor not found: {PHARMCAT_PREPROCESSOR}")

    def require_bcftools_version(self):
        res = subprocess.run(
            ["bcftools", "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        line = res.stdout.splitlines()[0].strip()
        parts = line.split()
        if len(parts) < 2:
            raise RuntimeError(f"Could not parse bcftools version: {line}")

        version = parts[1]
        nums = version.split(".")
        major = int(nums[0])
        minor = int(nums[1]) if len(nums) > 1 else 0

        if (major, minor) < (1, 18):
            raise RuntimeError(f"bcftools >= 1.18 required, found {version}")

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

    def get_ref_dict_path(self):
        return self.reference_fasta.with_suffix(".dict")

    def ensure_reference_indexes(self, log_file=None):
        fai = Path(str(self.reference_fasta) + ".fai")
        ref_dict = self.get_ref_dict_path()

        if not fai.exists():
            self.log(f"[INFO] Creating FASTA index: {fai}", log_file)
            self.run_cmd(["samtools", "faidx", str(self.reference_fasta)], log_file)

        if not ref_dict.exists():
            self.log(f"[INFO] Creating sequence dictionary: {ref_dict}", log_file)
            self.run_cmd(
                [
                    "gatk",
                    "CreateSequenceDictionary",
                    "-R",
                    str(self.reference_fasta),
                    "-O",
                    str(ref_dict),
                ],
                log_file,
            )

    def ensure_alignment_index(self, alignment, log_file=None):
        alignment = Path(alignment)

        if alignment.suffix == ".cram":
            idx1 = Path(str(alignment) + ".crai")
            idx2 = alignment.with_suffix(".crai")
            if not idx1.exists() and not idx2.exists():
                self.log(f"[INFO] Creating CRAI index: {alignment}", log_file)
                self.run_cmd(["samtools", "index", str(alignment)], log_file)

        elif alignment.suffix == ".bam":
            idx1 = Path(str(alignment) + ".bai")
            idx2 = alignment.with_suffix(".bai")
            if not idx1.exists() and not idx2.exists():
                self.log(f"[INFO] Creating BAI index: {alignment}", log_file)
                self.run_cmd(["samtools", "index", "-@", THREADS, str(alignment)], log_file)

        else:
            raise ValueError(f"Unsupported alignment type: {alignment}")

    def write_input_check(self):
        out = self.output_dir / f"{self.sample_id}.inputs.txt"
        with open(out, "w") as f:
            f.write(f"sample_id\t{self.sample_id}\n")
            f.write(f"alignment\t{self.alignment}\n")
            f.write(f"reference_fasta\t{self.reference_fasta}\n")
            f.write(f"pharmcat_jar\t{PHARMCAT_JAR}\n")
            f.write(f"pharmcat_preprocessor\t{PHARMCAT_PREPROCESSOR}\n")
            f.write(f"positions_sites_vcf\t{POSITIONS_SITES_VCF}\n")
            f.write(f"positions_ref_vcf\t{POSITIONS_REF_VCF}\n")
            f.write(f"threads\t{THREADS}\n")
            f.write(f"java_mem\t{JAVA_MEM}\n")

    def run_forcecall(self, log_file=None):
        raw_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.raw.vcf.gz"

        cmd = [
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
        ]
        self.run_cmd(cmd, log_file)

        self.log(f"[INFO] Indexing raw VCF: {raw_vcf}", log_file)
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(raw_vcf)], log_file)

        return raw_vcf

    def run_preprocessor(self, raw_vcf, log_file=None):
        cmd = [
            str(PHARMCAT_PREPROCESSOR),
            "-vcf", str(raw_vcf),
            "-o", str(self.output_dir),
            "-bf", self.sample_id,
            "-refVcf", str(POSITIONS_REF_VCF),
            "-refFna", str(self.reference_fasta),
        ]
        self.run_cmd(cmd, log_file)
        return self.find_preprocessed_vcf()

    def find_preprocessed_vcf(self):
        exact_bgz = self.output_dir / f"{self.sample_id}.preprocessed.vcf.bgz"
        exact_gz = self.output_dir / f"{self.sample_id}.preprocessed.vcf.gz"

        if exact_bgz.exists():
            return exact_bgz
        if exact_gz.exists():
            return exact_gz

        candidates = sorted(self.output_dir.glob(f"{self.sample_id}*.preprocessed.vcf.bgz"))
        if len(candidates) == 1:
            return candidates[0]

        candidates = sorted(self.output_dir.glob(f"{self.sample_id}*.preprocessed.vcf.gz"))
        if len(candidates) == 1:
            return candidates[0]

        raise FileNotFoundError(
            f"Could not uniquely find preprocessed VCF for prefix '{self.sample_id}' in {self.output_dir}"
        )

    def open_text_auto(self, path):
        path = Path(path)
        if path.suffix in [".gz", ".bgz"]:
            return gzip.open(path, "rt")
        return open(path, "r")

    def fix_inf_qual(self, preprocessed_vcf, log_file=None):
        fixed_vcf = self.output_dir / f"{self.sample_id}.preprocessed.fixed.vcf.gz"
        self.log(f"[INFO] Fixing QUAL=inf: {preprocessed_vcf} -> {fixed_vcf}", log_file)

        with self.open_text_auto(preprocessed_vcf) as in_fh, open(fixed_vcf, "wb") as out_fh:
            proc = subprocess.Popen(
                ["bgzip", "-c"],
                stdin=subprocess.PIPE,
                stdout=out_fh,
                stderr=subprocess.PIPE,
                text=True,
            )

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

        return fixed_vcf

    def index_fixed_vcf(self, fixed_vcf, log_file=None):
        self.log(f"[INFO] Indexing fixed VCF: {fixed_vcf}", log_file)
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(fixed_vcf)], log_file)

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
                raise FileExistsError(f"Refusing to overwrite existing file: {dest}")

            shutil.move(str(path), str(dest))

        shutil.rmtree(tmp_dir, ignore_errors=True)

    def run_pharmcat(self, fixed_vcf, log_file=None):
        tmp_dir = self.output_dir / f".{self.sample_id}.pharmcat_tmp"
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=False)

        cmd = [
            "java",
            "-jar",
            str(PHARMCAT_JAR),
            "-vcf", str(fixed_vcf),
            "-o", str(tmp_dir),
        ]
        self.run_cmd(cmd, log_file)
        self.flatten_pharmcat_outputs(tmp_dir)

    def run(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.require_tools()
        self.require_bcftools_version()

        self.ensure_file(self.reference_fasta, "Reference FASTA")
        self.ensure_file(self.alignment, "Alignment")
        self.ensure_file(PHARMCAT_JAR, "PharmCAT JAR")
        self.ensure_file(PHARMCAT_PREPROCESSOR, "PharmCAT preprocessor")
        self.ensure_file(POSITIONS_SITES_VCF, "Positions sites VCF")
        self.ensure_file(POSITIONS_REF_VCF, "Positions ref VCF")

        self.require_vcf_index(POSITIONS_SITES_VCF, "Positions sites VCF")
        self.require_vcf_index(POSITIONS_REF_VCF, "Positions ref VCF")

        self.log("=" * 80, self.log_file)
        self.log(f"[INFO] sample_id : {self.sample_id}", self.log_file)
        self.log(f"[INFO] alignment : {self.alignment}", self.log_file)
        self.log(f"[INFO] reference : {self.reference_fasta}", self.log_file)
        self.log(f"[INFO] outdir    : {self.output_dir}", self.log_file)

        self.write_input_check()
        self.ensure_reference_indexes(self.log_file)
        self.ensure_alignment_index(self.alignment, self.log_file)

        raw_vcf = self.run_forcecall(self.log_file)
        pre_vcf = self.run_preprocessor(raw_vcf, self.log_file)
        fixed_vcf = self.fix_inf_qual(pre_vcf, self.log_file)
        self.index_fixed_vcf(fixed_vcf, self.log_file)
        self.run_pharmcat(fixed_vcf, self.log_file)

        self.log(f"[INFO] Finished: {self.sample_id}", self.log_file)


def main():
    args = parse_args()
    PharmcatPipelineSingle(
        cram=args.cram,
        reference_fasta=args.reference,
        output_dir=args.outdir,
        sample_id=args.sample_id,
    ).run()


if __name__ == "__main__":
    main()