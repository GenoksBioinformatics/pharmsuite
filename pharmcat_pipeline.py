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

REF_GTS = {"0/0", "0|0", "0", "./.", ".|.", "."}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Single-sample PharmCAT pipeline. Generates a forcecalled PGx VCF, "
            "optionally refines forcecalled genotypes with a DRAGEN VCF, "
            "normalizes/splits the VCF, then runs PharmCAT."
        )
    )
    parser.add_argument("--cram", required=True, help="Input CRAM or BAM path")
    parser.add_argument("--reference", required=True, help="Reference FASTA path")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--sample-id", required=True, help="Sample ID prefix for all outputs")
    parser.add_argument(
        "--dragen-vcf",
        required=False,
        default=None,
        help=(
            "Optional DRAGEN VCF/VCF.GZ. If provided, PharmCAT forcecall genotypes are "
            "refined against DRAGEN: DRAGEN PASS non-ref genotypes are used; all other "
            "forcecalled sites are forced to 0/0."
        ),
    )
    return parser.parse_args()


class PharmcatPipeline:
    def __init__(self, cram, reference_fasta, output_dir, sample_id, dragen_vcf=None):
        self.alignment = Path(cram).expanduser().resolve()
        self.reference_fasta = Path(reference_fasta).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.sample_id = sample_id.strip()
        self.dragen_vcf = Path(dragen_vcf).expanduser().resolve() if dragen_vcf else None

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
                f"bcftools >= {MIN_BCFTOOLS_VERSION[0]}.{MIN_BCFTOOLS_VERSION[1]} "
                f"required, found {version_str}"
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
            self.run_cmd([
                "gatk",
                "CreateSequenceDictionary",
                "-R", str(self.reference_fasta),
                "-O", str(ref_dict),
            ])

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
            fh.write(f"dragen_vcf\t{self.dragen_vcf if self.dragen_vcf else 'NA'}\n")
            fh.write(f"pharmcat_jar\t{PHARMCAT_JAR}\n")
            fh.write(f"pharmcat_preprocessor\t{PHARMCAT_PREPROCESSOR}\n")
            fh.write(f"positions_sites_vcf\t{POSITIONS_SITES_VCF}\n")
            fh.write(f"positions_ref_vcf\t{POSITIONS_REF_VCF}\n")
            fh.write(f"threads\t{THREADS}\n")
            fh.write(f"java_mem\t{JAVA_MEM}\n")

    def open_text_auto(self, path):
        path = Path(path)
        suffixes = path.suffixes
        if path.suffix in {".gz", ".bgz"} or suffixes[-2:] in [[".vcf", ".gz"]]:
            return gzip.open(path, "rt")
        return open(path, "r")

    @staticmethod
    def normalize_gt(gt):
        if gt is None:
            return None, "/"
        sep = "|" if "|" in gt else "/"
        alleles = gt.replace("|", "/").split("/")
        return alleles, sep

    @staticmethod
    def gt_is_nonref(gt):
        if gt is None or gt in REF_GTS:
            return False
        alleles, _sep = PharmcatPipeline.normalize_gt(gt)
        if not alleles:
            return False
        return any(a not in {"0", "."} for a in alleles)

    @staticmethod
    def force_ref_gt_like_existing(_gt):
        """
        PharmCAT-oriented policy:
        Always force unsupported/non-DRAGEN PGx genotypes to diploid reference 0/0.
        Never emit ./. or haploid 0.
        """
        return "0/0"

    @staticmethod
    def variant_equivalence_key(chrom, pos, ref, alt):
        """
        Build a parsimonious allele-comparison key without rewriting the VCF.

        This does NOT split multiallelic PharmCAT records. It is only used internally
        to compare DRAGEN's often-minimal representation with PharmCAT/HaplotypeCaller
        forcecall representation.

        Examples:
          chr10:94942290 CG>TG       -> chr10:94942290 C>T
          chr2:233760233 CAT>CATAT   -> chr2:233760233 C>CAT
        """
        if alt in {".", ""} or alt.startswith("<") or ref in {".", ""}:
            return chrom, str(pos), ref, alt

        pos_int = int(pos)
        ref_s = str(ref)
        alt_s = str(alt)

        # Trim identical suffix while both alleles keep at least one base.
        while len(ref_s) > 1 and len(alt_s) > 1 and ref_s[-1] == alt_s[-1]:
            ref_s = ref_s[:-1]
            alt_s = alt_s[:-1]

        # Trim identical prefix while both alleles keep at least one base.
        while len(ref_s) > 1 and len(alt_s) > 1 and ref_s[0] == alt_s[0]:
            ref_s = ref_s[1:]
            alt_s = alt_s[1:]
            pos_int += 1

        return chrom, str(pos_int), ref_s, alt_s

    @staticmethod
    def split_info(info_str):
        info = {}
        order = []
        if info_str in {"", "."}:
            return info, order

        for item in info_str.split(";"):
            if not item:
                continue
            if "=" in item:
                key, value = item.split("=", 1)
            else:
                key, value = item, None
            info[key] = value
            order.append(key)

        return info, order

    @staticmethod
    def join_info(info, order):
        seen = set()
        parts = []

        for key in order:
            if key in seen or key not in info:
                continue
            seen.add(key)
            value = info[key]
            if value is None:
                parts.append(key)
            else:
                parts.append(f"{key}={value}")

        for key in sorted(info):
            if key in seen:
                continue
            value = info[key]
            if value is None:
                parts.append(key)
            else:
                parts.append(f"{key}={value}")

        return ";".join(parts) if parts else "."

    @staticmethod
    def format_af(value):
        if value == 0:
            return "0"
        if value == 1:
            return "1"
        return f"{value:.6g}"

    def update_info_counts_from_gt(self, info_str, gt, n_alts):
        """
        Keep AC/AF/AN and MLEAC/MLEAF consistent after genotype refinement.
        PharmCAT primarily needs GT, but stale INFO fields are confusing and risky.
        """
        info, order = self.split_info(info_str)

        alleles, _sep = self.normalize_gt(gt)
        if not alleles:
            alleles = ["0", "0"]

        called_alleles = [a for a in alleles if a != "."]
        an = len(called_alleles)
        ac = [0] * n_alts

        for allele in called_alleles:
            try:
                idx = int(allele)
            except ValueError:
                continue
            if idx > 0 and idx <= n_alts:
                ac[idx - 1] += 1

        af = [(x / an) if an else 0 for x in ac]

        info["AC"] = ",".join(str(x) for x in ac)
        info["AN"] = str(an)
        info["AF"] = ",".join(self.format_af(x) for x in af)

        if "MLEAC" in info:
            info["MLEAC"] = info["AC"]
        if "MLEAF" in info:
            info["MLEAF"] = info["AF"]

        for key in ["AC", "AF", "AN"]:
            if key not in order:
                order.append(key)

        return self.join_info(info, order)

    def parse_dragen_gt_by_equivalence_key(self, dragen_vcf):
        """
        Reads DRAGEN VCF and indexes PASS non-ref genotypes by minimal allele key.

        This avoids false 0/0 conversions when DRAGEN writes a minimal allele
        representation but the PharmCAT forcecall VCF uses a longer parsimonious
        multiallelic PGx representation.
        """
        self.ensure_file(dragen_vcf, "DRAGEN VCF")

        dragen_by_allele_key = {}
        samples = []
        total_records = 0
        kept_records = 0
        indexed_alleles = 0

        with self.open_text_auto(dragen_vcf) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue

                if line.startswith("#CHROM"):
                    header = line.rstrip("\n").split("\t")
                    samples = header[9:]
                    if not samples:
                        raise RuntimeError("DRAGEN VCF does not contain a sample column")
                    continue

                if not line.strip():
                    continue

                cols = line.rstrip("\n").split("\t")
                if len(cols) < 10:
                    continue

                total_records += 1

                chrom, pos, _id, ref, alt, _qual, filt = cols[:7]
                if filt not in {"PASS", "."}:
                    continue

                fmt_keys = cols[8].split(":")
                sample_values = cols[9].split(":")
                fmt = dict(zip(fmt_keys, sample_values))
                gt = fmt.get("GT")

                if not self.gt_is_nonref(gt):
                    continue

                alts = alt.split(",")
                allele_keys = [
                    self.variant_equivalence_key(chrom, pos, ref, alt_allele)
                    for alt_allele in alts
                ]

                record = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alts": alts,
                    "gt": gt,
                    "filter": filt,
                    "allele_keys": allele_keys,
                }

                # Index only alleles that are actually present in the DRAGEN genotype.
                dragen_alleles, _sep = self.normalize_gt(gt)
                used_alt_indexes = set()
                for allele in dragen_alleles or []:
                    if allele in {"0", "."}:
                        continue
                    try:
                        idx = int(allele) - 1
                    except ValueError:
                        continue
                    if 0 <= idx < len(allele_keys):
                        used_alt_indexes.add(idx)

                if not used_alt_indexes:
                    continue

                for idx in used_alt_indexes:
                    dragen_by_allele_key.setdefault(allele_keys[idx], []).append(record)
                    indexed_alleles += 1

                kept_records += 1

        if not samples:
            raise RuntimeError("Could not find #CHROM header line in DRAGEN VCF")

        self.log(f"[INFO] DRAGEN sample used for refinement: {samples[0]}")
        self.log(f"[INFO] DRAGEN records scanned: {total_records}")
        self.log(f"[INFO] DRAGEN PASS non-ref records loaded: {kept_records}")
        self.log(f"[INFO] DRAGEN non-ref allele keys indexed: {indexed_alleles}")

        return dragen_by_allele_key

    def translate_dragen_gt_to_forcecall_gt(self, dragen_record, forcecall_key_to_index, forcecall_ploidy):
        """
        Translate DRAGEN GT allele indexes to forcecall ALT indexes using equivalence keys.
        Returns None if the DRAGEN genotype cannot be represented in the forcecall record.
        """
        dragen_gt = dragen_record["gt"]
        dragen_alleles, sep = self.normalize_gt(dragen_gt)
        if not dragen_alleles:
            return None

        translated = []
        for allele in dragen_alleles:
            if allele == ".":
                return None

            if allele == "0":
                translated.append("0")
                continue

            try:
                dragen_alt_index = int(allele) - 1
            except ValueError:
                return None

            allele_keys = dragen_record["allele_keys"]
            if dragen_alt_index < 0 or dragen_alt_index >= len(allele_keys):
                return None

            dragen_allele_key = allele_keys[dragen_alt_index]
            forcecall_alt_index = forcecall_key_to_index.get(dragen_allele_key)
            if forcecall_alt_index is None:
                return None

            translated.append(str(forcecall_alt_index))

        if not any(a != "0" for a in translated):
            return None

        # DRAGEN can emit haploid GTs on chrX/chrY. PharmCAT input is safer as diploid.
        if len(translated) == 1 and forcecall_ploidy >= 2:
            translated = [translated[0]] * forcecall_ploidy
            sep = "/"

        return sep.join(translated)

    def refine_forcecall_with_dragen(self, forcecall_vcf, dragen_vcf):
        """
        Final PharmCAT genotype policy:
          - If DRAGEN has a PASS non-ref genotype equivalent to a PharmCAT forcecall ALT,
            use the DRAGEN genotype translated onto the forcecall ALT order.
          - If DRAGEN does not support a non-ref genotype at that PharmCAT forcecall record,
            force GT to 0/0.
          - Never emit ./. in the refined VCF.

        Important:
          We do NOT run bcftools norm -m -any here because PharmCAT expects certain
          PGx INDELs as combined multiallelic records. Equivalence matching is internal
          only and does not atomize or rewrite the PharmCAT-style VCF representation.
        """
        refined_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.dragen_refined.vcf.gz"

        self.log("[INFO] Refining forcecall VCF with DRAGEN VCF")
        self.log(f"[INFO] Forcecall VCF: {forcecall_vcf}")
        self.log(f"[INFO] DRAGEN VCF   : {dragen_vcf}")
        self.log(f"[INFO] Refined VCF  : {refined_vcf}")

        dragen_by_allele_key = self.parse_dragen_gt_by_equivalence_key(dragen_vcf)

        used_dragen_nonref = 0
        forced_to_ref = 0
        already_ref = 0
        untranslatable_dragen_sites = 0
        total_variants = 0

        with self.open_text_auto(forcecall_vcf) as in_fh, open(refined_vcf, "wb") as out_fh:
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
                if len(cols) < 10:
                    proc.stdin.write(line)
                    continue

                total_variants += 1

                chrom = cols[0]
                pos = cols[1]
                ref = cols[3]
                alt = cols[4]
                forcecall_alts = alt.split(",")

                fmt_keys = cols[8].split(":")
                sample_values = cols[9].split(":")

                if "GT" not in fmt_keys:
                    raise RuntimeError(f"FORMAT/GT not found in forcecall VCF at {chrom}:{pos}")

                gt_index = fmt_keys.index("GT")
                original_gt = sample_values[gt_index] if gt_index < len(sample_values) else None

                original_alleles, _original_sep = self.normalize_gt(original_gt)
                forcecall_ploidy = max(2, len(original_alleles or []))

                forcecall_key_to_index = {}
                candidate_records = []
                seen_record_ids = set()

                for idx, alt_allele in enumerate(forcecall_alts, start=1):
                    key = self.variant_equivalence_key(chrom, pos, ref, alt_allele)
                    forcecall_key_to_index[key] = idx

                    for rec in dragen_by_allele_key.get(key, []):
                        rec_id = id(rec)
                        if rec_id not in seen_record_ids:
                            seen_record_ids.add(rec_id)
                            candidate_records.append(rec)

                translated_gt = None
                for dragen_record in candidate_records:
                    translated_gt = self.translate_dragen_gt_to_forcecall_gt(
                        dragen_record=dragen_record,
                        forcecall_key_to_index=forcecall_key_to_index,
                        forcecall_ploidy=forcecall_ploidy,
                    )
                    if translated_gt is not None:
                        break

                if translated_gt is not None:
                    sample_values[gt_index] = translated_gt
                    used_dragen_nonref += 1
                else:
                    ref_gt = self.force_ref_gt_like_existing(original_gt)
                    if original_gt == ref_gt:
                        already_ref += 1
                    else:
                        forced_to_ref += 1

                    if candidate_records:
                        untranslatable_dragen_sites += 1

                    sample_values[gt_index] = ref_gt

                cols[7] = self.update_info_counts_from_gt(
                    info_str=cols[7],
                    gt=sample_values[gt_index],
                    n_alts=len(forcecall_alts),
                )
                cols[9] = ":".join(sample_values)
                proc.stdin.write("\t".join(cols) + "\n")

            proc.stdin.close()
            stderr = proc.stderr.read()
            ret = proc.wait()

            if ret != 0:
                raise RuntimeError(f"bgzip failed while refining VCF: {stderr.strip()}")

        self.log(f"[INFO] Indexing DRAGEN-refined VCF: {refined_vcf}")
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(refined_vcf)])

        self.log("[INFO] DRAGEN-refined VCF summary:")
        self.log(f"[INFO]   forcecall records scanned                : {total_variants}")
        self.log(f"[INFO]   DRAGEN non-ref genotypes applied         : {used_dragen_nonref}")
        self.log(f"[INFO]   forcecall genotypes forced to 0/0        : {forced_to_ref}")
        self.log(f"[INFO]   forcecall genotypes already 0/0          : {already_ref}")
        self.log(f"[INFO]   candidate DRAGEN records not translatable: {untranslatable_dragen_sites}")

        return refined_vcf

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
        self.log(f"[INFO] dragen_vcf: {self.dragen_vcf if self.dragen_vcf else 'NA'}")

        self.ensure_file(self.reference_fasta, "Reference FASTA")
        self.ensure_file(self.alignment, "Alignment")
        if self.dragen_vcf:
            self.ensure_file(self.dragen_vcf, "DRAGEN VCF")

        self.require_tools()
        self.require_bcftools_version()

        self.write_input_check()
        self.ensure_reference_indexes()
        self.ensure_alignment_index()

        raw_vcf = self.run_forcecall()

        if self.dragen_vcf:
            pharmcat_input_vcf = self.refine_forcecall_with_dragen(
                forcecall_vcf=raw_vcf,
                dragen_vcf=self.dragen_vcf,
            )
        else:
            pharmcat_input_vcf = raw_vcf

        pre_vcf = self.run_preprocessor(pharmcat_input_vcf)
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
            dragen_vcf=args.dragen_vcf,
        ).run()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()