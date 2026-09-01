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

# The only true multi-base (MNP) sites in the PharmCAT panel: both are RYR1
# dinucleotide substitutions and require GATK's --max-mnp-distance 1 to be
# called as a single 2bp event instead of two separate 1bp SNPs.
# (chrom, anchor_pos, ref_len_in_bp)
MNP_SPLICE_SITES = [
    ("chr19", 38572266, 2),  # rs193922862, TC>CT
    ("chr19", 38580039, 2),  # unnamed, TT>AA
]

# The only sex-chromosome gene in the PharmCAT panel (G6PD) is on chrX,
# outside the pseudoautosomal regions, so it is hemizygous (ploidy 1) in XY
# samples. All chrX positions in the panel belong to this gene.
CHRX_HAPLOID_CHROM = "chrX"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Single-sample PharmCAT pipeline. Generates a PharmCAT PGx-site forcecall VCF "
            "with GATK HaplotypeCaller, optionally uses a DRAGEN VCF only as a presence/absence "
            "support filter for forcecalled non-ref genotypes, then runs PharmCAT."
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
            "Optional DRAGEN VCF/VCF.GZ. Used only to check whether DRAGEN has any PASS "
            "non-ref call at a forcecalled position. If forcecall GT is non-ref and DRAGEN "
            "has no PASS non-ref call at that CHROM:POS, only GT is changed to 0/0. "
            "INFO/AD/DP/GQ/PL are not modified."
        ),
    )
    parser.add_argument(
        "--ploidy-metrics",
        required=False,
        default=None,
        help=(
            "Optional DRAGEN ploidy_estimation_metrics.csv for this sample. If present and "
            "the 'Ploidy estimation' row reads XY, chrX (G6PD) positions are forcecalled "
            "with --sample-ploidy 1 instead of the default diploid model. Any other value "
            "(XX, etc.) or a missing file leaves chrX diploid, unchanged."
        ),
    )
    return parser.parse_args()


class PharmcatPipeline:
    def __init__(self, cram, reference_fasta, output_dir, sample_id, dragen_vcf=None, ploidy_metrics=None):
        self.alignment = Path(cram).expanduser().resolve()
        self.reference_fasta = Path(reference_fasta).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.sample_id = sample_id.strip()
        self.dragen_vcf = Path(dragen_vcf).expanduser().resolve() if dragen_vcf else None
        self.ploidy_metrics = Path(ploidy_metrics).expanduser().resolve() if ploidy_metrics else None
        self.sample_sex = None

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

    @staticmethod
    def ensure_file(path, label):
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")

    @staticmethod
    def require_vcf_index(vcf_path, label):
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

    @staticmethod
    def require_bcftools_version():
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

    def detect_sample_sex(self):
        """
        Read the 'Ploidy estimation' row out of a DRAGEN ploidy_estimation_metrics.csv.

        Returns "XY", "XX", some other DRAGEN-reported value (e.g. "XO", "XXY"), or
        None if --ploidy-metrics wasn't given. Only an exact "XY" triggers haploid
        chrX forcecalling; everything else keeps the default diploid behavior.
        """
        if self.ploidy_metrics is None:
            self.log("[INFO] No --ploidy-metrics given: chrX (G6PD) will be forcecalled as diploid")
            return None

        self.ensure_file(self.ploidy_metrics, "Ploidy metrics CSV")

        ploidy_value = None
        with open(self.ploidy_metrics, "r") as fh:
            for line in fh:
                fields = line.rstrip("\n").split(",")
                if len(fields) >= 4 and fields[2].strip() == "Ploidy estimation":
                    ploidy_value = fields[3].strip()
                    break

        if ploidy_value is None:
            raise RuntimeError(
                f"Could not find a 'Ploidy estimation' row in {self.ploidy_metrics}"
            )

        self.log(f"[INFO] Ploidy metrics report: Ploidy estimation = {ploidy_value}")
        return ploidy_value

    def write_input_check(self):
        out = self.output_dir / f"{self.sample_id}.inputs.txt"
        with open(out, "w") as fh:
            fh.write(f"sample_id\t{self.sample_id}\n")
            fh.write(f"alignment\t{self.alignment}\n")
            fh.write(f"reference_fasta\t{self.reference_fasta}\n")
            fh.write(f"dragen_vcf\t{self.dragen_vcf if self.dragen_vcf else 'NA'}\n")
            mnp_sites_str = ";".join(f"{c}:{p}" for c, p, _ in MNP_SPLICE_SITES)
            fh.write("mnp_calling_policy\tdefault_forcecall_genomewide_plus_mnp1_spliced_at_known_dinucleotide_sites\n")
            fh.write(f"mnp_spliced_sites\t{mnp_sites_str}\n")
            fh.write(f"ploidy_metrics\t{self.ploidy_metrics if self.ploidy_metrics else 'NA'}\n")
            fh.write(f"detected_sample_sex\t{self.sample_sex if self.sample_sex else 'NA'}\n")
            fh.write(f"chrx_ploidy_policy\thaploid_sample_ploidy_1_if_XY_else_diploid_default\t(chrom={CHRX_HAPLOID_CHROM}, gene=G6PD)\n")
            fh.write("dragen_mode\tposition_support_filter_only\n")
            fh.write("dragen_policy\tforcecall_nonref_without_dragen_position_support_to_0/0\n")
            fh.write("ad_balance_policy\t0/x_ALT_fraction_ge_0.85_to_x/x__0/x_ALT_fraction_le_0.15_to_0/0\n")
            fh.write("ad_balance_fraction\tAD[x]/(AD[0]+AD[x])\n")
            fh.write("fields_modified_by_dragen_filter\tGT_only\n")
            fh.write(f"pharmcat_jar\t{PHARMCAT_JAR}\n")
            fh.write(f"pharmcat_preprocessor\t{PHARMCAT_PREPROCESSOR}\n")
            fh.write(f"positions_sites_vcf\t{POSITIONS_SITES_VCF}\n")
            fh.write(f"positions_ref_vcf\t{POSITIONS_REF_VCF}\n")
            fh.write(f"threads\t{THREADS}\n")
            fh.write(f"java_mem\t{JAVA_MEM}\n")

    @staticmethod
    def open_text_auto(path):
        path = Path(path)
        suffixes = path.suffixes
        if path.suffix in {".gz", ".bgz"} or suffixes[-2:] == [".vcf", ".gz"]:
            return gzip.open(path, "rt")
        return open(path, "r")

    @staticmethod
    def normalize_chrom_keys(chrom, pos):
        """
        Make chr/no-chr naming differences less fragile.
        We still write the original forcecall chromosome unchanged.
        """
        chrom = str(chrom)
        pos = str(pos)
        keys = {(chrom, pos)}

        if chrom.startswith("chr"):
            keys.add((chrom[3:], pos))
        else:
            keys.add((f"chr{chrom}", pos))

        return keys

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
    def force_ref_gt(_original_gt):
        """
        Requested policy:
        If forcecall is non-ref but DRAGEN has no PASS non-ref call at that position,
        modify only GT and set it to diploid 0/0.
        """
        return "0/0"

    @staticmethod
    def split_gt_alleles(gt):
        """
        Split GT like 0/1, 0|1, 1/1.

        Returns:
            tuple[list[int] | None, str]: allele list and separator.
        """
        if gt is None:
            return None, "/"

        gt = str(gt).strip()
        if gt in {"", ".", "./.", ".|."}:
            return None, "/"

        sep = "|" if "|" in gt else "/"
        parts = gt.split(sep)

        if not parts or any(part in {"", "."} for part in parts):
            return None, sep

        try:
            alleles = [int(part) for part in parts]
        except ValueError:
            return None, sep

        return alleles, sep

    @staticmethod
    def force_hom_alt_gt(original_gt, alt_index):
        """
        Convert GT to x/x or x|x while preserving ploidy and separator.
        Example: 0/2 -> 2/2, 0|2 -> 2|2.
        """
        alleles, sep = PharmcatPipeline.split_gt_alleles(original_gt)
        if alleles is None:
            return original_gt

        return sep.join([str(alt_index)] * len(alleles))

    @staticmethod
    def get_ad_based_gt_override(original_gt, fmt_keys, sample_values, high_ab=0.85, low_ab=0.15):
        """
        Apply AD-based genotype correction only for simple diploid 0/x calls.

        Rule:
            alt_fraction = AD[x] / (AD[0] + AD[x])

            alt_fraction >= high_ab -> x/x
            alt_fraction <= low_ab  -> 0/0

        Returns:
            tuple[str | None, str | None]: new GT and reason.
        """
        if original_gt is None or "AD" not in fmt_keys:
            return None, None

        alleles, _sep = PharmcatPipeline.split_gt_alleles(original_gt)
        if alleles is None:
            return None, None

        # Only diploid 0/x genotypes. Leave 1/2, 1/1, 2/2, haploid etc. unchanged.
        if len(alleles) != 2:
            return None, None

        unique_alleles = set(alleles)
        if 0 not in unique_alleles:
            return None, None

        alt_alleles = [allele for allele in unique_alleles if allele != 0]
        if len(alt_alleles) != 1:
            return None, None

        alt_index = alt_alleles[0]
        ad_index = fmt_keys.index("AD")
        if ad_index >= len(sample_values):
            return None, None

        ad_raw = sample_values[ad_index]
        if ad_raw in {"", "."}:
            return None, None

        try:
            ad_values = [int(value) if value not in {"", "."} else 0 for value in ad_raw.split(",")]
        except ValueError:
            return None, None

        if len(ad_values) <= alt_index:
            return None, None

        ref_ad = ad_values[0]
        alt_ad = ad_values[alt_index]
        denom = ref_ad + alt_ad
        if denom <= 0:
            return None, None

        alt_fraction = alt_ad / denom

        if alt_fraction >= high_ab:
            return PharmcatPipeline.force_hom_alt_gt(original_gt, alt_index), "AD_HIGH_TO_HOM_ALT"

        if alt_fraction <= low_ab:
            return PharmcatPipeline.force_ref_gt(original_gt), "AD_LOW_TO_REF"

        return None, None

    def load_dragen_supported_positions(self, dragen_vcf):
        """
        Return set of CHROM:POS positions where DRAGEN has a PASS non-ref genotype.

        Important:
        - No allele matching.
        - No GT translation.
        - No INFO/FORMAT usage except GT and FILTER.
        - Multiallelic/indel/SNV are treated the same: if DRAGEN called any PASS
          non-ref at the same position, forcecall record is left untouched.
        """
        self.ensure_file(dragen_vcf, "DRAGEN VCF")

        supported_positions = set()
        samples = []
        total_records = 0
        pass_nonref_records = 0

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

                chrom, pos = cols[0], cols[1]
                filt = cols[6]

                if filt not in {"PASS", "."}:
                    continue

                fmt_keys = cols[8].split(":")
                sample_values = cols[9].split(":")
                fmt = dict(zip(fmt_keys, sample_values))
                gt = fmt.get("GT")

                if not self.gt_is_nonref(gt):
                    continue

                for key in self.normalize_chrom_keys(chrom, pos):
                    supported_positions.add(key)

                pass_nonref_records += 1

        if not samples:
            raise RuntimeError("Could not find #CHROM header line in DRAGEN VCF")

        self.log(f"[INFO] DRAGEN sample used for support check: {samples[0]}")
        self.log(f"[INFO] DRAGEN records scanned: {total_records}")
        self.log(f"[INFO] DRAGEN PASS non-ref records loaded: {pass_nonref_records}")
        self.log(f"[INFO] DRAGEN supported CHROM:POS keys stored: {len(supported_positions)}")

        return supported_positions

    def filter_forcecall_nonref_by_dragen_position(self, forcecall_vcf, dragen_vcf):
        """
        Core policy:
          - If forcecall GT is reference, leave record unchanged.
          - If forcecall GT is 0/x and AD[x] / (AD[0] + AD[x]) >= 0.85,
            change only GT to x/x.
          - If forcecall GT is 0/x and AD[x] / (AD[0] + AD[x]) <= 0.15,
            change only GT to 0/0.
          - Otherwise, if forcecall GT is non-ref and DRAGEN has any PASS non-ref call
            at same CHROM:POS, leave record unchanged.
          - Otherwise, if forcecall GT is non-ref and DRAGEN has no PASS non-ref call
            at same CHROM:POS, change only GT to 0/0.
          - Do not modify INFO, FILTER, QUAL, AD, DP, GQ, PL or any other FORMAT value.
        """
        filtered_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.dragen_position_filtered.vcf.gz"

        self.log("[INFO] Filtering forcecall non-ref genotypes using AD balance and DRAGEN position support")
        self.log(f"[INFO] Forcecall VCF: {forcecall_vcf}")
        self.log(f"[INFO] DRAGEN VCF   : {dragen_vcf}")
        self.log(f"[INFO] Filtered VCF : {filtered_vcf}")

        supported_positions = self.load_dragen_supported_positions(dragen_vcf)

        total_records = 0
        forcecall_ref_unchanged = 0
        forcecall_nonref_ad85_to_hom_alt = 0
        forcecall_nonref_ad15_to_ref = 0
        forcecall_nonref_ad_no_decision = 0
        forcecall_nonref_supported_unchanged = 0
        forcecall_nonref_forced_to_ref = 0

        with self.open_text_auto(forcecall_vcf) as in_fh, open(filtered_vcf, "wb") as out_fh:
            proc = subprocess.Popen(
                ["bgzip", "-c"],
                stdin=subprocess.PIPE,
                stdout=out_fh,
                stderr=subprocess.PIPE,
                text=True,
            )

            if proc.stdin is None or proc.stderr is None:
                raise RuntimeError("Failed to open bgzip subprocess streams")

            try:
                for line in in_fh:
                    if line.startswith("#"):
                        proc.stdin.write(line)
                        continue

                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 10:
                        proc.stdin.write(line)
                        continue

                    total_records += 1

                    chrom, pos = cols[0], cols[1]
                    fmt_keys = cols[8].split(":")
                    sample_values = cols[9].split(":")

                    if "GT" not in fmt_keys:
                        raise RuntimeError(f"FORMAT/GT not found in forcecall VCF at {chrom}:{pos}")

                    gt_index = fmt_keys.index("GT")
                    if gt_index >= len(sample_values):
                        raise RuntimeError(f"Sample FORMAT value missing GT at {chrom}:{pos}")

                    original_gt = sample_values[gt_index]

                    if not self.gt_is_nonref(original_gt):
                        forcecall_ref_unchanged += 1
                        proc.stdin.write("\t".join(cols) + "\n")
                        continue

                    ad_based_gt, ad_reason = self.get_ad_based_gt_override(
                        original_gt=original_gt,
                        fmt_keys=fmt_keys,
                        sample_values=sample_values,
                        high_ab=0.85,
                        low_ab=0.15,
                    )

                    if ad_based_gt is not None and ad_based_gt != original_gt:
                        sample_values[gt_index] = ad_based_gt
                        cols[9] = ":".join(sample_values)

                        if ad_reason == "AD_HIGH_TO_HOM_ALT":
                            forcecall_nonref_ad85_to_hom_alt += 1
                        elif ad_reason == "AD_LOW_TO_REF":
                            forcecall_nonref_ad15_to_ref += 1

                        proc.stdin.write("\t".join(cols) + "\n")
                        continue

                    forcecall_nonref_ad_no_decision += 1

                    has_dragen_support = any(
                        key in supported_positions
                        for key in self.normalize_chrom_keys(chrom, pos)
                    )

                    if has_dragen_support:
                        forcecall_nonref_supported_unchanged += 1
                    else:
                        sample_values[gt_index] = self.force_ref_gt(original_gt)
                        cols[9] = ":".join(sample_values)
                        forcecall_nonref_forced_to_ref += 1

                    proc.stdin.write("\t".join(cols) + "\n")

            except Exception:
                try:
                    proc.stdin.close()
                except Exception:
                    pass
                try:
                    proc.kill()
                except Exception:
                    pass
                raise

            proc.stdin.close()
            stderr = proc.stderr.read()
            ret = proc.wait()

            if ret != 0:
                raise RuntimeError(f"bgzip failed while filtering forcecall VCF: {stderr.strip()}")

        self.log(f"[INFO] Indexing DRAGEN-position-filtered VCF: {filtered_vcf}")
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(filtered_vcf)])

        self.log("[INFO] DRAGEN position support filtering summary:")
        self.log(f"[INFO]   forcecall records scanned                         : {total_records}")
        self.log(f"[INFO]   forcecall reference genotypes left unchanged      : {forcecall_ref_unchanged}")
        self.log(f"[INFO]   forcecall 0/x with ALT fraction >= 0.85 -> x/x   : {forcecall_nonref_ad85_to_hom_alt}")
        self.log(f"[INFO]   forcecall 0/x with ALT fraction <= 0.15 -> 0/0   : {forcecall_nonref_ad15_to_ref}")
        self.log(f"[INFO]   forcecall non-ref AD no decision                  : {forcecall_nonref_ad_no_decision}")
        self.log(f"[INFO]   forcecall non-ref with DRAGEN support unchanged   : {forcecall_nonref_supported_unchanged}")
        self.log(f"[INFO]   forcecall non-ref without DRAGEN support -> 0/0   : {forcecall_nonref_forced_to_ref}")

        return filtered_vcf

    def _run_haplotypecaller(self, output_vcf, extra_args=None):
        cmd = [
            "gatk",
            "--java-options", f"-Xmx{JAVA_MEM}",
            "HaplotypeCaller",
            "--native-pair-hmm-threads", THREADS,
            "--alleles", str(POSITIONS_SITES_VCF),
            "-R", str(self.reference_fasta),
            "-I", str(self.alignment),
            "-O", str(output_vcf),
            "-L", str(POSITIONS_SITES_VCF),
            "-ip", "20",
            "--output-mode", "EMIT_ALL_CONFIDENT_SITES",
        ]
        if extra_args:
            cmd.extend(extra_args)

        self.run_cmd(cmd)

        self.log(f"[INFO] Indexing VCF: {output_vcf}")
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(output_vcf)])

    @staticmethod
    def _mnp_removal_positions():
        """
        All (chrom, pos) keys spanned by each MNP_SPLICE_SITES anchor, e.g. a
        2bp site at pos N spans N and N+1.
        """
        positions = set()
        for chrom, anchor, ref_len in MNP_SPLICE_SITES:
            for offset in range(ref_len):
                positions.add((chrom, str(anchor + offset)))
        return positions

    def _extract_mnp_anchor_lines(self, mnp_vcf):
        anchors = {(chrom, str(pos)) for chrom, pos, _ in MNP_SPLICE_SITES}
        lines = {}

        with self.open_text_auto(mnp_vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                key = (cols[0], cols[1])
                if key in anchors:
                    lines[key] = line if line.endswith("\n") else line + "\n"

        missing = anchors - set(lines.keys())
        if missing:
            raise RuntimeError(
                f"MNP-mode forcecall is missing expected record(s) at: {sorted(missing)}"
            )

        return lines

    def splice_mnp_target_sites(self, default_vcf, mnp_vcf):
        """
        Use the default (non-MNP-merged) forcecall for the whole panel, except
        at the known RYR1 dinucleotide sites, where the MNP-merged (2bp)
        record is spliced in instead. This avoids GATK's --max-mnp-distance
        from silently entangling unrelated nearby SNPs elsewhere in the panel
        while still correctly calling the two true MNP sites.
        """
        spliced_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.raw.vcf.gz"

        self.log("[INFO] Splicing MNP-mode calls into default forcecall for known dinucleotide sites")
        self.log(f"[INFO] Default forcecall VCF : {default_vcf}")
        self.log(f"[INFO] MNP-mode forcecall VCF: {mnp_vcf}")
        self.log(f"[INFO] Spliced target sites  : {[(c, p) for c, p, _ in MNP_SPLICE_SITES]}")

        remove_positions = self._mnp_removal_positions()
        anchor_lines = self._extract_mnp_anchor_lines(mnp_vcf)
        anchor_positions_pending = {(c, str(p)) for c, p, _ in MNP_SPLICE_SITES}

        removed_count = 0
        spliced_count = 0

        with self.open_text_auto(default_vcf) as in_fh, open(spliced_vcf, "wb") as out_fh:
            proc = subprocess.Popen(
                ["bgzip", "-c"],
                stdin=subprocess.PIPE,
                stdout=out_fh,
                stderr=subprocess.PIPE,
                text=True,
            )

            if proc.stdin is None or proc.stderr is None:
                raise RuntimeError("Failed to open bgzip subprocess streams")

            try:
                for line in in_fh:
                    if line.startswith("#"):
                        proc.stdin.write(line)
                        continue

                    cols = line.rstrip("\n").split("\t")
                    key = (cols[0], cols[1])

                    if key in remove_positions:
                        removed_count += 1
                        if key in anchor_positions_pending:
                            proc.stdin.write(anchor_lines[key])
                            anchor_positions_pending.discard(key)
                            spliced_count += 1
                        continue

                    proc.stdin.write(line)
            except Exception:
                try:
                    proc.stdin.close()
                except Exception:
                    pass
                try:
                    proc.kill()
                except Exception:
                    pass
                raise

            proc.stdin.close()
            stderr = proc.stderr.read()
            ret = proc.wait()

            if ret != 0:
                raise RuntimeError(f"bgzip failed while splicing MNP sites: {stderr.strip()}")

        if anchor_positions_pending:
            raise RuntimeError(
                f"MNP anchor position(s) not found in default forcecall VCF: {sorted(anchor_positions_pending)}"
            )

        self.log(f"[INFO] Removed {removed_count} default-run record(s) within MNP target site span(s)")
        self.log(f"[INFO] Inserted {spliced_count} MNP-mode record(s) at target anchor position(s)")

        self.log(f"[INFO] Indexing spliced VCF: {spliced_vcf}")
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(spliced_vcf)])

        return spliced_vcf

    def splice_chrx_haploid_calls(self, base_vcf, chrx_vcf):
        """
        Replace all chrX records in base_vcf (called diploid by default) with the
        --sample-ploidy 1 records from chrx_vcf. Only used for XY samples. chrX is
        the last chromosome block in the PharmCAT panel, so the haploid records can
        simply be appended after every non-chrX record has been streamed through.
        """
        adjusted_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.chrx_adjusted.vcf.gz"

        self.log("[INFO] Splicing haploid (--sample-ploidy 1) chrX calls into forcecall (sample detected as XY)")
        self.log(f"[INFO] Base forcecall VCF: {base_vcf}")
        self.log(f"[INFO] chrX haploid VCF  : {chrx_vcf}")

        chrx_lines = []
        with self.open_text_auto(chrx_vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                chrom = line.split("\t", 1)[0]
                if chrom == CHRX_HAPLOID_CHROM:
                    chrx_lines.append(line if line.endswith("\n") else line + "\n")

        if not chrx_lines:
            raise RuntimeError(f"No {CHRX_HAPLOID_CHROM} records found in haploid forcecall VCF: {chrx_vcf}")

        removed_count = 0

        with self.open_text_auto(base_vcf) as in_fh, open(adjusted_vcf, "wb") as out_fh:
            proc = subprocess.Popen(
                ["bgzip", "-c"],
                stdin=subprocess.PIPE,
                stdout=out_fh,
                stderr=subprocess.PIPE,
                text=True,
            )

            if proc.stdin is None or proc.stderr is None:
                raise RuntimeError("Failed to open bgzip subprocess streams")

            try:
                for line in in_fh:
                    if line.startswith("#"):
                        proc.stdin.write(line)
                        continue

                    chrom = line.split("\t", 1)[0]
                    if chrom == CHRX_HAPLOID_CHROM:
                        removed_count += 1
                        continue

                    proc.stdin.write(line)

                for line in chrx_lines:
                    proc.stdin.write(line)
            except Exception:
                try:
                    proc.stdin.close()
                except Exception:
                    pass
                try:
                    proc.kill()
                except Exception:
                    pass
                raise

            proc.stdin.close()
            stderr = proc.stderr.read()
            ret = proc.wait()

            if ret != 0:
                raise RuntimeError(f"bgzip failed while splicing chrX haploid calls: {stderr.strip()}")

        self.log(f"[INFO] Removed {removed_count} diploid {CHRX_HAPLOID_CHROM} record(s) from default forcecall")
        self.log(f"[INFO] Inserted {len(chrx_lines)} haploid {CHRX_HAPLOID_CHROM} record(s)")

        self.log(f"[INFO] Indexing chrX-adjusted VCF: {adjusted_vcf}")
        self.run_cmd(["tabix", "-f", "-p", "vcf", str(adjusted_vcf)])

        return adjusted_vcf

    def run_forcecall(self):
        default_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.default.vcf.gz"
        mnp_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.mnp_sites.vcf.gz"

        self.log("[INFO] Running default (non-MNP-merged) forcecall for the full PharmCAT panel")
        self._run_haplotypecaller(default_vcf)

        self.log("[INFO] Running MNP-merged forcecall to capture known dinucleotide PharmCAT sites (RYR1)")
        self._run_haplotypecaller(mnp_vcf, extra_args=["--max-mnp-distance", "1"])

        working_vcf = default_vcf
        sample_sex = self.sample_sex

        if sample_sex == "XY":
            chrx_vcf = self.output_dir / f"{self.sample_id}.pharmcat.forcecall.chrx_haploid.vcf.gz"
            self.log(f"[INFO] Sample detected as XY: forcecalling {CHRX_HAPLOID_CHROM} (G6PD) at --sample-ploidy 1")
            self._run_haplotypecaller(
                chrx_vcf,
                extra_args=[
                    "-L", CHRX_HAPLOID_CHROM,
                    "--interval-set-rule", "INTERSECTION",
                    "--sample-ploidy", "1",
                ],
            )
            working_vcf = self.splice_chrx_haploid_calls(working_vcf, chrx_vcf)
        elif sample_sex is not None:
            self.log(f"[INFO] Sample ploidy estimation is '{sample_sex}' (not XY): {CHRX_HAPLOID_CHROM} (G6PD) stays diploid")

        return self.splice_mnp_target_sites(working_vcf, mnp_vcf)

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
            "-reporterHtml",
            "-reporterJson",
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
            self.require_vcf_index(self.dragen_vcf, "DRAGEN VCF")

        if self.ploidy_metrics:
            self.ensure_file(self.ploidy_metrics, "Ploidy metrics CSV")

        self.require_tools()
        self.require_bcftools_version()

        self.sample_sex = self.detect_sample_sex()

        self.write_input_check()
        self.ensure_reference_indexes()
        self.ensure_alignment_index()

        raw_vcf = self.run_forcecall()

        if self.dragen_vcf:
            pharmcat_input_vcf = self.filter_forcecall_nonref_by_dragen_position(
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
            ploidy_metrics=args.ploidy_metrics,
        ).run()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()