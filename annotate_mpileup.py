#!/usr/bin/env python3
import argparse
import csv
import subprocess
from pathlib import Path

# Fixed input TSV path inside the Docker image
INPUT_TSV = Path("/opt/mpileup/pharmvariants_wfeatures.tsv")


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate loci TSV with mpileup strand/depth summary for a single CRAM."
    )
    p.add_argument(
        "--cram",
        required=True,
        help="Input CRAM path",
    )
    p.add_argument(
        "--sample-id",
        required=True,
        help="Sample ID to use as the output file name prefix",
    )
    p.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA for samtools mpileup",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    return p.parse_args()


def clean_pileup_bases(bases: str) -> str:
    """
    Remove mpileup control tokens so current-position base counting is reliable.
    Handles:
      ^ + following mapq char
      $
      +/-<num><seq>
      * # < >
    """
    cleaned = []
    i = 0
    n = len(bases)

    while i < n:
        ch = bases[i]

        if ch == "^":
            i += 2
            continue

        if ch == "$":
            i += 1
            continue

        if ch in "+-":
            i += 1
            num = []
            while i < n and bases[i].isdigit():
                num.append(bases[i])
                i += 1
            if not num:
                continue
            indel_len = int("".join(num))
            i += indel_len
            continue

        if ch in "*#<>":
            i += 1
            continue

        cleaned.append(ch)
        i += 1

    return "".join(cleaned)


def count_ref_alt_from_bases(bases: str, ref: str, alt: str):
    """
    mpileup bases meanings:
      .  = ref on forward
      ,  = ref on reverse
      A/C/G/T/N = alt/non-ref on forward
      a/c/g/t/n = alt/non-ref on reverse
    """
    ref = ref.upper()
    alt = alt.upper()

    ref_f = 0
    ref_r = 0
    alt_f = 0
    alt_r = 0

    cleaned = clean_pileup_bases(bases)

    for ch in cleaned:
        if ch == ".":
            ref_f += 1
        elif ch == ",":
            ref_r += 1
        elif ch.upper() == alt:
            if ch.isupper():
                alt_f += 1
            else:
                alt_r += 1

    return ref_f, ref_r, alt_f, alt_r


def run_mpileup(reference: str, cram: str, chrom: str, pos: str):
    region = f"{chrom}:{pos}-{pos}"
    cmd = [
        "samtools",
        "mpileup",
        "--reference",
        reference,
        "-r",
        region,
        cram,
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)

    if res.returncode != 0:
        raise RuntimeError(
            f"samtools mpileup failed for {cram} at {region}\nSTDERR:\n{res.stderr}"
        )

    stdout = res.stdout.strip()
    if not stdout:
        return None

    lines = [x for x in stdout.splitlines() if x.strip()]
    if not lines:
        return None

    return lines[-1]


def summarize_locus(reference: str, cram: str, chrom: str, pos: str, refalt: str):
    line = run_mpileup(reference, cram, chrom, pos)
    if line is None:
        return "NOCALL"

    parts = line.split("\t")
    if len(parts) < 5:
        return "NOCALL"

    try:
        depth = int(parts[3])
    except ValueError:
        return "NOCALL"

    if depth == 0:
        return "NOCALL"

    bases = parts[4] if len(parts) > 4 else ""

    try:
        ref, alt = refalt.split(">")
        ref = ref.strip().upper()
        alt = alt.strip().upper()
    except ValueError:
        return f"DP={depth}"

    if len(ref) != 1 or len(alt) != 1:
        return f"DP={depth}"

    ref_f, ref_r, alt_f, alt_r = count_ref_alt_from_bases(bases, ref, alt)
    alt_dp = alt_f + alt_r
    vaf = alt_dp / depth if depth > 0 else None
    vaf_str = f"{vaf:.4f}" if vaf is not None else "NA"

    return (
        f"DP={depth}"
        f"|REF_F={ref_f}"
        f"|REF_R={ref_r}"
        f"|ALT_F={alt_f}"
        f"|ALT_R={alt_r}"
        f"|ALT_DP={alt_dp}"
        f"|VAF={vaf_str}"
    )


def main():
    args = parse_args()

    if not INPUT_TSV.exists():
        raise SystemExit(f"ERROR: Fixed input TSV not found in container: {INPUT_TSV}")

    cram = Path(args.cram)
    if not cram.exists():
        raise SystemExit(f"ERROR: CRAM not found: {cram}")

    reference = Path(args.reference)
    if not reference.exists():
        raise SystemExit(f"ERROR: Reference FASTA not found: {reference}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    required_cols = {"chr", "location", "refalt"}
    rows = []

    with open(INPUT_TSV, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit("ERROR: Could not read input TSV header.")

        missing = required_cols - set(reader.fieldnames)
        if missing:
            raise SystemExit(
                f"ERROR: Required column(s) missing in input TSV: {', '.join(sorted(missing))}"
            )

        fieldnames = reader.fieldnames[:]

        for row in reader:
            rows.append(row)

    if not rows:
        raise SystemExit("ERROR: Input TSV appears to be empty.")

    total = len(rows)
    print(f"[START] sample_id={args.sample_id} | loci={total}", flush=True)

    for i, row in enumerate(rows, start=1):
        chrom = row["chr"].strip()
        pos = row["location"].strip()
        refalt = row["refalt"].strip()

        row[args.sample_id] = summarize_locus(
            str(reference),
            str(cram),
            chrom,
            pos,
            refalt,
        )

        if i % 50 == 0 or i == total:
            print(f"[{args.sample_id}] {i}/{total}", flush=True)

    final_cols = fieldnames + [args.sample_id]
    outfile = outdir / f"{args.sample_id}.mpileup.tsv"

    with open(outfile, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=final_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[DONE] Output: {outfile}", flush=True)


if __name__ == "__main__":
    main()