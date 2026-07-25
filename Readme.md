# PharmSuite Docker Pipeline

A lightweight Docker-based pipeline for single-sample PharmCAT preprocessing, PharmCAT reporting, DRAGEN-assisted genotype refinement, and targeted mpileup annotation.

## What it does

This container runs the pipeline for **one sample at a time** using:

- a single `CRAM` file
- a reference genome FASTA
- a user-defined `sample_id`
- a single output directory
- optionally, a DRAGEN VCF for further genotype refinement
- optionally, a DRAGEN ploidy estimation metrics CSV for correct chrX (G6PD) hemizygous calling

The image contains the required PharmCAT resources internally, so they do not need to be passed at runtime.

## Included steps

The pipeline performs:

1. **Force-calling PharmCAT positions** from the input CRAM using GATK HaplotypeCaller
   - Two HaplotypeCaller passes are run and spliced together: a default (non-MNP-merged)
     pass for the whole panel, plus a `--max-mnp-distance 1` pass used only for the two
     true dinucleotide (MNP) RYR1 sites in the panel. This avoids GATK entangling unrelated
     nearby SNPs elsewhere in the panel into false MNP calls while still correctly calling
     the two RYR1 sites that actually need it.
   - If `--ploidy-metrics` is given and reports the sample as XY, chrX (the G6PD gene,
     the only sex-chromosome gene in the panel) is forcecalled separately with
     `--sample-ploidy 1` and spliced in, so PharmCAT reports G6PD as hemizygous instead
     of a diploid homozygous call. Without this flag, or for any non-XY result, chrX stays
     diploid (default GATK behavior).
2. **Optional DRAGEN VCF-based refinement** of force-called genotypes
3. **VCF preprocessing** using the PharmCAT VCF preprocessor
4. **QUAL=inf fixing** for downstream compatibility
5. **PharmCAT report generation**
6. **mpileup-based annotation** for the predefined pharmacogenomic loci TSV bundled inside the image

## DRAGEN refinement logic

When a DRAGEN VCF is provided, force-called PharmCAT genotypes are refined using DRAGEN-supported calls.

The refinement is position-based and is intended to reduce unsupported non-reference force-calls. If a force-called non-reference genotype is not supported by a PASS non-reference DRAGEN call at the same genomic position, the genotype is converted to reference. Other FORMAT fields such as `AD`, `DP`, `GQ`, and `PL` are kept unchanged.

## Build

Run this in the repository directory:

```bash
docker build --no-cache -t pharmsuite:latest .
```

## Run

```bash
docker run --rm \
  -v /path/to/input_data:/data \
  -v /path/to/output_dir:/out \
  mertcdll/pharmsuite:0.1.2 \
    python3 /usr/local/bin/pharmcat_pipeline.py \
      --cram /data/sample.cram \
      --reference /data/reference.fa \
      --outdir /out \
      --sample-id SAMPLE_ID \
      --dragen-vcf /data/sample.hard-filtered.vcf.gz \
      --ploidy-metrics /data/sample.ploidy_estimation_metrics.csv
```

`--dragen-vcf` and `--ploidy-metrics` are both optional:

- `--dragen-vcf`: DRAGEN `hard-filtered.vcf.gz` for the same sample. See "DRAGEN
  refinement logic" above for the exact policy. Omitting it skips this filter.
- `--ploidy-metrics`: DRAGEN `ploidy_estimation_metrics.csv` for the same sample. Used
  only to correctly forcecall chrX (G6PD) as hemizygous for XY samples, see above.
  Omitting it, or any non-XY value, keeps chrX diploid.

## Run mpileup annotation only

Use this to check each predefined variant's depth and allele fraction.

```bash
docker run --rm \
  -v /path/to/input_data:/data \
  -v /path/to/output_dir:/out \
  mertcdll/pharmsuite:0.1.2 \
    python3 /usr/local/bin/annotate_mpileup.py \
      --cram /data/sample.cram \
      --reference /data/reference.fa \
      --outdir /out \
      --sample-id SAMPLE_ID
```
