Tabii, hafif güncellenmiş hali şöyle:

````markdown
# PharmSuite Docker Pipeline

A lightweight Docker-based pipeline for single-sample PharmCAT preprocessing, PharmCAT reporting, DRAGEN-assisted genotype refinement, and targeted mpileup annotation.

## What it does

This container runs the pipeline for **one sample at a time** using:

- a single `CRAM` file
- a reference genome FASTA
- a user-defined `sample_id`
- a single output directory
- optionally, a DRAGEN VCF for further genotype refinement

The image contains the required PharmCAT resources internally, so they do not need to be passed at runtime.

## Included steps

The pipeline performs:

1. **Force-calling PharmCAT positions** from the input CRAM using GATK HaplotypeCaller
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
````

## Run

```bash
docker run --rm \
  -v /path/to/input_data:/data \
  -v /path/to/output_dir:/out \
  pharmsuite:latest \
    python3 /usr/local/bin/pharmcat_pipeline.py \
      --cram /data/sample.cram \
      --reference /data/reference.fa \
      --outdir /out \
      --sample-id SAMPLE_ID
```

## Run with DRAGEN VCF refinement

```bash
docker run --rm \
  -v /path/to/input_data:/data \
  -v /path/to/output_dir:/out \
  pharmsuite:latest \
    python3 /usr/local/bin/pharmcat_pipeline.py \
      --cram /data/sample.cram \
      --reference /data/reference.fa \
      --dragen-vcf /data/sample.dragen.vcf.gz \
      --outdir /out \
      --sample-id SAMPLE_ID
```

## Run mpileup annotation only

Use this to check each predefined variant's depth and allele fraction.

```bash
docker run --rm \
  -v /path/to/input_data:/data \
  -v /path/to/output_dir:/out \
  pharmsuite:latest \
    python3 /usr/local/bin/annotate_mpileup.py \
      --cram /data/sample.cram \
      --reference /data/reference.fa \
      --outdir /out \
      --sample-id SAMPLE_ID
```

````
