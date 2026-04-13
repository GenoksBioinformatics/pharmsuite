# PharmSuite Docker Pipeline

A lightweight Docker-based pipeline for single-sample PharmCAT preprocessing, PharmCAT reporting, and targeted mpileup annotation.

## What it does

This container runs the pipeline for **one sample at a time** using:

- a single `CRAM` file
- a reference genome FASTA
- a user-defined `sample_id`
- a single output directory

The image contains the required PharmCAT resources internally, so they do not need to be passed at runtime.

## Included steps

The pipeline performs:

1. **Force-calling PharmCAT positions** from the input CRAM using GATK HaplotypeCaller
2. **VCF preprocessing** using the PharmCAT VCF preprocessor
3. **QUAL=inf fixing** for downstream compatibility
4. **PharmCAT report generation**
5. **mpileup-based annotation** for the predefined pharmacogenomic loci TSV bundled inside the image

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
  pharmsuite:latest \
  bash -lc '
    python3 /usr/local/bin/pharmcat_pipeline.py \
      --cram /data/sample.cram \
      --reference /data/reference.fa \
      --outdir /out \
      --sample-id SAMPLE_ID
  '
```

## Run mpileup annotation only (Check each variant's depth & AF)

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