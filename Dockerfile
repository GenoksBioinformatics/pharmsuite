FROM python:3.10.20-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PYTHONPATH=/opt/pharmcat

ARG HTSLIB_VERSION=1.23.1
ARG SAMTOOLS_VERSION=1.23.1
ARG BCFTOOLS_VERSION=1.23.1
ARG GATK_VERSION=4.6.1.0

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    wget \
    unzip \
    build-essential \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir \
    colorama==0.4.6 \
    numpy==2.2.6 \
    packaging \
    pandas==2.3.3 \
    python-dateutil==2.9.0.post0 \
    pytz \
    six==1.17.0 \
    tzdata

WORKDIR /tmp/build

RUN wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure \
    && make -j"$(nproc)" \
    && make install

RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make -j"$(nproc)" \
    && make install

RUN wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && make -j"$(nproc)" \
    && make install

RUN wget -q -O /tmp/gatk.zip https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
    && unzip -q /tmp/gatk.zip -d /opt \
    && ln -s /opt/gatk-${GATK_VERSION}/gatk /usr/local/bin/gatk \
    && rm -f /tmp/gatk.zip

RUN mkdir -p /opt/mpileup /opt/pharmcat /work

COPY pharmvariants_wfeatures.tsv /opt/mpileup/input_loci.tsv

COPY annotate_mpileup.py /usr/local/bin/annotate_mpileup.py
COPY pharmcat_pipeline.py /usr/local/bin/pharmcat_pipeline.py

COPY pharmcat-3.2.0-all.jar /opt/pharmcat/pharmcat-3.2.0-all.jar
COPY pharmcat_vcf_preprocessor /opt/pharmcat/pharmcat_vcf_preprocessor
COPY pcat /opt/pharmcat/pcat
COPY pharmcat_positions.sites.vcf.gz /opt/pharmcat/pharmcat_positions.sites.vcf.gz
COPY pharmcat_positions.sites.vcf.gz.tbi /opt/pharmcat/pharmcat_positions.sites.vcf.gz.tbi
COPY pharmcat_positions.vcf.bgz /opt/pharmcat/pharmcat_positions.vcf.bgz
COPY pharmcat_positions.vcf.bgz.csi /opt/pharmcat/pharmcat_positions.vcf.bgz.csi

RUN chmod +x \
    /usr/local/bin/annotate_mpileup.py \
    /usr/local/bin/pharmcat_pipeline.py \
    /opt/pharmcat/pharmcat_vcf_preprocessor \
    && ldconfig \
    && rm -rf /tmp/build

WORKDIR /work

CMD ["bash"]