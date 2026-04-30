# Single image for beCasKAS pipeline. Built and published to ghcr.io on tagged release.

FROM mambaorg/micromamba:1.5.8 AS base

LABEL org.opencontainers.image.source="https://github.com/GreenleafLab/beCasKAS"
LABEL org.opencontainers.image.description="beCasKAS Nextflow pipeline tools"
LABEL org.opencontainers.image.licenses="MIT"

USER root

# System deps for KAS-Analyzer wrapper, picard, gatk, etc.
RUN apt-get update && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        curl \
        wget \
        unzip \
        git \
        openjdk-17-jre-headless \
        procps \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libxt-dev \
        libcairo2-dev \
    && rm -rf /var/lib/apt/lists/*

USER mambauser

# Conda env with all bioinformatics tools (R, BWA, samtools, MACS2, etc.)
COPY --chown=mambauser:mambauser containers/env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes

# Activate env on shell entry
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# KAS-Analyzer 1.0 (clone + add to PATH)
USER root
RUN git clone --depth 1 --branch v1.0 https://github.com/Ruitulyu/KAS-Analyzer.git /opt/KAS-Analyzer \
    && chmod +x /opt/KAS-Analyzer/KAS-Analyzer
ENV PATH="/opt/KAS-Analyzer:${PATH}"

# Install R/BioConductor packages (DiffBind + deps)
COPY --chown=mambauser:mambauser containers/install-r-packages.R /tmp/install-r-packages.R
USER mambauser
RUN Rscript /tmp/install-r-packages.R

# Sanity check: verify each tool is on PATH
RUN bash -c "trim_galore --version && bowtie2 --version | head -1 && KAS-Analyzer --help >/dev/null && samtools --version | head -1 && macs2 --version && picard MarkDuplicates --version 2>&1 | head -1 && gencore --version 2>&1 | head -1 && gatk --version | head -1 && bcftools --version | head -1 && bedtools --version && bedGraphToBigWig 2>&1 | head -1 && Rscript -e 'library(DiffBind); cat(as.character(packageVersion(\"DiffBind\")), \"\\n\")'"

WORKDIR /workspace
ENTRYPOINT []
CMD ["/bin/bash"]
