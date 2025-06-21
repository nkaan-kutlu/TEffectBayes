# Use a base image with Python and R (Bioconductor)
FROM ghcr.io/bioconductor/bioconductor_docker:RELEASE_3_18

# Metadata Labels (must come after FROM)
LABEL maintainer="Necati Kaan Kutlu <nkaan.kutlu@gmail.com>"
LABEL description="SalmonTE container for bioinformatics workflows"
LABEL created="2025-04-21T12:00:00Z"
LABEL licenses="MIT"

# Silence debconf warnings
ARG DEBIAN_FRONTEND=noninteractive
ARG DEBCONF_NONINTERACTIVE_SEEN=true

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        unzip \
        bowtie2 \
        samtools \
        salmon \
        git \
        python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Python packages (Snakemake and dependencies)
RUN pip3 install --no-cache-dir snakemake docopt pandas

# Install R packages (tidyverse, DESeq2, etc.)
RUN Rscript -e \
    "install.packages(c('tidyverse', 'scales', 'WriteXLS', 'BiocManager'), repos='https://cloud.r-project.org/')" && \
    Rscript -e \
    "BiocManager::install('DESeq2', ask = FALSE)"

RUN echo 'Successfully built the container.'

# Install SalmonTE
WORKDIR /opt
RUN git clone https://github.com/hyunhwan-jeong/SalmonTE.git
WORKDIR /opt/SalmonTE
RUN chmod +x *.py *.R

# Set default working directory for user data
WORKDIR /data

# Set PATH (if needed)
ENV PATH="/opt/SalmonTE:${PATH}"

