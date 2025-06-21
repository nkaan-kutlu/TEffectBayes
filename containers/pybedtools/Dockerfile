# Use a base image with Python and R (Bioconductor)
FROM python:3.10-slim

# Metadata Labels (must come after FROM)
LABEL maintainer="Necati Kaan Kutlu <nkaan.kutlu@gmail.com>"
LABEL description="pybedtools container for bioinformatics workflows"
LABEL created="2025-04-22T12:00:00Z"
LABEL licenses="MIT"

# Silence debconf warnings
ARG DEBIAN_FRONTEND=noninteractive
ARG DEBCONF_NONINTERACTIVE_SEEN=true

# Install system dependencies and bedtools
RUN apt-get update && \
    apt-get install -y apt-utils && \	
    apt-get install -y --no-install-recommends \
        wget \
        g++ \
        make \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev && \
    # Install bedtools 2.30.0
    wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
    tar -zxvf bedtools-2.30.0.tar.gz && \
    cd bedtools2 && \
    make -j12 && \
    mv bin/* /usr/local/bin/ && \
    cd .. && \
    rm -rf bedtools2 bedtools-2.30.0.tar.gz && \
    # Install Python packages with build dependencies
    pip install --no-cache-dir \
        cython \
        numpy && \
    pip install --no-cache-dir \
        pybedtools==0.11.0 \
        pandas && \
    # Clean up build tools to reduce image size
    apt-get remove -y \
        build-essential \
        wget && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*


# Install Python packages
RUN pip install --no-cache-dir \
    pybedtools==0.11.0 \
    pandas \
    numpy

# Verify installations
RUN bedtools --version && \
    python -c "import pybedtools; print(f'pybedtools version: {pybedtools.__version__}')" && \
    python -c "import pandas; print(f'pandas version: {pandas.__version__}')" && \
    python -c "import numpy; print(f'numpy version: {numpy.__version__}')"



# Set working directory
WORKDIR /data

# Default command
CMD ["python"]
