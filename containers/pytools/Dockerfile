# Use a base image with Python and R (Bioconductor)
FROM python:3.10-slim

# Metadata Labels (must come after FROM)
LABEL maintainer="Necati Kaan Kutlu <nkaan.kutlu@gmail.com>"
LABEL description="pytools container for bioinformatics workflows"
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
    # Install Python packages with build dependencies
    pip install --no-cache-dir \
        pandas \
        numpy \
        bnlearn \
        scikit-learn\
        matplotlib \
        networkx \
        tqdm && \
        pip install --no-cache-dir \
        pyranges && \
    # Clean up build tools to reduce image size
    apt-get remove -y \
        build-essential \
        wget && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Verify installations
RUN python -c "import pandas; print(f'pandas version: {pandas.__version__}')" && \
    python -c "import numpy; print(f'numpy version: {numpy.__version__}')" && \
    python -c "import bnlearn; print(f'bnlearn version: {bnlearn.__version__}')" && \
    python -c "import sklearn; print(f'scikit-learn version: {sklearn.__version__}')" && \
    python -c "import matplotlib; print(f'matplotlib version: {matplotlib.__version__}')" && \
    python -c "import networkx; print(f'networkx version: {networkx.__version__}')" &&\
    python -c "import pyranges; print(f'pyranges version: {pyranges.__version__}')" && \
    python -c "import tqdm; print(f'tqdm version: {tqdm.__version__}')"


# Set working directory
WORKDIR /data

# Default command
CMD ["python"]
