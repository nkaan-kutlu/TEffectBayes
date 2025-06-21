# Use a base image with Python and R (Bioconductor)
FROM r-base:4.5.0

# Metadata Labels (must come after FROM)
LABEL maintainer="Necati Kaan Kutlu <nkaan.kutlu@gmail.com>"
LABEL description="Container with several standart R packages for bioinformatics workflows"
LABEL created="2025-04-22T12:00:00Z"
LABEL licenses="MIT"

# Silence debconf warnings
ARG DEBIAN_FRONTEND=noninteractive
ARG DEBCONF_NONINTERACTIVE_SEEN=true

# First fix the filesystem structure
RUN apt-get update && \
    apt-get install -y usrmerge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libpng-dev \
        libudunits2-dev \
        libgdal-dev && \
    rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN R -e "install.packages(c('optparse', 'stringr', 'tidyr', 'dplyr', 'purrr', 'readr', 'data.table'), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
RUN R -e "\
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager'); \
    BiocManager::install(); \
    BiocManager::install(c('rtracklayer', 'Rsubread', 'edgeR'))"

# Verify installations
RUN R -e "\
    packages <- c('optparse', 'stringr', 'tidyr', 'dplyr', 'purrr', 'readr', 'data.table', 'rtracklayer', 'Rsubread', 'edgeR'); \
    versions <- sapply(packages, function(pkg) packageVersion(pkg)); \
    print(versions)"

# Set working directory
WORKDIR /data

# Default command
CMD ["R"]
