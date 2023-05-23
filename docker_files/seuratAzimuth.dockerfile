# Use Ubuntu 20.04 as base image
FROM ubuntu:20.04

# Non_interactive Install
ENV DEBIAN_FRONTEND=noninteractive 

# Metadata
LABEL base.image="ubuntu:20.04"
LABEL description="Seurat in R Docker"
LABEL maintainer="<spriyansh29@gmail.com>"

# System dependencies for Seurat
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libmysqlclient-dev \
    libpng-dev \
    libmysqlclient-dev \
    libncurses-dev \
    libssl-dev \
    python3 \
    python3-pip \
    python3-dev \
    gfortran \
    perl \
    liblzma-dev \
    libpcre2-dev \
    libcurl4-openssl-dev \
    libreadline-gplv2-dev \
    libncursesw5-dev \
    libncurses5-dev \
    libssl-dev \
    libsqlite3-dev \
    tk-dev \
    libgdbm-dev \
    libc6-dev \
    libbz2-dev \
    libffi-dev \
    zlib1g-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    wget \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
    libcairo2-dev \
    libxt-dev

# Clean up to reduce Docker image size
RUN apt-get clean && rm -rf /var/lib/apt/lists/*
    
# Set environment variables for locale and PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
    
# Installing R
ENV R_VERSION=R-4.2.1
WORKDIR /usr/local
RUN wget https://cran.r-project.org/src/base/R-4/$R_VERSION.tar.gz && \
        tar xvf $R_VERSION.tar.gz && \
        rm $R_VERSION.tar.gz && \
        cd $R_VERSION && \
./configure --with-readline=no --with-x=no && make -j 4 && make install
ENV PATH=$PATH:~/.local/bin/ 
RUN R --version

# Installing dependencies
RUN Rscript -e 'install.packages("BiocManager", repos = "https://cran.us.r-project.org")'
RUN R -e "install.packages('Seurat', repos = 'http://cran.us.r-project.org')"

# Installing Additional Packages
RUN Rscript -e 'install.packages("tidyverse", repos = "https://cran.us.r-project.org", Ncpus = 4)'
RUN Rscript -e 'BiocManager::install("rhdf5", force = TRUE)'
RUN R -e "install.packages('remotes', repos = 'http://cran.us.r-project.org')"
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"
RUN R -e "install.packages('coop', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('Matrix', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e 'BiocManager::install("SingleR", force = TRUE)'
RUN Rscript -e 'BiocManager::install("celldex", force = TRUE)'

RUN R -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e 'devtools::install_github("satijalab/seurat-data", "seurat5")'
RUN Rscript -e 'devtools::install_github("satijalab/azimuth", "seurat5")'
RUN Rscript -e "remotes::install_github('satijalab/azimuth', ref = 'master')"


# Set the working directory to root
WORKDIR /

# Create Required Folder
RUN mkdir -p input output script

# Set additional Parameters
ENV BATCH_MEMORY=32000
ENV BATCH_CPU=6

# Run when a container is launched
CMD ["/bin/bash"]
