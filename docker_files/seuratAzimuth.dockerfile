# Use Ubuntu 20.04 as base image
FROM rocker/r-base:4.2.1

# Non_interactive Install
ENV DEBIAN_FRONTEND=noninteractive 

# Metadata
LABEL base.image="ubuntu:20.04"
LABEL description="Seurat-Azimuth in R Docker"
LABEL maintainer="<spriyansh29@gmail.com>"

RUN apt-get update && apt-get install -y --no-install-recommends \
        xvfb xorg \
        xauth \
        xfonts-base \
        libx11-dev \
        systemd \
        libssl-dev \
        libhdf5-dev \
        libfontconfig1-dev \
	libxml2-dev libharfbuzz-dev \
	libfribidi-dev \
	libfreetype6-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev
        
# Installing cran packages
RUN Rscript -e 'install.packages("tidyverse", dependencies = T, Ncpus = 16)'
RUN Rscript -e 'install.packages("BiocManager", dependencies = T, Ncpus = 16)'
RUN Rscript -e 'install.packages("remotes", dependencies = T, Ncpus = 16)'
RUN Rscript -e 'install.packages("devtools", dependencies = T, Ncpus = 16)'
RUN Rscript -e 'install.packages("Seurat", dependencies = T, Ncpus = 16)'

# Installing From Bioconductor Packages
RUN Rscript -e 'BiocManager::install("rhdf5", force = TRUE, ask = FALSE)'

# Install 
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"
RUN R -e "remotes::install_github('satijalab/azimuth')"

# Set the working directory to root
WORKDIR /

# Create Required Folder
RUN mkdir -p input output script

# Set additional Parameters
#ENV BATCH_MEMORY=32000
#ENV BATCH_CPU=6

# Run when a container is launched
CMD ["/bin/bash"]

