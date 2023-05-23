# Use Ubuntu 20.04 as base image
FROM ubuntu:20.04

# Non_interactive Install
ENV DEBIAN_FRONTEND=noninteractive 

# Metadata
LABEL base.image="ubuntu:20.04"
LABEL description="Cell Ranger Version 7.1.0 Docker"
LABEL maintainer="<spriyansh29@gmail.com>"

# Set Version
ENV CELLRANGER_VERSION=7.1.0

# Install necessary packages
RUN apt-get update

# Clean up to reduce Docker image size
RUN apt-get clean && rm -rf /var/lib/apt/lists/*
    
# Set environment variables for locale and PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Download the cellranger-7.1.0 and keep it in the tools directory
COPY tools/cellranger-$CELLRANGER_VERSION.tar.gz /usr/local/

# Set Working Directory to Local
WORKDIR /usr/local/

# Extract Cell Ranger tar file
RUN tar -xzvf cellranger-$CELLRANGER_VERSION.tar.gz

# Remove the file
RUN rm cellranger-$CELLRANGER_VERSION.tar.gz

# Export
ENV PATH /usr/local/cellranger-$CELLRANGER_VERSION/:$PATH

# Create Required Folder
RUN mkdir -p input output script

# Set the working directory to root
WORKDIR /

# Set additional Parameters
ENV BATCH_MEMORY=48000
ENV BATCH_CPU=12

# Run when a container is launched
CMD ["/bin/bash"]
