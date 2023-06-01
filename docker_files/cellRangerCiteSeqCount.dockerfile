# Use Ubuntu 20.04 as base image and python 3.8
FROM python:3.8-slim-buster

# Metadata
LABEL base.image="python3.8_ubuntu:20.04"
LABEL description="CITE-Seq"
LABEL maintainer="<spriyansh29@gmail.com>"

# Set a working directory.
WORKDIR /app

# Install dependencies.
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    git

# Upgrade pip.
RUN pip install --upgrade pip

# Create Required Folder
RUN mkdir -p input output script input2

# Install CITE-seq-Count.
RUN pip install CITE-seq-Count

# Set additional Parameters
#ENV BATCH_MEMORY=48000
#ENV BATCH_CPU=12

# Run when a container is launched
CMD ["/bin/bash"]