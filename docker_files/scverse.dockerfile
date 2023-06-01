# Use Ubuntu 20.04 as base image
FROM python:3.8

# Non_interactive Install
ENV DEBIAN_FRONTEND=noninteractive 

# Metadata
LABEL base.image="ubuntu:20.04"
LABEL description="Scverse with Python 3.8"
LABEL maintainer="<spriyansh29@gmail.com>"

# Install Updates
RUN apt-get update && apt-get install -y --no-install-recommends

# Update pip
RUN python -m ensurepip --upgrade

# Set the working directory to root
WORKDIR /

# Create Required Folder
RUN mkdir -p input output script

# Install necessary packages
RUN pip install numpy pandas scipy scanpy anndata

# Set additional Parameters
#ENV BATCH_MEMORY=32000
#ENV BATCH_CPU=6

# Run when a container is launched
CMD ["/bin/bash"]
