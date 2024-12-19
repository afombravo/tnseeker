# Use an official Python runtime as a parent image
FROM python:3.10

# Set the working directory in the container
WORKDIR /data

# Install system dependencies
RUN apt-get update && \
    apt-get install -y ncbi-blast+ && \
    apt-get install -y bowtie2 && \
    rm -rf /var/lib/apt/lists/*

# Copy the current directory contents into the container at /app
COPY . /data

# Install tnseeker from PyPI
RUN pip install tnseeker

CMD ["bash"]
