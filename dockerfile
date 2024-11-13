# Use an official Python runtime as a parent image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt update && \
    apt install -y ncbi-blast+ bowtie2 && \
    rm -rf /var/lib/apt/lists/*

# Copy the current directory contents into the container at /app
COPY . /app

# Install Python dependencies
RUN pip install tnseeker --no-cache-dir

# wait for cmd input
CMD ["bash"]
