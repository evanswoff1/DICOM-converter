# syntax = docker/dockerfile:1
# Lightweight, reproducible container for the dicom-to-nifti CLI
#
# Build:
#   docker build -t dicom-to-nifti:1.0.1 .
#
# Run (example):
#   docker run --rm -it \
#     -v "/path/to/dicom:/data/in:ro" \
#     -v "/path/to/out:/data/out" \
#     dicom-to-nifti:1.0.1 \
#     --input_root /data/in --output_root /data/out --workers 6 --qc warn

FROM python:3.11-slim
# Avoid Python writing .pyc files and enable unbuffered logs
ENV PYTHONDONTWRITEBYTECODE = 1 \
    PYTHONUNBUFFERED = 1
WORKDIR /opt/app

# Install OS-level basics (kept minimal). Add more only if you discover runtime needs.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Copy project metadata first for better layer caching
COPY pyproject.toml README.md LICENSE ./
# Copy the source (single-file module or package dir)
# If you later refactor into a package directory, add it here too.
COPY dicom_to_nifti.py ./
# If your CLI entry module is named differently, also copy it:
# COPY dicom_to_nifti.py ./

# Install your project into the image
RUN python -m pip install --upgrade pip \
 && python -m pip install .

# Default entrypoint to the CLI
ENTRYPOINT ["dicom-to-nifti"]
