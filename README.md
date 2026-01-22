# dicom-to-nifti (BIDS-oriented)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18332432.svg)](https://doi.org/10.5281/zenodo.18332432)

a parallel DICOM -> NIfTI converter designed for ML / research pipelines.
It scans a DICOM directory tree, groups series, applies basic QC, and writes BIDS-style outputs (`.nii.gz` + `.json`) into modality folders (e.g., `anat/`, `func/`, `dwi/`, `ct/`, `pet`).

This project is intended to be **reproducible** (versioned, installable, containerizable) and **practical** on real-world mixed datasets (MR/CT/PET + localizers + SR/SEG noise).

---

## Features

- Parallel conversion with configurable worker count
- BIDS-ish naming and folder routing (`anat`, `func`, `dwi`, `ct`, `pet`)
- Always writes **compressed** NIfTI (`.nii.gz`) when enabled in the program
- Optional resampling via dicom2nifti
- QC modes: `off`, `warn`, `strict`
- Safe DICOM header reads (uses specific tags where possible)
- Staging optimized for local disks (hardlink-first with fallback when supported)

---

## Installation

### Virtual environment (recommended)

From the directory containing `pyproject.toml` and `dicom_to_nifti.py`:

#### Windows (PowerShell)
```powershell
python -m venv .venv
.\.venv\Scripts\activate
python -m pip install -U pip
python -m pip install -e .
```

#### macOS / Linux
```bash
python -m venv .venv
source .venv/bin/active
python -m pip install -U pip
python -m pip install -e .
```

### Verify installation
```bash
dicom-to-nifti --help
dicom-to-nifti --version
```

If `dicom-to-nifti` is not found, ensure your venv is activated. On Windows you can also run the executable directly:
```powershell
.\.venv\Scripts\dicom-to-nifti.exe --help
```

## Containerized usage (Docker / Singularity)

Pre-built containers allow fully reporducible execution without managing local Python environments.

### Docker

Build the image from the repository root:
```bash
docker build -t dicom-to-nifti:1.0.1 .
```

Run the converter:
```bash
docker run --rm -it \
  -v /path/to/dicom:/data/in:ro \
  -v /path/to/output:/data/out |
  dicom-to-nifti:1.0.1 \
  --input_root /data/in --output_root /data/out --workers 6 --qc warn
```

### Singularity

Build the image:
```bash
singularity build dicom-to-nifti_1.0.0.sif Singularity.def
```

Run:
```bash
singularity run dicom-to-nifti_1.0.1.sif \
  --input_root /path/to/dicom --output_root /path/to/output
```

---

## Quick start

```bash
dicom-to-nifti \ 
  --input_root /path/to/dicom \
  --output_root /path/to/bids \
  --workers auto \
  --qc warn
```

Windows example:
```powershell
dicom-to-nifti --input_root .\dicom --output_root .\bids_out --workers 6 --qc warn
```

---

## Command-line options

Run:
```bash
dicom-to-nifti --help
```

Key options:
- `--input_root`: Root directory containing DICOMs
- `--output_root`: Output directory
- `--task`: Default task name for fMRI (optional)
- `--workers`: Integer or `auto`
- `--qc`: `off`, `warn`, or `strict`
- `--resample`: Enable dicom2nifti resampling (slower)
- `--stack-echoes`: Optional multi-echo handling
- `-v / -vv`: Verbose logging

---

## Output structure

```
output_root/
  sub-<id>/
    anat/
    func/
    dwi/
    ct/
    pet/
```

Series that are not usable image stacks (SEG, SR, localizers, etc.) are skipped.

---

## Performance notes

- NVMe SSD + modern CPU: start with 4-8 workers
- Antivirus or OneDrive syncing can affect runtime
- Resampling significantly increases runtime

---

## Citation

If you use this software in academic work, please tag a release and cite the repository version.

- DOI: https://doi.org/10.5281/zenodo.18332432 

a `CITATION.cff` file is included in this repository.

---

## License

MIT License. See `LICENSE`.