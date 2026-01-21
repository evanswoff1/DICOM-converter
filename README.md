# dicom-to-nifti (BIDS-oriented)

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

## License

MIT License. See `LICENSE`.

---

## Citation

If you use this software in academic work, please tag a release and cite the repository version.