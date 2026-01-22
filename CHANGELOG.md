# Changelog

## v1.0.1 - Initial Release

### Added
- Parallel DICOM -> BIDS-style NIfTI conversion pipeline
- Command-line interface (`dicom-to-nifti`)
- Automatic series discovery and grouping
- BIDS-style directory structure (`anat`, `func`, `dwi`, `ct`, `pet`)
- JSON sidecar generation
- Quality control modes: `off`, `warn`, `strict`
- Optional dicom2nifti resampling
- Optimized local-disk staging (hardlink-first with safe fallback)
- Support for MRI, CT, PET (including enhanced/multiframe when possible)

### Performance
- Tuned worker scheduling for SSD/NVMe storage
- Reduced redundant DICOM header reads
- Safe parallel execution across series

### Packaging
- MIT License
- `pyproject.toml`-based build system
- Wheel (`.whl`) and source (`.tar.gz`) distributions
- Zenodo DOI integration
- `CITATION.cff` for academic citation

---

Future changes will be documented here following semantic versioning.