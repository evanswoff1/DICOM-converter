from __future__ import annotations
import argparse
import json
import logging  
import os
import re
import shutil
import tempfile
from pathlib import Path 
from typing import Dict, Tuple, List, Optional  
from concurrent.futures import ThreadPoolExecutor, as_completed 
from collections import defaultdict
import datetime
import io
from contextlib import redirect_stdout, redirect_stderr

import dicom2nifti
import dicom2nifti.settings as dcm_settings 
import nibabel as nib 
import numpy as np 
import pydicom

import time

from dicom2nifti.exceptions import ConversionValidationError
from dataclasses import dataclass

__version__ = "1.0.1"

# ============================================
# Logging
# ============================================
LOGGER = logging.getLogger("dicom_to_bids")

def setup_logging(verbosity: int):
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    
    logging.basicConfig(level = level, format = "[%(levelname)s] %(message)s")
    
# ============================================
# dicom2nifti configuration
# ============================================
def configure_dicom2nifti(*, resample: bool = False):
    dcm_settings.disable_validate_slice_increment()
    dcm_settings.disable_validate_orientation()
    
    # Resampling can be a major runtime cost
    if resample:
        dcm_settings.enable_resampling()
        if hasattr(dcm_settings, "set_resample_padding"):
            dcm_settings.set_resample_padding(-1024)
        elif hasattr(dcm_settings, "set_resampling_padding"):
            dcm_settings.set_resampling_padding(-1024)
        else:
            LOGGER.warning("dicom2nifti has no resample padding setter; continuing without it.")
    else:
        # If available, explicitly disable resampling to ensure default behavior
        if hasattr(dcm_settings, "disable_resampling"):
            try:
                dcm_settings.disable_resampling()
            except Exception:
                pass

# ============================================
# Helpers
# ============================================
SKIP_MODALITIES = {"SEG", "SR", "PR", "KO", "OT", "SM", "ANN"}
ALLOW_MODALITIES = {"MR", "CT", "PT"}

def looks_like_localizer(ds) -> bool:
    sd = (getattr(ds, "SeriesDescription", "") or "").upper()
    it = getattr(ds, "ImageType", None)
    it_str = " ".join(it).upper() if isinstance(it, (list, tuple)) else str(it or "").upper()
    return ("LOCALIZER" in sd) or ("LOCALIZER" in it_str)

def subject_from_path(path: Path) -> str:
    for part in path.parts:
        if part.startswith('sub-'):
            return part
    return 'sub-unk'

def ensure_bids_subfolders(subject_root: Path) -> None:
    for name in ("anat","dwi","func","ct","pet"):
        (subject_root / name).mkdir(parents = True, exist_ok = True)

def sidecar_path_for_nifti(nifti_path: Path) -> Path:
    if nifti_path.name.endswith(".nii.gz"):
        return nifti_path.with_name(nifti_path.name[:-7] + ".json")
    else:
        return nifti_path.with_suffix(".json")

def _sec(val):
    if val is None:
        return None
    try:
        return float(val) / 1000.0
    except Exception:
        return None

def robust_rmtree(path: Path, retries: int = 12, base_sleep: float = 0.2):
    for i in range(retries):
        try:
            shutil.rmtree(path)
            return
        except PermissionError:
            time.sleep(base_sleep * (1.6 ** i))
        except FileNotFoundError:
            return
    
    LOGGER.warning(f"Could not delete temp dir (file lock persists): {path}")

def make_work_dir(output_root: Path) -> Path:
    tmp_root = output_root / ".tmp"
    tmp_root.mkdir(parents = True, exist_ok = True)
    return Path(tempfile.mkdtemp(dir = str(tmp_root)))

# ============================================
# Modality detection 
# ============================================
def _series_text(ds) -> str:
    parts = []
    for attr in ("SeriesDescription","ProtocolName","SequenceName","ScanningSequence","SequenceVariant"):
        try:
            v = getattr(ds, attr, None)
        except Exception:
            v = None
        if v:
            parts.append(str(v))
    try:
        it = getattr(ds, "ImageType", None)
        if it:
            if isinstance(it, (list, tuple)):
                parts.extend([str(x) for x in it if x])
            else:
                parts.append(str(it))
    except Exception:
        pass
    return " ".join(parts).upper()

def _get_first_number(x):
    try:
        if x is None:
            return None
        # pydicom MultiValue -> list-like
        if isinstance(x, (list, tuple)) and x:
            x = x[0]
        return float(x)
    except Exception:
        return None

def _looks_like_epi(text: str) -> bool:
    # EPI / BOLD / fMRI cues (vendors vary a lot)
    epi_keys = (" EPI","EPI","BOLD","FMRI","REST","TASK","FUNC","MB","MULTIBAND")
    return any(k in text for k in epi_keys)

def _looks_like_diffusion(ds, text: str) -> bool:
    # prefer explicit tags if present; otherwise fall back to text cues
    b = None
    for attr in ("DiffusionBValue","DiffusionBvalue","BValue"):
        b = _get_first_number(getattr(ds, attr, None))
        if b is not None:
            break
    if b is not None and b > 0:
        return True
    # some datasets encode diffusion in the description/protocol
    return any(k in text for k in ("DIFF","DWI","DTI","TRACE","ADC"))

def _norm_text(*parts: object) -> str:
    s = " ".join([str(p) for p in parts if p is not None])
    return " ".join(s.lower().replace("_"," ").replace("-", " ").split())

def infer_task_from_series(ds) -> str:
    desc = _norm_text(getattr(ds, "SeriesDescription", ""), getattr(ds, "ProtocolName", ""))
    if any(k in desc for k in ["rest", "rsfmri","rs fmri","resting"]):
        return "rest"
    
    m = re.search(r"\btask\s*([a-z0-9]+)\b", desc)
    if m:
        return m.group(1)
    if any(k in desc for k in ["bold", "fmri", "func"]):
        return "rest"
    return ""

def detect_modality(ds) -> Tuple[str, str]:
    mod = str(getattr(ds, "Modality", "") or "").upper().strip()
    
    # -- CT / PET --
    if mod == "CT":
        return ("ct", "ct")
    if mod in {"PT", "PET"}:
        return ("pet", "pet")
    
    # only MR is routed to func/dwi/anat with MR-specific heuristics
    if mod != "MR":
        return ("anat", "unknown")
    
    # MR heuristics
    # use multiple fields because vendors vary
    desc = _norm_text(getattr(ds, "SeriesDescription", ""), getattr(ds, "ProtocolName", ""))
    imgtype = _norm_text(getattr(ds, "ImageType", ""))
    seq = _norm_text(getattr(ds, "SequenceName", ""), getattr(ds, "ScanningSequence", ""), getattr(ds, "SequenceVariant", ""))
    
    # Pull TR/TE/TI if present (usually in ms in DICOM)
    te = getattr(ds, "EchoTime", None)
    tr = getattr(ds, "RepetitionTime", None)
    try:
        tr = float(tr) if tr is not None else None
    except Exception:
        tr = None
    try:
        te = float(te) if te is not None else None
    except Exception:
        te = None
    
    tr_s = None
    if tr is not None:
        tr_s = tr / 1000.0 if tr > 50 else tr
    
    has_time = False
    for tag in ("NumberOfTemporalPositions", "NumberOfFrames"):
        v = getattr(ds, tag, None)
        try:
            if v is not None and int(v) > 1:
                has_time = True
                break
        except Exception:
            pass
    
    # --- DWI (diffusion) ---
    bval = getattr(ds, "DiffusionBValue", None)
    if bval is not None:
        try:
            if float(bval) > 0:
                return ("dwi", "dwi")
        except Exception:
            pass
    if any(k in desc for k in ["dwi","diff","diffusion","trace","adc","b0"]) or "diffusion" in imgtype:
        return "dwi", "dwi"
    
    # --- fMRI/BOLD ---
    bold_keywords = ["bold","fmri","rsfmri","rest","task","func"]
    epi_keywords = ["epi","ep2d","ge epi","se epi"]
    looks_bold = any(k in desc for k in bold_keywords) or any(k in imgtype for k in ["bold","fmri"])
    looks_epi = any(k in desc for k in epi_keywords) or any(k in seq for k in epi_keywords) or "epi" in imgtype
    
    if looks_bold:
        return ("func", "bold")
    if looks_epi and tr_s is not None and 0.3 <= tr_s <= 6.0 and (has_time or "mossaic" in imgtype or "fmri" in imgtype):
        return ("func", "bold")
    
    if "flair" in desc:
        return ("anat", "FLAIR")
    
    if te is not None:
        te_ms = te if te > 1 else te * 1000.0
        if te_ms < 20:
            return ("anat", "T1w")
        if te_ms >= 70:
            return ("anat", "T2w")
    if any(k in desc for k in ["t1","mprage","spgr","bravo","vibe","tf1"]):
        return ("anat", "T1w")
    if any(k in desc for k in ["t2","tse","fs","blade","space"]):
        return ("anat", "T2w")
    
    # Evidence is weak: keep unknown
    return "anat", "unknown"

# ============================================
# Echo handling
# ============================================
def extract_te(ds):
    try:
        return float(ds.EchoTime)
    except Exception:
        return None

def group_dicoms_by_te(dicom_files, tol = 1e-3):
    groups = {}
    
    for f in dicom_files:
        try:
            if hasattr(f, "is_file") and not f.is_file():
                continue
        except Exception:
            continue
        ds = safe_dcmread(f, stop_before_pixels = True, force = True, specific_tags = ["EchoeTime"])
        if ds is None:
            continue
        
        te = extract_te(ds)
        matched = False
        
        for k in groups:
            if te is not None and k is not None and abs(k - te) < tol:
                groups[k].append(f)
                matched = True
                break
    
        if not matched:
            groups[te] = [f]
    
    return groups

def validate_te_groups(groups, qc_mode = "warn"):
    if len(groups) <= 1:
        return
    
    counts = [len(v) for v in groups.values()]
    if len(set(counts)) != 1:
        msg = f"Inconsistent volume counts across echoes: {counts}"
        if qc_mode == "strict":
            raise RuntimeError(msg)
        LOGGER.warning(msg)

def is_true_multi_echo(dicom_files: List[Path], *, min_slices: int = 3, sample_limit: int = 300) -> bool:
    # Quick TE scan (header only)
    te_counts: Dict[float, int] = {}
    slice_to_tes: Dict[str, set] = {}
    scanned = 0
    
    for fp in dicom_files:
        if scanned >= sample_limit:
            break
        scanned += 1
        ds = safe_dcmread(fp, stop_before_pixels = True)
        if ds is None:
            continue
        
        te = extract_te(ds)
        if te is None:
            continue
        # Round to avoid tiny floating-point/vendor differences
        te_r = round(float(te), 2)
        te_counts[te_r] = te_counts.get(te_r, 0) + 1
        
        # slice key: prefer ImagePosition Patient z, else SliceLocation, else InstanceNumber
        z_key = None
        ipp = getattr(ds, "ImagePositionPatient", None)
        if ipp is not None and isinstance(ipp, (list, tuple)) and len(ipp) >= 3:
            try:
                z_key = float(ipp[2])
            except Exception:
                try:
                    z_key = float(ipp[2])
                except Exception:
                    z_key = None
        if z_key is None:
            z_key = getattr(ds, "SliceLocation", None)
        if z_key is None:
            z_key = getattr(ds, "InstanceNumber", None)
        
        sk = str(z_key)
        slice_to_tes.setdefault(sk, set()).add(te_r)
    
    if len(te_counts) < 2:
        return False
    
    # require enough data per TE group (using the sampled headers)
    if any(c < int(min_slices) for c in te_counts.values()):
        return False
    
    # require repeated slice positions with >1 TE 
    multi_slice = sum(1 for tes in slice_to_tes.values() if len(tes) >= 2)
    total_slices = len(slice_to_tes)
    # heurestic: at least 30% of slice positions show multiple echoes
    if total_slices == 0:
        return False
    if multi_slice < max(int(min_slices), int(0.3 * total_slices)):
        return False
    
    return True

def sorted_echo_groups(groups):
    return sorted(groups.items(), key = lambda x: (x[0] is None, x[0]))

# ============================================
# Slice timing and PE direction
# ============================================
def extract_slice_timing(ds):
    val = getattr(ds, "SliceTiming", None)
    if val is None:
        return None
    try:
        return [float(v) / 1000.0 for v in val]
    except Exception:
        return None

def extract_phase_encoding(ds):
    direction = getattr(ds, "InPlanePhaseEncodingDirection", None)
    readout = getattr(ds, "TotalReadoutTime", None)
    return direction, readout

# ============================================
# JSON metadata
# ============================================
def extract_json_metadata(ds) -> Dict:
    meta = {
        'Modality': getattr(ds, 'Modality', None),
        'Manufacturer': getattr(ds, 'Manufacturer', None),
        'MagneticFieldStrength': getattr(ds, 'MagneticFieldStrength', None),
        'RepetitionTime': _sec(getattr(ds, 'RepetitionTime', None)),
        'EchoTime': _sec(getattr(ds, 'EchoTime', None)),
        'InversionTime': _sec(getattr(ds, 'InversionTime', None)),
        'FlipAngle': getattr(ds, 'FlipAngle', None),
    }
    
    st = extract_slice_timing(ds)
    if st is not None:
        meta['SliceTiming'] = st
    
    pe, readout = extract_phase_encoding(ds)
    
    if pe:
        meta["PhaseEncodingDirection"] = pe 
    if readout:
        meta["TotalReadoutTime"] = readout
    
    return meta

# ============================================
# DWI validation
# ============================================
def validate_dwi_gradients(bvals, bvecs, qc_mode = 'warn'):
    if bvals is None or bvecs is None:
        msg = "Missing bvals or bvecs"
        if qc_mode == 'strict':
            raise RuntimeError(msg)
        LOGGER.warning(msg)
        return
    
    bvals = np.asarray(bvals).astype(float)
    bvecs = np.asarray(bvecs).astype(float)
    
    if bvecs.ndim != 2 or bvecs.shape[0] != 3:
        msg = f"Invalid bvec shape: {bvecs.shape}"
        if qc_mode =='strict':
            raise RuntimeError(msg)
        LOGGER.warning(msg)
    
    if bvecs.shape[1] != len(bvals):
        msg = "Mismatch between number of bvals and bvecs"
        if qc_mode == 'strict':
            raise RuntimeError(msg)
        LOGGER.warning(msg)
    
    norms = np.linalg.norm(bvecs, axis = 0)
    bad = (norms > 0) & ((norms < 0.5) | (norms > 1.5))
    
    if np.any(bad):
        msg = "Non-unit diffesion gradient vectors detected"
        if qc_mode == 'strict':
            raise RuntimeError(msg)
        LOGGER.warning(msg)
            
# ============================================
# Echo stacking
# ============================================
def stack_echoes(nifti_paths, output_path):
    imgs = [nib.load(str(p), mmap = False) for p in nifti_paths]
    data = [img.get_fdata(dtype = np.float32) for img in imgs]
    
    shapes = [d.shape[:3] for d in data]
    if len(set(shapes)) != 1:
        raise RuntimeError("Echo volumes have mismatched spatial shapes")
    
    stacked = np.stack(data, axis = -1)
    out = nib.Nifti1Image(stacked, imgs[0].affine, imgs[0].header)
    nib.save(out, str(output_path))

# ============================================
# QC
# ============================================
def run_qc(nifti_path: Path, qc_mode: str = "warn") -> nib.Nifti1Image:
    
    qc_mode = (qc_mode or "warn").lower().strip()
    img = nib.load(str(nifti_path))
    
    if qc_mode == "off":
        return img
    
    def _fail(msg: str):
        if qc_mode == "strict":
            raise RuntimeError(msg)
        LOGGER.warning(msg)
    
    # shape checks
    shape = tuple(int(x) for x in img.shape)
    if len(shape) < 3 or any(d < 1 for d in shape[:3]):
        _fail(f"QC: invalid NIfTI shape {shape} for {nifti_path}")
    
    zooms = img.header.get_zooms()
    if len(zooms) >= 3:
        if any((z is None) or (float(z) <= 0) for z in zooms[:3]):
            _fail(f"QC: non-positive zooms {zooms} for {nifti_path}")
    
    # affine checks
    try:
        aff = img.affine
        if aff.shape != (4, 4):
            _fail(f"QC: affine has shape {aff.shape} for {nifti_path}")
        else:
            import numpy as _np

            if not _np.isfinite(aff).all():
                _fail(f"QC: affine contains non-finite values for {nifti_path}")
    except Exception as e:
        _fail(f"QC: failed affine inspection for {nifti_path}: {e}")
    
    # data fitness check
    try:
        import numpy as _np
        
        dataobj = img.dataobj
        if len(shape) >= 3:
            xs = slice(0, min(shape[0], 4))
            ys = slice(0, min(shape[1], 4))
            zs = slice(0, min(shape[2], 4))
            sample = _np.asanyarray(dataobj[xs, ys, zs])
        else:
            sample = _np.asanyarray(dataobj)
        if not _np.isfinite(sample).all():
            _fail(f"QC: non-finite voxel values detected in {nifti_path}")
    except Exception as e:
        _fail(f"QC: failed data sampling for {nifti_path}: {e}")
    
    try:
        qcode = int(img.header.get_qform(coded = True)[1])
        scode = int(img.header.get_sform(coded = True)[1])
        if qcode == 0 and scode == 0:
            _fail(f"QC: both qform and sform codes are 0 for {nifti_path}")
    except Exception:
        pass
    
    return img

# ============================================
# DWI export
# ============================================
def nifti_stem_without_nii_gz(p: Path) -> Path:
    if p.name.endswith(".nii.gz"):
        return p.with_name(p.name[:-7])
    if p.suffix == ".nii":
        return p.with_suffix("")
    return p.with_suffix("")

def export_bvals_bvecs(nifti_path: Path, out_base: Path, qc_mode = 'warn'):
    base = nifti_stem_without_nii_gz(nifti_path)
    for ext in ['.bval', '.bvec']:
        src = base.with_suffix(ext)
        if src.exists():
            dst = out_base.with_suffix(ext)
            try:
                src.replace(dst)
            except PermissionError:
                shutil.copy2(src, dst)
                try:
                    src.unlink()
                except Exception:
                    pass
    
    bval = out_base.with_suffix('.bval')
    bvec = out_base.with_suffix('.bvec')
    if bval.exists() and bvec.exists():
        validate_dwi_gradients(np.loadtxt(bval), np.loadtxt(bvec), qc_mode)

# ============================================
# Conversion core
# ============================================
@dataclass
class OneSeriesSummary:
    dicom_dir: str
    converted: int
    skipped_uids: int
    written: List[str]
    errors: List[str]

def convert_one_series(dicom_dir: Path, output_root: Path, task: str, qc_mode: str, stack_echoes_flag: bool) -> OneSeriesSummary:
    # Collect DICOM files in this directory
    files = [p for p in dicom_dir.glob("**/*") if p.is_file()]
    if not files:
        LOGGER.warning(f"No files in {dicom_dir}")
        return OneSeriesSummary(str(dicom_dir), converted = 0, skipped_uids = 0, written = [], errors = [])
    
    # 1) group by SeriesInstanceUID (prevent mixed-series conversion failures)
    series_groups: Dict[str, List[Path]] = defaultdict(list)
    series_meta: Dict[str, pydicom.Dataset] = {} 
    errors: List[str] = []  
    first_unreadable_err: Optional[str] = None
    
    for f in files:
        ds = safe_dcmread(f, stop_before_pixels = True, force = True, specific_tags = ["Modality","SeriesInstanceUID","SeriesDescription","ImageType","EchoTime","RepetitionTime","FlipAngle"])
        if ds is None:
            continue
        
        modality = str(getattr(ds, "Modality", "")).upper()
        if modality in {"SEG","SR","ANN","SM","OT","PR","RTSTRUCT","RTDOSE","RTPLAN","REG"}:
            continue
        
        uid = getattr(ds, "SeriesInstanceUID", None)
        if not uid:
            continue
        
        uid = str(uid)
        series_groups[uid].append(f)
        
        # Keep best representative DICOM
        if uid not in series_meta:
            series_meta[uid] = ds
        else:
            score_new = sum(1 for tag in ("EchoTime", "RepetitionTime", "FlipAngle") if hasattr(ds, tag))
            score_old = sum(1 for tag in ("EchoTime", "RepetitionTime", "FlipAngle") if hasattr(series_meta[uid], tag))
            if score_new > score_old:
                series_meta[uid] = ds
    
    if not series_groups:
        LOGGER.warning(f"No readable/usable DICOM image slices in {dicom_dir} (often SEG/SR/ANN/localizers or nested files)")
        return OneSeriesSummary(str(dicom_dir), converted = 0, skipped_uids = 0, written = [], errors = [])
    
    subject = subject_from_path(dicom_dir)
    # always create the expected BIDS-ish subfolders under the subject root
    ensure_bids_subfolders(output_root / subject)
    
    converted_any = 0
    skipped_uids = 0
    written_all: List[str] = []
    
    # 2) Convert each true series separately
    for uid, files_all in series_groups.items():
        # skip localizers / single-slice / scout series (dicom2nifti requires >=3 slices)
        if len(files_all) < 3:
            LOGGER.info(f"Skip UID {uid}: too few slices (n={len(files_all)})")
            skipped_uids += 1
            continue
        
        ds_rep = series_meta[uid]
        try:
            bids_folder, suffix = detect_modality(ds_rep)
        except Exception as e:
            LOGGER.warning(f"Could not detect modality for UID {uid} in {dicom_dir}: {e}")
            skipped_uids += 1
            continue
        
        # Make run label stable per series
        series_number = getattr(ds_rep, "SeriesNumber", None)
        if series_number is not None:
            try:
                run_label = f"run-{int(series_number):02d}"
            except Exception:
                run_label = f"run-{uid[-6:]}"
        else:
            run_label = f"run-{uid[-6:]}"
        
        out_dir = output_root / subject / bids_folder
        out_dir.mkdir(parents = True, exist_ok = True)
        
        # 3) echo grouping
        sorted_groups: List[Tuple[Optional[float], List[Path]]]
        do_split_echoes = (bids_folder == "anat") and is_true_multi_echo(files_all, min_slices = 3)
        if do_split_echoes:
            try:
                te_groups = group_dicoms_by_te(files_all)
                validate_te_groups(te_groups, qc_mode)
                sorted_groups = sorted_echo_groups(te_groups)
            except Exception as e:
                LOGGER.info(f"Skip UID {uid}: TE grouping/validation failed ({e})")
                skipped_uids += 1
                continue
        else:
            sorted_groups = [(None, files_all)]
        
        written_this_uid: List[Path] = []

        for echo_idx, (_te, files_echo) in enumerate(sorted_groups, start = 1):
            work_dir = make_work_dir(output_root)
            try:
                for f in files_echo:
                    try:
                        # Hardlink when possible
                        # Falls back to copy if hardlinking isn't supported
                        link_or_copy(f, work_dir / f.name)
                    except Exception:
                        continue
            
                # Convert echo group
                try:
                    _buf_out = io.StringIO()
                    _buf_err = io.StringIO()
                    try:
                        with redirect_stdout(_buf_out), redirect_stderr(_buf_err):
                            dicom2nifti.convert_directory(str(work_dir), str(work_dir), compression = True)
                    except ConversionValidationError as e:
                        LOGGER.info(f"Skip UID {uid} echo {echo_idx}: dicom2nifti validation ({e})")
                        continue
                except Exception as e:
                    LOGGER.warning(f"Failed conversion UID {uid} echo {echo_idx} in {dicom_dir}: {e}")
                    errors.append(f"{dicom_dir} UID {uid} echo {echo_idx}: CONVERT {type(e).__name__}: {e}")
                    continue
            
                nifti_files = sorted(work_dir.glob("*.nii*"), key = lambda p: p.stat().st_mtime, reverse = True)
            
                if not nifti_files:
                    log_txt = (_buf_out.getvalue() + "\n" + _buf_err.getvalue()).strip()
                    if ("TOO_FEW_SLICES" in log_txt) or ("At least 3 slices" in log_txt) or ("LOCALIZER" in log_txt):
                        LOGGER.info(f"Skip UID {uid} echo {echo_idx}: too few distinct slices/localizer (dir={dicom_dir})")
                    else:
                        LOGGER.warning(f"No NIfTI produced for UID {uid} echo {echo_idx} in {dicom_dir}" + (f" | dicom2nifti says: {log_txt[-300:]}" if log_txt else ""))
                    continue
            
                nifti = nifti_files[0]
                
                try:
                    img = run_qc(nifti, qc_mode)
                except Exception as e:
                    LOGGER.warning(f"QC failed for UID {uid} echo {echo_idx} in {dicom_dir}: {e}")
                    errors.append(f"{dicom_dir} UID {uid} echo {echo_idx}: QC {type(e).__name__}: {e}")
                    continue
            
                # 4) BIDS-ish naming
                base = f"{subject}_{run_label}"
                
                # ensure func series always get a reasonable task label
                task_use = task
                if bids_folder == "func":
                    task_use = (task or infer_task_from_series(ds_rep) or "rest")
                    base += f"_task-{task_use}"
            
                if len(sorted_groups) > 1:
                    base += f"_echo-{echo_idx}"
                base += f"_{suffix}"
            
                out_path = out_dir / f"{base}.nii.gz"
                nib.save(img, str(out_path))
            
                # sidecar JSON from a DICOM that belongs to echo group
                json_path = sidecar_path_for_nifti(out_path)
                ds_echo = safe_dcmread(files_echo[0], stop_before_pixels = True, force = True)
                if ds_echo is None:
                    ds_echo - ds_rep
                
                try:
                    with open(json_path, "w", encoding = "utf-8") as jf:
                        json.dump(extract_json_metadata(ds_echo), jf, indent = 2)
                except Exception as e:
                    LOGGER.warning(f"Failed writing JSON for {out_path}: {e}")
                    errors.append(f"{out_path}: JSON {type(e).__name__}: {e}")
            
                written_this_uid.append(out_path)
                written_all.append(str(out_path))
                converted_any += 1
                LOGGER.info(f"Wrote: {out_path}")
        
            finally:
                robust_rmtree(work_dir)
    
        if stack_echoes_flag and len(written_this_uid) > 1: 
            stack_path = out_dir / f"{subject}_{run_label}_{suffix}_echostack.nii.gz"
            try:
                stack_echoes(written_this_uid, stack_path)
                written_all.append(str(stack_path))
                LOGGER.info(f"Wrote: {stack_path}")
            except Exception as e:
                LOGGER.warning(f"Failed stacking echoes for UID {uid} in {dicom_dir}: {e}")
                errors.append(f"{dicom_dir} UID {uid}: stack {type(e).__name__}: {e}")
    
    if not converted_any:
        LOGGER.warning(f"No convertible series found in {dicom_dir} (likely localizers/too few slices).")
        
    return OneSeriesSummary(dicom_dir = str(dicom_dir), converted = converted_any, skipped_uids = skipped_uids, written = written_all, errors = errors)

# ============================================
# Dataset metadata
# ============================================
def write_dataset_description(out_root: Path):
    desc = { 
            "Name": "DICOM to BIDS Converted Dataset",
            "BIDSVersion": "1.9.0",
            "GeneratedBy": [{
                "Name": "dicom_to_nifti",
                "version": "1.0.1",
                "Description": "Automated DICOM -> BIDS converter with QC"
            }],
            "DatasetType": "raw",
            "HowToAcknowledge": "Please cite the software repository.",
            "generatedOn": datetime.datetime.now(datetime.UTC).isoformat()
            }
    
    path = out_root / "dataset_description.json"
    if not path.exists():
        with open(path, "w", encoding = "utf-8") as f:
            json.dump(desc, f, indent = 2)

# ============================================
# Batch execution
# ============================================
def safe_dcmread(fp: Path, *, stop_before_pixels: bool = True, force: bool = True, specific_tags = None, **kwargs):
    try:
        if specific_tags is not None:
            try:
                return pydicom.dcmread(str(fp), stop_before_pixels = stop_before_pixels, force = force, specific_tags = specific_tags, **kwargs)
            except TypeError:
                return pydicom.dcmread(str(fp), stop_before_pixels = stop_before_pixels, force = force, **kwargs)
        return pydicom.dcmread(str(fp), stop_before_pixels = stop_before_pixels, force = force, **kwargs)
    except Exception:
        return None

def link_or_copy(src: Path, dst: Path):
    try:
        os.link(src, dst)
    except Exception:
        shutil.copy2(src, dst)

def _slice_distinct_key(ds):
    ipp = getattr(ds, "ImagePositionPatient", None)
    if ipp is not None and isinstance(ipp, (list, tuple)) and len(ipp) >= 3:
        try:
            return float(ipp[2])
        except Exception:
            try:
                return str(ipp[2])
            except Exception:
                pass
    sl = getattr(ds, "SliceLocation", None)
    if sl is not None:
        return sl
    inst = getattr(ds, "InstanceNumber", None)
    if inst is not None:
        return inst
    return getattr(ds, "SOPInstanceUID", None)

def find_series_dirs(dicom_root: Path, *, modalities: tuple[str, ...] = ("MR","CT","pt"), min_slices: int = 3, max_files_to_scan: int = 200) -> list[Path]:
    dicom_root = Path(dicom_root)
    if not dicom_root.exists():
        return []
    
    # normalize modalities to a set
    want = {str(m).upper().strip() for m in modalities if str(m).strip()}
    if not want:
        want = {"MR", "CT", "PT"}
    
    out: list[Path] = []
    
    for root, subdirs, files in os.walk(dicom_root):
        subdirs.sort()
        files.sort()
        
        if not files:
            continue
        
        root_p = Path(root)
        
        distinct: set[object] = set()
        series_uid: Optional[str] = None
        convertible_seen = 0
        mixed = False
        
        scanned = 0
        for fn in files:
            if scanned >= int(max_files_to_scan):
                break
            fp = root_p / fn
            if not fp.is_file() or not is_dicom_file(fp):
                continue
            
            scanned += 1
            
            # IMPORTANT: call safe_dcmread without keyword to avoid "unexpected keyword argument 'stop_before_pixels'"
            ds = safe_dcmread(fp, stop_before_pixels = True)
            if ds is None:
                continue
            
            
            modality = str(getattr(ds, "Modality", "") or "").upper().strip()
            if modality not in want:
                continue
            
            if not is_convertible_image_slice(ds, modalities = want):
                continue
            
            convertible_seen += 1
            
            suid = str(getattr(ds, "SeriesInstanceUID", "") or "")
            if not suid:
                continue
            
            if series_uid is None:
                series_uid = suid
            elif suid != series_uid:
                mixed = True
                break
            
            # multi-frame DICOM can represent a whole volume in one file
            try:
                nframes = int(getattr(ds, "NumberOfFrames", 0) or 0)
            except Exception:
                nframes = 0
            if nframes >= int(min_slices):
                distinct = {0, 1, 2}
                break
            
            distinct.add(_slice_distinct_key(ds))
            
            if len(distinct) >= int(min_slices):
                break
        
        if (not mixed) and series_uid and convertible_seen >= int(min_slices) and len(distinct) >= int(min_slices):
            out.append(root_p)
    
    # de-dupe while preserving order, then sort for stable output
    seen: set[str] = set()
    unique: list[Path] = []
    for p in out:
        rp = str(p.resolve())
        if rp not in seen:
            seen.add(rp)
            unique.append(p)
    return sorted(unique)

def is_dicom_file(fp: Path) -> bool:
    if not fp.is_file():
        return False
    try:
        from pydicom.misc import is_dicom
        return bool(is_dicom(str(fp)))
    except Exception:
        try:
            with open(fp, "rb") as f:
                f.seek(128)
                return f.read(4) == b"DICM"
        except Exception:
            return False

def is_convertible_image_slice(ds, *, modalities: Optional[set[str]] = None) -> bool:
    if ds is None:
        return False
    
    mod = str(getattr(ds, "Modality","") or "").upper().strip()
    if not mod:
        return False
    
    if mod in SKIP_MODALITIES:
        return False
    
    allowed = modalities if modalities is not None else ALLOW_MODALITIES
    if mod not in allowed:
        return False
    
    # exclude obvious localizers/scouts
    if looks_like_localizer(ds):
        return False
    
    # must look like an image (rows/cols) or be multiframe
    rows = getattr(ds, "Rows", None)
    cols = getattr(ds, "Columns", None)
    if rows is None or cols is None:
        # allow multiframe volumes
        nframes = getattr(ds, "NumberOfFrames", None)
        try:
            if nframes is not None and int(nframes) > 1:
                return True
        except Exception:
            pass
        return False
    
    return True

def batch_convert(input_root: Path, output_root: Path, task: str = "rest", workers: int = 0, qc_mode: str = "warn", stack_echoes_flag: bool = False, *, max_workers: Optional[int] = None, report_jsonl: Optional[Path] = None) -> List[OneSeriesSummary]:
    
    dicom_root = Path(input_root)
    output_root = Path(output_root)
    output_root.mkdir(parents = True, exist_ok = True)
    
    if max_workers is None:
        if max_workers is None or int(workers) <= 0:
            max_workers = max(1, (os.cpu_count() or 4))
        else:
            max_workers = max(1, int(workers))
    
    candidates = find_series_dirs(dicom_root, modalities = ("MR","CT","PT"), min_slices = 3)
    LOGGER.info(f"Found {len(candidates)} candidate series directories under {dicom_root}")
    
    summaries: List[OneSeriesSummary] = []
    
    jf = open(report_jsonl, "w", encoding = "utf-8") if report_jsonl else None
    try:
        if max_workers and max_workers > 1:
            with ThreadPoolExecutor(max_workers = max_workers) as ex:
                futs = {ex.submit(convert_one_series, d, output_root, task, qc_mode, stack_echoes_flag): d for d in candidates}
                for fut in as_completed(futs):
                    d = futs[fut]
                    try:
                        summ = fut.result()
                    except Exception as e:
                        LOGGER.exception(f"Unhandled error converting {d}: {e}")
                        summ = OneSeriesSummary(str(d), converted = 0, skipped_uids = 0, written = [], errors = [f"UNHANDLED: {e}"])
                        
                    summaries.append(summ)
                    if jf:
                        jf.write(json.dumps(summ.__dict__) + "\n")
                        jf.flush()
        else:
            for d in candidates:
                try:
                    summ = convert_one_series(d, output_root, task, qc_mode, stack_echoes_flag)
                except Exception as e:
                    LOGGER.exception(f"Unhandled error converting {d}: {e}")
                    summ = OneSeriesSummary(str(d), converted = 0, skipped_uids = 0, written = [], errors = [f"UNHANDLED: {e}"])
                summaries.append(summ)
                if jf:
                    jf.write(json.dumps(summ.__dict__) + "\n")
    
    finally:
        if jf:
            jf.close()
    
    n_conv = sum(s.converted > 0 for s in summaries)
    LOGGER.info(f"Batch done. Converted dirs: {n_conv}/{len(summaries)}")
    
    return summaries

# ============================================
# Entry point
# ============================================
def main():
    parser = argparse.ArgumentParser(description="Parallel DICOM -> BIDS NIfTI with QC + DWI support")
    parser.add_argument("--input_root", required = True, type = Path)
    parser.add_argument("--output_root", required = True, type = Path)
    parser.add_argument("--task", default = 'rest')
    parser.add_argument("--workers", default = 'auto')
    parser.add_argument("--qc", default = 'warn', choices = ['off', 'warn', 'strict'])
    parser.add_argument("--resample", action = "store_true", help = "Enable dicom2nifti resampling (slower). Default: off")
    parser.add_argument("--stack-echoes", action = "store_true")
    parser.add_argument("-v", "--verbose", action = "count", default = 0, help = "Increase verbosity (-v, -vv)")
    parser.add_argument("--version", action = "version", version = f"%(prog)s {__version__}")
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    configure_dicom2nifti(resample = args.resample)
    workers = -1 if args.workers == 'auto' else int(args.workers)
    
    write_dataset_description(args.output_root)
    batch_convert(args.input_root, args.output_root, args.task, workers, args.qc, args.stack_echoes)
    
    LOGGER.info("Conversion complete")

if __name__ == "__main__":
    main()