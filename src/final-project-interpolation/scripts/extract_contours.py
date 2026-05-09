#!/usr/bin/env python3
# =============================================================================
# Final project: Geometric Interpolation of Tumor Contours
# Author: Abel Albuez Sanchez
#
# Stage 1 — extract 2D ET contours from BraTS 2024 GLI.
#
#  1. Take the first 100 cases from the BraTS zip (or the already-extracted
#     copy next to it) and place each case's *-seg.nii.gz under data/raw/.
#  2. For every axial slice with the ET label (label == 3), trace the largest
#     outer contour and write it to data/contours/<case>/slice_XXXX.obj.
#  3. Write data/contours/index.csv listing every (case, slice) pair so the
#     C++ stage knows which consecutive pairs exist.
# =============================================================================
from __future__ import annotations

import csv
import shutil
import sys
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np
from skimage import measure
from tqdm import tqdm

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

ET_LABEL          = 3
MIN_PIXELS        = 10
N_CASES           = 100
CASE_PREFIX       = "BraTS-GLI-"
SEG_SUFFIX        = "-seg.nii.gz"

PROJECT_ROOT      = Path(__file__).resolve().parents[1]
RAW_DIR           = PROJECT_ROOT / "data" / "raw"
CONTOURS_DIR      = PROJECT_ROOT / "data" / "contours"
INDEX_CSV         = CONTOURS_DIR / "index.csv"

DATASET_DIR       = Path(
    "/Users/abelalbuez/Documents/Maestria/Tercer Semestre/Tesis datasets"
)
DATASET_ZIP       = DATASET_DIR / "BraTS2024-BraTS-GLI-TrainingData.zip"


# -----------------------------------------------------------------------------
# Stage 1.1 — make sure data/raw/ has the first N_CASES BraTS cases.
# -----------------------------------------------------------------------------

def list_cases_in_zip(zip_path: Path) -> list[str]:
    """Return the sorted, deduplicated list of case ids inside the zip."""
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()

    cases: set[str] = set()
    for name in names:
        # Entries look like "BraTS-GLI-00000-000/BraTS-GLI-00000-000-seg.nii.gz"
        parts = name.split("/")
        for part in parts:
            if part.startswith(CASE_PREFIX) and not part.endswith(".nii.gz"):
                cases.add(part)
                break
    return sorted(cases)


def list_cases_extracted(root: Path) -> list[str]:
    """Return cases already sitting next to the zip as plain folders."""
    if not root.exists():
        return []
    return sorted(
        p.name for p in root.iterdir()
        if p.is_dir() and p.name.startswith(CASE_PREFIX)
    )


def copy_case_from_disk(src_case_dir: Path, dst_case_dir: Path) -> bool:
    """Copy only the segmentation file from an already-extracted case."""
    seg = src_case_dir / f"{src_case_dir.name}{SEG_SUFFIX}"
    if not seg.exists():
        return False
    dst_case_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(seg, dst_case_dir / seg.name)
    return True


def extract_case_from_zip(zf: zipfile.ZipFile, case_id: str,
                          dst_case_dir: Path) -> bool:
    """Extract only the segmentation file of one case from the open zip."""
    target_name = f"{case_id}{SEG_SUFFIX}"
    member = None
    for name in zf.namelist():
        if name.endswith(target_name):
            member = name
            break
    if member is None:
        return False

    dst_case_dir.mkdir(parents=True, exist_ok=True)
    with zf.open(member) as src, open(dst_case_dir / target_name, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return True


def stage_raw_cases() -> list[Path]:
    """Populate data/raw/ with the first N_CASES segmentation files."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    extracted_root = DATASET_DIR / "BraTS2024-BraTS-GLI-TrainingData"
    on_disk = list_cases_extracted(extracted_root)

    if on_disk:
        cases = on_disk[:N_CASES]
        print(f"Found {len(on_disk)} extracted cases — staging first "
              f"{len(cases)} from {extracted_root}")
        for case_id in tqdm(cases, desc="copying"):
            dst = RAW_DIR / case_id
            seg_dst = dst / f"{case_id}{SEG_SUFFIX}"
            if seg_dst.exists():
                continue
            if not copy_case_from_disk(extracted_root / case_id, dst):
                tqdm.write(f"  skip (no seg): {case_id}")
        return [RAW_DIR / c for c in cases]

    if not DATASET_ZIP.exists():
        sys.exit(f"BraTS source not found. Looked for:\n"
                 f"  {extracted_root}\n  {DATASET_ZIP}")

    cases = list_cases_in_zip(DATASET_ZIP)[:N_CASES]
    print(f"Staging {len(cases)} cases from zip {DATASET_ZIP.name}")
    with zipfile.ZipFile(DATASET_ZIP) as zf:
        for case_id in tqdm(cases, desc="extracting"):
            dst = RAW_DIR / case_id
            seg_dst = dst / f"{case_id}{SEG_SUFFIX}"
            if seg_dst.exists():
                continue
            if not extract_case_from_zip(zf, case_id, dst):
                tqdm.write(f"  skip (no seg in zip): {case_id}")
    return [RAW_DIR / c for c in cases]


# -----------------------------------------------------------------------------
# Stage 1.2 — extract per-slice ET contours and write them as .obj files.
# -----------------------------------------------------------------------------

def write_contour_obj(path: Path, contour_xy: np.ndarray, case_id: str,
                      slice_z: int) -> None:
    """Write one closed 2D polyline as a Wavefront .obj file."""
    n = contour_xy.shape[0]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(f"# Contorno tumor ET - caso {case_id} - slice axial {slice_z}\n")
        f.write(f"# {n} vertices\n")
        for x, y in contour_xy:
            f.write(f"v {x:.6f} {y:.6f}\n")
        for i in range(1, n):
            f.write(f"l {i} {i + 1}\n")
        f.write(f"l {n} 1\n")


def extract_case_contours(case_dir: Path) -> tuple[int, list[dict]]:
    """Return (n_slices_with_et, index rows) for one case."""
    case_id = case_dir.name
    seg_path = case_dir / f"{case_id}{SEG_SUFFIX}"
    if not seg_path.exists():
        return 0, []

    seg = nib.load(str(seg_path)).get_fdata().astype(np.int16)
    out_dir = CONTOURS_DIR / case_id

    rows: list[dict] = []
    n_with_et = 0

    n_slices = seg.shape[2]
    for z in range(n_slices):
        mask = seg[:, :, z] == ET_LABEL
        if mask.sum() < MIN_PIXELS:
            continue
        n_with_et += 1

        contours = measure.find_contours(mask.astype(np.float32), level=0.5)
        if not contours:
            continue

        # Largest contour by point count (the dominant outer boundary).
        contour_rc = max(contours, key=lambda c: c.shape[0])
        # find_contours returns (row, col); flip to (x, y) = (col, row).
        contour_xy = np.column_stack([contour_rc[:, 1], contour_rc[:, 0]])

        # find_contours often duplicates the first point at the end — drop it
        # so the explicit "l N 1" closing edge is not redundant.
        if (contour_xy.shape[0] >= 2
                and np.allclose(contour_xy[0], contour_xy[-1])):
            contour_xy = contour_xy[:-1]

        if contour_xy.shape[0] < 3:
            continue

        obj_path = out_dir / f"slice_{z:04d}.obj"
        write_contour_obj(obj_path, contour_xy, case_id, z)

        rows.append({
            "case_id": case_id,
            "slice_z": z,
            "n_vertices": contour_xy.shape[0],
            "obj_path": str(obj_path.relative_to(PROJECT_ROOT)),
        })

    return n_with_et, rows


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

def main() -> int:
    case_dirs = stage_raw_cases()
    CONTOURS_DIR.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    cases_with_et    = 0
    cases_without_et = 0

    for case_dir in tqdm(case_dirs, desc="contours"):
        n_with_et, rows = extract_case_contours(case_dir)
        if n_with_et > 0 and rows:
            cases_with_et += 1
            all_rows.extend(rows)
        else:
            cases_without_et += 1

    with open(INDEX_CSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["case_id", "slice_z", "n_vertices", "obj_path"]
        )
        writer.writeheader()
        writer.writerows(all_rows)

    total_processed = cases_with_et + cases_without_et
    print()
    print(f"OK Casos procesados:   {total_processed}/{N_CASES}")
    print(f"OK Contornos exportados: {len(all_rows)} archivos .obj")
    print(f"OK Index guardado en:    {INDEX_CSV.relative_to(PROJECT_ROOT)}")
    print(f"!! Casos sin ET:         {cases_without_et}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
