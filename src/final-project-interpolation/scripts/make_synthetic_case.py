#!/usr/bin/env python3
# =============================================================================
# Generate a SYNTHETIC BraTS-format dataset so the whole pipeline can be run
# end-to-end without the real (large, restricted) BraTS 2024 GLI download.
#
# It writes, for each case:
#   <out>/BraTS2024-BraTS-GLI-TrainingData/<case>/<case>-seg.nii.gz
# containing an ET (label 3) blob shaped like a lobed ellipsoid, with the
# NIfTI affine carrying the requested voxel spacing (so the spacing path in
# extract_contours.py is exercised). Point extract_contours.py at <out> via
# --dataset or $BRATS_DIR. This is the ONLY substitute for real BraTS data;
# everything downstream is identical.
# =============================================================================
from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np

CASE_PREFIX = "BraTS-GLI-"
EXTRACTED   = "BraTS2024-BraTS-GLI-TrainingData"
ET_LABEL    = 3


def make_seg(shape, a, c, lobes, phase) -> np.ndarray:
    """A lobed ellipsoid of ET voxels centred in the volume."""
    nx, ny, nz = shape
    seg = np.zeros(shape, dtype=np.uint8)
    cx, cy, cz = nx / 2.0, ny / 2.0, nz / 2.0
    zz = np.arange(nz)
    for z in zz:
        u = (z - cz) / c
        if abs(u) >= 1.0:
            continue
        base = a * np.sqrt(1.0 - u * u)
        if base < 2.0:
            continue
        yy, xx = np.meshgrid(np.arange(ny), np.arange(nx))
        ang = np.arctan2(yy - cy, xx - cx)
        rad = base * (1.0 + 0.15 * np.cos(lobes * ang + phase))
        dist = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
        seg[:, :, z][dist <= rad] = ET_LABEL
    return seg


def main() -> int:
    p = argparse.ArgumentParser(description="Synthetic BraTS-format dataset.")
    p.add_argument("--out", required=True, help="output dataset root")
    p.add_argument("--cases", type=int, default=2)
    p.add_argument("--spacing", type=float, nargs=3, default=[1.0, 1.0, 1.0],
                   metavar=("DX", "DY", "DZ"), help="voxel spacing in mm")
    p.add_argument("--size", type=int, nargs=3, default=[80, 80, 70],
                   metavar=("NX", "NY", "NZ"))
    args = p.parse_args()

    root = Path(args.out) / EXTRACTED
    affine = np.diag(args.spacing + [1.0])
    for i in range(args.cases):
        case = f"{CASE_PREFIX}{i:05d}-000"
        cdir = root / case
        cdir.mkdir(parents=True, exist_ok=True)
        seg = make_seg(tuple(args.size),
                       a=20.0 + 4 * i, c=26.0,
                       lobes=3 + i, phase=0.4 * i)
        nib.save(nib.Nifti1Image(seg, affine), str(cdir / f"{case}-seg.nii.gz"))
        print(f"wrote {case}  ET voxels={int((seg == ET_LABEL).sum())}")
    print(f"dataset root: {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
