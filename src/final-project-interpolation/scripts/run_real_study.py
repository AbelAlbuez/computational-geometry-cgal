#!/usr/bin/env python3
# =============================================================================
# Real BraTS end-to-end study.
#   1. extract ET contours from N cases
#   2. pick the cases with the most ET slices
#   3. per case: Poisson-reconstruct (linear + spline), render 3D PNGs, compute
#      surface metrics, and compare mesh volume to the segmentation's voxel
#      ground-truth volume
#   4. run leave-one-slice-out per case and aggregate across cases
# Outputs everything to <out> (default data/results/, gitignored).
# =============================================================================
from __future__ import annotations

import argparse
import collections
import csv
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
METHODS = ["linear", "spline", "sdf"]


def run(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def run_ok(cmd):
    """Non-fatal run: returns (returncode, stdout, stderr)."""
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def grab(txt, key):
    for line in txt.splitlines():
        if key in line:
            return line.split(":", 1)[1].strip()
    return ""


def voxel_volume_mm3(seg_path: Path, label: int = 3) -> float:
    img = nib.load(str(seg_path))
    a = np.asanyarray(img.dataobj)
    zx, zy, zz = [float(z) for z in img.header.get_zooms()[:3]]
    return float((a == label).sum()) * zx * zy * zz


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=os.environ.get("BRATS_DIR"))
    ap.add_argument("--bin", required=True)
    ap.add_argument("--extract", type=int, default=25)
    ap.add_argument("--topk", type=int, default=3)
    ap.add_argument("--min-slices", type=int, default=15)
    ap.add_argument("--out", default=str(PROJ / "data" / "results"))
    args = ap.parse_args()
    if not args.dataset:
        sys.exit("Set --dataset or $BRATS_DIR")
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    print(f"== extracting {args.extract} cases ==")
    run([sys.executable, str(HERE / "extract_contours.py"),
         "--dataset", args.dataset, "--cases", str(args.extract)])

    contours = PROJ / "data" / "contours"
    rows = list(csv.DictReader(open(contours / "index.csv")))
    by = collections.defaultdict(list)
    for r in rows:
        by[r["case_id"]].append(int(r["slice_z"]))
    dz = float(rows[0].get("dz_mm", 1.0))
    cand = sorted(((c, len(z)) for c, z in by.items()), key=lambda x: -x[1])
    selected = [c for c, n in cand if n >= args.min_slices][:args.topk]
    print("selected cases:", selected)

    recon = []     # per case/method reconstruction summary
    loo_raw = []   # all raw LOO rows across cases
    ds = Path(args.dataset)

    for case in selected:
        cdir = str(contours / case)
        seg = ds / case / f"{case}-seg.nii.gz"
        gt_vol = voxel_volume_mm3(seg) if seg.exists() else float("nan")

        for method in ("linear", "spline"):
            mesh = str(out / f"mesh_{case}_{method}.off")
            rc, _, err = run_ok([args.bin, "reconstruct", "--contours", cdir,
                                 "--method", method, "--out", mesh,
                                 "--dz", str(dz), "--upsample", "3"])
            closed = (rc == 0)            # rc==4: mesh written but not closed
            mesh_exists = Path(mesh).exists()
            if not mesh_exists:
                print(f"  {case} {method}: reconstruction FAILED (rc={rc})")
                recon.append(dict(case=case, method=method, status=f"failed(rc={rc})",
                                  mesh_volume_mm3="", voxel_gt_volume_mm3=f"{gt_vol:.1f}",
                                  surface_area_mm2="", stack_volume="",
                                  symmetry="", mean_curvature=""))
                continue
            run([sys.executable, str(HERE / "render_mesh.py"), mesh,
                 str(out / f"mesh_{case}_{method}.png"), "--title", f"{case} {method}"])
            vol = sym = area = stack = curv = ""
            if closed:
                m = run([args.bin, "metrics", "--mesh", mesh, "--contours", cdir,
                         "--dz", str(dz), "--grid", "36"]).stdout
                vol = grab(m, "volume         "); area = grab(m, "surface area")
                stack = grab(m, "stack volume"); sym = grab(m, "symmetry")
                curv = grab(m, "mean curvature")
            recon.append(dict(
                case=case, method=method,
                status="closed" if closed else "open(not volume-bounding)",
                mesh_volume_mm3=vol, voxel_gt_volume_mm3=f"{gt_vol:.1f}",
                surface_area_mm2=area, stack_volume=stack,
                symmetry=sym, mean_curvature=curv))
            print(f"  {case} {method}: {'closed' if closed else 'OPEN'}  "
                  f"mesh_vol={vol or 'n/a'}  voxel_GT={gt_vol:.0f} mm^3")

        # leave-one-slice-out for this case
        od = out / f"loo_{case}"
        run([sys.executable, str(HERE / "evaluate.py"), "--contours", cdir,
             "--bin", args.bin, "--out", str(od)])
        for rr in csv.DictReader(open(od / "results_raw.csv")):
            rr["case"] = case
            loo_raw.append(rr)
        print(f"  {case}: LOO done")

    # -- write reconstruction summary
    with open(out / "reconstruction_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(recon[0].keys()))
        w.writeheader(); w.writerows(recon)

    # -- aggregate LOO across all cases/slices, per method
    agg = []
    for method in METHODS:
        mr = [r for r in loo_raw if r["method"] == method]
        if not mr:
            continue
        f = lambda k: np.mean([float(r[k]) for r in mr])
        agg.append(dict(method=method, n=len(mr), dice=f("dice"), iou=f("iou"),
                        hausdorff=f("hausdorff"), mean_dist=f("mean_dist"),
                        rel_area_err=f("rel_area_err")))
    with open(out / "loo_aggregate.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["method", "n", "dice", "iou",
                                          "hausdorff", "mean_dist", "rel_area_err"])
        w.writeheader(); w.writerows(agg)

    # -- aggregate comparison plot
    if agg:
        labels = [a["method"] for a in agg]
        fig, axs = plt.subplots(1, 3, figsize=(12, 4))
        axs[0].bar(labels, [a["dice"] for a in agg], color="tab:blue")
        axs[0].set_title("Dice (higher better)"); axs[0].set_ylim(0, 1)
        axs[1].bar(labels, [a["hausdorff"] for a in agg], color="tab:orange")
        axs[1].set_title("Hausdorff mm (lower better)")
        axs[2].bar(labels, [a["rel_area_err"] for a in agg], color="tab:green")
        axs[2].set_title("Rel. area error (lower better)")
        fig.suptitle(f"Leave-one-slice-out on real BraTS ({len(selected)} cases)")
        fig.tight_layout(); fig.savefig(out / "loo_aggregate.png", dpi=120)
        plt.close(fig)

    print("\n=== Aggregated LOO (real BraTS) ===")
    print(f"{'method':<8}{'n':>5}{'Dice':>8}{'IoU':>8}{'Hausd':>9}{'meanD':>8}{'areaErr':>9}")
    for a in agg:
        print(f"{a['method']:<8}{a['n']:>5}{a['dice']:>8.3f}{a['iou']:>8.3f}"
              f"{a['hausdorff']:>9.3f}{a['mean_dist']:>8.3f}{a['rel_area_err']:>9.3f}")
    print(f"\nresults in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
