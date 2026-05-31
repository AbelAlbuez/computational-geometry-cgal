#!/usr/bin/env python3
# =============================================================================
# Real BraTS end-to-end study.
#   1. extract ET contours from many cases; pick those with the most ET slices
#   2. leave-one-slice-out on ALL selected cases (linear/spline/sdf) -> stats
#      (mean, median, std, % Dice>0.8) aggregated across cases
#   3. for a few "visualisation" cases: Poisson-reconstruct linear/spline/sdf,
#      build a marching-cubes GROUND-TRUTH surface from the voxels, render PNGs,
#      and compare mesh volumes to the voxel ground-truth volume
# Outputs to <out> (default data/results/, gitignored). A separate dashboard
# script turns the meshes + stats into one interactive HTML.
# =============================================================================
from __future__ import annotations

import argparse
import collections
import csv
import glob
import os
import shutil
import statistics as st
import subprocess
import sys
from pathlib import Path

import numpy as np
import nibabel as nib
from skimage import measure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
METHODS = ["linear", "spline", "sdf"]


def run(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def run_ok(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def grab(txt, key):
    for line in txt.splitlines():
        if key in line:
            return line.split(":", 1)[1].strip()
    return ""


def voxel_volume_mm3(seg_path, label=3):
    img = nib.load(str(seg_path))
    a = np.asanyarray(img.dataobj)
    zx, zy, zz = [float(z) for z in img.header.get_zooms()[:3]]
    return float((a == label).sum()) * zx * zy * zz


def write_off(path, verts, faces):
    with open(path, "w") as f:
        f.write("OFF\n")
        f.write(f"{len(verts)} {len(faces)} 0\n")
        for v in verts:
            f.write(f"{v[0]:.5f} {v[1]:.5f} {v[2]:.5f}\n")
        for t in faces:
            f.write(f"3 {t[0]} {t[1]} {t[2]}\n")


def gt_surface(seg_path, off_out, label=3, step=2):
    """Marching-cubes ground-truth surface from the ET voxel mask (mm units),
    reordered to (x=col, y=row, z=slice) to match our reconstructions."""
    img = nib.load(str(seg_path))
    a = np.asanyarray(img.dataobj)
    zx, zy, zz = [float(z) for z in img.header.get_zooms()[:3]]
    mask = (a == label)
    if mask.sum() < 10:
        return False
    verts, faces, _, _ = measure.marching_cubes(
        mask.astype(np.float32), level=0.5, spacing=(zx, zy, zz), step_size=step)
    verts = verts[:, [1, 0, 2]]            # (axis0,axis1,axis2) -> (x,y,z)
    write_off(off_out, verts, faces)
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=os.environ.get("BRATS_DIR"))
    ap.add_argument("--bin", required=True)
    ap.add_argument("--extract", type=int, default=40)
    ap.add_argument("--topk", type=int, default=10, help="cases for LOO/stats")
    ap.add_argument("--n-vis", type=int, default=3, help="cases for 3D + GT")
    ap.add_argument("--min-slices", type=int, default=15)
    ap.add_argument("--out", default=str(PROJ / "data" / "results"))
    args = ap.parse_args()
    if not args.dataset:
        sys.exit("Set --dataset or $BRATS_DIR")
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    ds = Path(args.dataset)

    # Clean stale per-run artefacts so a smaller rerun does not leave old cases.
    for p in glob.glob(str(out / "mesh_*")) + glob.glob(str(out / "loo_*")):
        shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)

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
    vis = selected[:args.n_vis]
    print("LOO/stats cases:", selected)
    print("3D/GT cases    :", vis)

    # ---- 3D reconstruction + GT for the visualisation cases ----------------
    recon = []
    for case in vis:
        cdir = str(contours / case)
        seg = ds / case / f"{case}-seg.nii.gz"
        gt_vol = voxel_volume_mm3(seg) if seg.exists() else float("nan")
        if seg.exists() and gt_surface(seg, str(out / f"mesh_{case}_gt.off")):
            run([sys.executable, str(HERE / "render_mesh.py"),
                 str(out / f"mesh_{case}_gt.off"),
                 str(out / f"mesh_{case}_gt.png"), "--title", f"{case} ground truth"])
        for method in METHODS:
            mesh = str(out / f"mesh_{case}_{method}.off")
            rc, _, _ = run_ok([args.bin, "reconstruct", "--contours", cdir,
                               "--method", method, "--out", mesh,
                               "--dz", str(dz), "--upsample", "3"])
            if not Path(mesh).exists():
                recon.append(dict(case=case, method=method, status=f"failed(rc={rc})",
                                  mesh_volume_mm3="", voxel_gt_volume_mm3=f"{gt_vol:.1f}"))
                print(f"  {case} {method}: FAILED rc={rc}"); continue
            run([sys.executable, str(HERE / "render_mesh.py"), mesh,
                 str(out / f"mesh_{case}_{method}.png"), "--title", f"{case} {method}"])
            vol = ""
            if rc == 0:
                m = run([args.bin, "metrics", "--mesh", mesh, "--contours", cdir,
                         "--dz", str(dz), "--grid", "36"]).stdout
                vol = grab(m, "volume         ")
            recon.append(dict(case=case, method=method,
                              status="closed" if rc == 0 else "open",
                              mesh_volume_mm3=vol, voxel_gt_volume_mm3=f"{gt_vol:.1f}"))
            print(f"  {case} {method}: {'closed' if rc==0 else 'OPEN'} vol={vol or 'n/a'}")

    with open(out / "reconstruction_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case", "method", "status",
                                          "mesh_volume_mm3", "voxel_gt_volume_mm3"])
        w.writeheader(); w.writerows(recon)

    # ---- leave-one-slice-out on all selected cases -------------------------
    loo_raw = []
    for case in selected:
        cdir = str(contours / case)
        od = out / f"loo_{case}"
        run([sys.executable, str(HERE / "evaluate.py"), "--contours", cdir,
             "--bin", args.bin, "--out", str(od)])
        for rr in csv.DictReader(open(od / "results_raw.csv")):
            rr["case"] = case; loo_raw.append(rr)
        print(f"  {case}: LOO done")

    # ---- aggregate stats (mean / median / std / %Dice>0.8) -----------------
    agg = []
    for method in METHODS:
        mr = [r for r in loo_raw if r["method"] == method]
        if not mr:
            continue
        col = lambda k: [float(r[k]) for r in mr]
        d = col("dice")
        agg.append(dict(
            method=method, n=len(mr),
            dice_mean=st.mean(d), dice_median=st.median(d),
            dice_pct_gt80=100.0 * sum(x > 0.8 for x in d) / len(d),
            iou_mean=st.mean(col("iou")),
            hausdorff_mean=st.mean(col("hausdorff")),
            hausdorff_median=st.median(col("hausdorff")),
            area_err_mean=st.mean(col("rel_area_err")),
            area_err_median=st.median(col("rel_area_err"))))
    with open(out / "loo_aggregate.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(agg[0].keys()))
        w.writeheader(); w.writerows(agg)

    if agg:
        labels = [a["method"] for a in agg]
        fig, axs = plt.subplots(1, 3, figsize=(12, 4))
        x = np.arange(len(labels))
        axs[0].bar(x - 0.2, [a["dice_mean"] for a in agg], 0.4, label="mean", color="tab:blue")
        axs[0].bar(x + 0.2, [a["dice_median"] for a in agg], 0.4, label="median", color="tab:cyan")
        axs[0].set_xticks(x); axs[0].set_xticklabels(labels); axs[0].set_ylim(0, 1)
        axs[0].set_title("Dice"); axs[0].legend()
        axs[1].bar(labels, [a["hausdorff_mean"] for a in agg], color="tab:orange")
        axs[1].set_title("Hausdorff mean (mm, lower better)")
        axs[2].bar(labels, [a["dice_pct_gt80"] for a in agg], color="tab:green")
        axs[2].set_title("% slices with Dice > 0.8"); axs[2].set_ylim(0, 100)
        fig.suptitle(f"Leave-one-slice-out on real BraTS ({len(selected)} cases, "
                     f"{agg[0]['n']} slices/method)")
        fig.tight_layout(); fig.savefig(out / "loo_aggregate.png", dpi=120); plt.close(fig)

    print(f"\n=== Aggregated LOO (real BraTS, {len(selected)} cases) ===")
    print(f"{'method':<8}{'n':>5}{'Dice mean':>11}{'Dice med':>10}{'%>0.8':>8}"
          f"{'Hausd mean':>12}{'areaErr med':>13}")
    for a in agg:
        print(f"{a['method']:<8}{a['n']:>5}{a['dice_mean']:>11.3f}{a['dice_median']:>10.3f}"
              f"{a['dice_pct_gt80']:>7.1f}%{a['hausdorff_mean']:>12.3f}{a['area_err_median']:>13.3f}")
    print(f"\nresults in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
