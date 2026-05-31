#!/usr/bin/env python3
# =============================================================================
# Final project: leave-one-slice-out (LOO) comparison of interpolation methods.
#
# For every interior slice k of a case, hide C_k and predict it from its
# neighbours C_{k-1}, C_{k+1} with each method (linear / sdf at t=0.5; spline
# over a small z-window), then score the prediction against the true C_k with:
#   - Dice / IoU of the enclosed regions (rasterized),
#   - symmetric Hausdorff distance and mean contour distance,
#   - relative area error.
# Aggregates per method -> results.csv + comparison bar charts.
#
# Predictions are produced by the C++ `contour_interpolator` binary so the
# study uses exactly the same geometry as the rest of the pipeline.
# =============================================================================
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

METHODS = ["linear", "spline", "sdf"]


# ----------------------------------------------------------------------------
# .obj contour I/O and polygon geometry
# ----------------------------------------------------------------------------
def read_obj(path: str) -> np.ndarray:
    pts = []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                _, x, y = line.split()[:3]
                pts.append((float(x), float(y)))
    return np.asarray(pts, dtype=float)


def shoelace_area(p: np.ndarray) -> float:
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def resample_closed(p: np.ndarray, n: int = 300) -> np.ndarray:
    d = np.sqrt(((np.roll(p, -1, 0) - p) ** 2).sum(1))
    s = np.concatenate([[0], np.cumsum(d)])
    total = s[-1]
    if total <= 0:
        return p
    targets = np.linspace(0, total, n, endpoint=False)
    out = np.empty((n, 2))
    j = 0
    for i, t in enumerate(targets):
        while j + 1 < len(s) and s[j + 1] < t:
            j += 1
        seg = s[j + 1] - s[j] if j + 1 < len(s) else 1.0
        f = (t - s[j]) / seg if seg > 0 else 0.0
        a = p[j % len(p)]
        b = p[(j + 1) % len(p)]
        out[i] = a + f * (b - a)
    return out


def dice_iou(a: np.ndarray, b: np.ndarray):
    lo = np.minimum(a.min(0), b.min(0)) - 2
    hi = np.maximum(a.max(0), b.max(0)) + 2
    ext = (hi - lo).max()
    h = max(0.25, ext / 200.0)
    xs = np.arange(lo[0], hi[0], h)
    ys = np.arange(lo[1], hi[1], h)
    gx, gy = np.meshgrid(xs, ys)
    grid = np.column_stack([gx.ravel(), gy.ravel()])
    ina = MplPath(a).contains_points(grid)
    inb = MplPath(b).contains_points(grid)
    inter = np.logical_and(ina, inb).sum()
    union = np.logical_or(ina, inb).sum()
    sa, sb = ina.sum(), inb.sum()
    dice = (2.0 * inter / (sa + sb)) if (sa + sb) > 0 else 0.0
    iou = (inter / union) if union > 0 else 0.0
    return dice, iou


def hausdorff_and_mean(a: np.ndarray, b: np.ndarray):
    A = resample_closed(a, 300)
    B = resample_closed(b, 300)
    d = np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))
    ab = d.min(1)
    ba = d.min(0)
    hausdorff = max(ab.max(), ba.max())
    mean_d = 0.5 * (ab.mean() + ba.mean())
    return hausdorff, mean_d


# ----------------------------------------------------------------------------
# slice discovery + prediction via the C++ binary
# ----------------------------------------------------------------------------
def stem_index(name: str) -> int:
    m = re.search(r"(\d+)\D*$", Path(name).stem)
    return int(m.group(1)) if m else -1


def list_slices(contours: str):
    files = [str(p) for p in Path(contours).glob("*.obj")]
    files.sort(key=lambda f: (stem_index(f), f))
    return [(stem_index(f), f) for f in files]


def predict_pairwise(binp, method, fa, fb, out, t=0.5):
    subprocess.run([binp, method, fa, fb, f"{t:.6f}", out],
                   check=True, capture_output=True)


def predict_spline(binp, slices, k, out_dir):
    # window of up to 2 neighbours each side, excluding k
    idxs = [j for j in (k - 2, k - 1, k + 1, k + 2) if 0 <= j < len(slices)]
    zs = [slices[j][0] for j in idxs]
    with tempfile.TemporaryDirectory() as win:
        for j in idxs:
            z, f = slices[j]
            shutil.copy(f, os.path.join(win, f"slice_{z:05d}.obj"))
        K = 2
        subprocess.run([binp, "series", "--kind", "spline", "--window", win,
                        "--upsample", str(K), out_dir], check=True, capture_output=True)
        # rebuild the query list to find the file at z_k
        q = []
        for a in range(len(zs) - 1):
            for jj in range(K):
                q.append(zs[a] + (zs[a + 1] - zs[a]) * (jj / K))
        q.append(zs[-1])
        zk = slices[k][0]
        qi = int(np.argmin([abs(v - zk) for v in q]))
        return os.path.join(out_dir, f"interp_{qi:04d}.obj")


# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contours", required=True, help="dir of slice_*.obj (one case)")
    ap.add_argument("--bin", required=True, help="path to contour_interpolator")
    ap.add_argument("--out", default="loo_out", help="output directory")
    ap.add_argument("--max-gap", type=int, default=3,
                    help="max z-gap (in slices) of a usable LOO triplet")
    args = ap.parse_args()

    slices = list_slices(args.contours)
    if len(slices) < 3:
        print("Need >= 3 slices for leave-one-out.", file=sys.stderr)
        return 2
    os.makedirs(args.out, exist_ok=True)

    rows = []                         # per (method, slice)
    skipped = 0
    with tempfile.TemporaryDirectory() as tmp:
        for k in range(1, len(slices) - 1):
            z0, z1, z2 = slices[k - 1][0], slices[k][0], slices[k + 1][0]
            # Real data has z-gaps: only use triplets whose neighbours straddle
            # the held-out slice within MAX_GAP, and set t from the true heights.
            if z2 - z0 <= 0 or z2 - z0 > args.max_gap or not (z0 < z1 < z2):
                skipped += 1
                continue
            t = (z1 - z0) / (z2 - z0)
            zk, true_f = slices[k]
            true = read_obj(true_f)
            if true.shape[0] < 3:
                continue
            fa, fb = slices[k - 1][1], slices[k + 1][1]
            for method in METHODS:
                try:
                    if method in ("linear", "sdf"):
                        pred_f = os.path.join(tmp, f"{method}_{k}.obj")
                        predict_pairwise(args.bin, method, fa, fb, pred_f, t)
                    else:
                        od = os.path.join(tmp, f"spline_{k}")
                        pred_f = predict_spline(args.bin, slices, k, od)
                    pred = read_obj(pred_f)
                    if pred.shape[0] < 3:
                        continue
                    dice, iou = dice_iou(pred, true)
                    hd, md = hausdorff_and_mean(pred, true)
                    a_true = shoelace_area(true)
                    rel_area = abs(shoelace_area(pred) - a_true) / a_true if a_true > 0 else 0
                    rows.append(dict(method=method, slice=zk, dice=dice, iou=iou,
                                     hausdorff=hd, mean_dist=md, rel_area_err=rel_area))
                except subprocess.CalledProcessError as e:
                    print(f"  {method} k={k} failed: {e}", file=sys.stderr)

    # aggregate per method
    agg = []
    for method in METHODS:
        mr = [r for r in rows if r["method"] == method]
        if not mr:
            continue
        agg.append(dict(
            method=method, n=len(mr),
            dice=np.mean([r["dice"] for r in mr]),
            iou=np.mean([r["iou"] for r in mr]),
            hausdorff=np.mean([r["hausdorff"] for r in mr]),
            mean_dist=np.mean([r["mean_dist"] for r in mr]),
            rel_area_err=np.mean([r["rel_area_err"] for r in mr]),
        ))

    csv_path = os.path.join(args.out, "results.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["method", "n", "dice", "iou",
                                          "hausdorff", "mean_dist", "rel_area_err"])
        w.writeheader()
        for a in agg:
            w.writerow(a)

    # raw per-slice rows (so several cases can be aggregated downstream)
    raw_path = os.path.join(args.out, "results_raw.csv")
    with open(raw_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["method", "slice", "dice", "iou",
                                          "hausdorff", "mean_dist", "rel_area_err"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # comparison plots
    if agg:
        labels = [a["method"] for a in agg]
        fig, axs = plt.subplots(1, 3, figsize=(12, 4))
        axs[0].bar(labels, [a["dice"] for a in agg], color="tab:blue");  axs[0].set_title("Dice (higher better)"); axs[0].set_ylim(0, 1)
        axs[1].bar(labels, [a["hausdorff"] for a in agg], color="tab:orange"); axs[1].set_title("Hausdorff (lower better)")
        axs[2].bar(labels, [a["rel_area_err"] for a in agg], color="tab:green"); axs[2].set_title("Rel. area error (lower better)")
        fig.suptitle("Leave-one-slice-out comparison")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "comparison.png"), dpi=120)

    # console table
    print(f"\n=== Leave-one-slice-out comparison (mean over slices; "
          f"{skipped} triplets skipped for z-gap > {args.max_gap}) ===")
    print(f"{'method':<8} {'n':>4} {'Dice':>7} {'IoU':>7} {'Hausd':>8} {'meanD':>7} {'areaErr':>8}")
    for a in agg:
        print(f"{a['method']:<8} {a['n']:>4} {a['dice']:>7.3f} {a['iou']:>7.3f} "
              f"{a['hausdorff']:>8.3f} {a['mean_dist']:>7.3f} {a['rel_area_err']:>8.3f}")
    print(f"\nwrote {csv_path} and comparison.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
