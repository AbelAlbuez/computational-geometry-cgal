#!/usr/bin/env python3
# =============================================================================
# Build ONE interactive HTML dashboard from a results folder:
#   - per visualisation case, a row of rotatable 3D surfaces:
#       ground truth (marching cubes) | linear | spline | sdf
#   - the leave-one-slice-out stats table (mean / median / %Dice>0.8)
#   - the reconstruction-vs-ground-truth volume table
# Self-contained (plotly.js embedded); open in any browser.
#
#   python3 build_dashboard.py --results data/results --out dashboard.html
# =============================================================================
from __future__ import annotations

import argparse
import csv
import glob
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

COLS = ["gt", "linear", "spline", "sdf"]
COL_TITLE = {"gt": "ground truth", "linear": "linear",
             "spline": "spline", "sdf": "sdf"}
COL_COLOR = {"gt": "#8fce8f", "linear": "#9ec5fe",
             "spline": "#f5c2a0", "sdf": "#c5a3e0"}


def read_off(path):
    with open(path) as f:
        toks = f.read().split()
    p = 1 if toks[0].upper() == "OFF" else 0
    nv, nf = int(toks[p]), int(toks[p + 1])
    vals = toks[p + 3:]
    verts = np.array(vals[:nv * 3], dtype=float).reshape(nv, 3)
    rest = vals[nv * 3:]
    faces, q = [], 0
    for _ in range(nf):
        kk = int(rest[q]); ids = list(map(int, rest[q + 1:q + 1 + kk])); q += 1 + kk
        for t in range(1, kk - 1):
            faces.append((ids[0], ids[t], ids[t + 1]))
    return verts, np.asarray(faces, dtype=int)


def table_html(rows, headers, title):
    th = "".join(f"<th>{h}</th>" for h in headers)
    trs = ""
    for r in rows:
        trs += "<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>"
    return (f"<h3>{title}</h3><table>"
            f"<thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table>")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    res = Path(args.results)

    cases = sorted(os.path.basename(p)[5:-7]
                   for p in glob.glob(str(res / "mesh_*_gt.off")))
    if not cases:
        cases = sorted({os.path.basename(p).split("_")[1]
                        for p in glob.glob(str(res / "mesh_*_linear.off"))})
    print("dashboard cases:", cases)

    nrows = max(1, len(cases))
    titles = [f"{c} — {COL_TITLE[m]}" for c in cases for m in COLS]
    fig = make_subplots(rows=nrows, cols=4,
                        specs=[[{"type": "scene"}] * 4 for _ in range(nrows)],
                        subplot_titles=titles, horizontal_spacing=0.01,
                        vertical_spacing=0.04)
    for ri, case in enumerate(cases, start=1):
        for ci, m in enumerate(COLS, start=1):
            off = res / f"mesh_{case}_{m}.off"
            if not off.exists():
                continue
            v, f = read_off(str(off))
            if len(f) == 0:
                continue
            fig.add_trace(go.Mesh3d(
                x=v[:, 0], y=v[:, 1], z=v[:, 2],
                i=f[:, 0], j=f[:, 1], k=f[:, 2],
                color=COL_COLOR[m], opacity=1.0,
                lighting=dict(ambient=0.5, diffuse=0.85),
                showscale=False, hoverinfo="skip"),
                row=ri, col=ci)
    fig.update_scenes(aspectmode="data", xaxis_visible=False,
                      yaxis_visible=False, zaxis_visible=False)
    fig.update_layout(height=330 * nrows, margin=dict(l=0, r=0, t=40, b=0),
                      title_text="BraTS tumor reconstruction — ground truth vs "
                                 "linear / spline / SDF interpolation")
    plot_div = fig.to_html(full_html=False, include_plotlyjs=True)

    # tables
    tables = ""
    agg = res / "loo_aggregate.csv"
    if agg.exists():
        rd = list(csv.DictReader(open(agg)))
        headers = ["method", "n", "Dice mean", "Dice median", "% Dice>0.8",
                   "IoU mean", "Hausdorff mean (mm)", "area-err median"]
        rows = [[r["method"], r["n"], f"{float(r['dice_mean']):.3f}",
                 f"{float(r['dice_median']):.3f}", f"{float(r['dice_pct_gt80']):.1f}%",
                 f"{float(r['iou_mean']):.3f}", f"{float(r['hausdorff_mean']):.2f}",
                 f"{float(r['area_err_median']):.3f}"] for r in rd]
        tables += table_html(rows, headers,
                             "Leave-one-slice-out accuracy (held-out real slices)")
    rec = res / "reconstruction_summary.csv"
    if rec.exists():
        rd = list(csv.DictReader(open(rec)))
        headers = ["case", "method", "status", "mesh vol (mm³)", "voxel GT (mm³)", "ratio"]
        rows = []
        for r in rd:
            mv = r["mesh_volume_mm3"].replace(" mm^3", "").strip()
            gt = r["voxel_gt_volume_mm3"].strip()
            try:
                ratio = f"{float(mv) / float(gt):.2f}"
            except Exception:
                ratio = "—"
            rows.append([r["case"], r["method"], r.get("status", ""), mv or "—", gt, ratio])
        tables += table_html(rows, headers, "Reconstruction volume vs voxel ground truth")

    style = """<style>
      body{font-family:system-ui,Arial,sans-serif;margin:18px;color:#222}
      h1{font-size:20px} h3{margin-top:26px}
      table{border-collapse:collapse;margin:8px 0 20px}
      th,td{border:1px solid #ccc;padding:4px 10px;text-align:right;font-size:13px}
      th{background:#f0f0f0} td:first-child,th:first-child{text-align:left}
      .note{color:#555;font-size:13px;max-width:900px}
    </style>"""
    note = ("<p class='note'>Each 3D panel is interactive: drag to rotate, scroll "
            "to zoom. Ground truth is a marching-cubes surface of the segmented "
            "voxels; the others are Poisson reconstructions of the interpolated "
            "contour stack. Dice/IoU are region overlap (higher better); Hausdorff "
            "and area error are lower-better. Median &gt;&gt; mean because a few "
            "hard slices (tiny or splitting tumor) score near zero.</p>")
    page = (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>BraTS interpolation dashboard</title>{style}</head><body>"
            f"<h1>BraTS glioma — interpolation &amp; reconstruction dashboard</h1>"
            f"{note}{plot_div}{tables}</body></html>")
    Path(args.out).write_text(page, encoding="utf-8")
    print(f"wrote {args.out}  ({len(cases)} cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
