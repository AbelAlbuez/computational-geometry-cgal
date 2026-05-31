#!/usr/bin/env python3
# =============================================================================
# Build ONE interactive HTML dashboard with a CASE SELECTOR (dropdown):
#   - pick a tumor from the dropdown; only that case is shown.
#   - per case: rotatable 3D surfaces  ground truth | linear | spline | sdf
#     plus a metrics table (LOO accuracy + reconstruction volume vs ground truth).
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


def load_loo(res, case):
    """per-method LOO row (means for this case) from loo_<case>/results.csv."""
    p = res / f"loo_{case}" / "results.csv"
    out = {}
    if p.exists():
        for r in csv.DictReader(open(p)):
            out[r["method"]] = r
    return out


def load_recon(res):
    out = {}
    p = res / "reconstruction_summary.csv"
    if p.exists():
        for r in csv.DictReader(open(p)):
            out[(r["case"], r["method"])] = r
    return out


def case_table(case, loo, recon):
    methods = ["linear", "spline", "sdf"]
    def vol(m):
        r = recon.get((case, m), {})
        return r.get("mesh_volume_mm3", "").replace(" mm^3", "").strip()
    gt = ""
    for m in methods:
        r = recon.get((case, m), {})
        if r.get("voxel_gt_volume_mm3"):
            gt = r["voxel_gt_volume_mm3"].strip(); break
    def f(m, k, nd=3):
        try:
            return f"{float(loo.get(m, {}).get(k, '')):.{nd}f}"
        except Exception:
            return "—"
    def ratio(m):
        try:
            return f"{float(vol(m)) / float(gt):.2f}"
        except Exception:
            return "—"
    header = ["method", "Dice", "IoU", "Hausdorff (mm)", "area err",
              "mesh vol (mm³)", "GT vol (mm³)", "ratio"]
    cols = [
        methods,
        [f(m, "dice") for m in methods],
        [f(m, "iou") for m in methods],
        [f(m, "hausdorff", 2) for m in methods],
        [f(m, "rel_area_err") for m in methods],
        [vol(m) or "—" for m in methods],
        [gt or "—"] * 3,
        [ratio(m) for m in methods],
    ]
    return header, cols


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
    if not cases:
        print("no meshes found in", res); return 1
    print("dashboard cases:", cases)
    recon = load_recon(res)

    fig = make_subplots(
        rows=2, cols=4,
        specs=[[{"type": "scene"}] * 4,
               [{"type": "table", "colspan": 4}, None, None, None]],
        row_heights=[0.74, 0.26], vertical_spacing=0.06,
        horizontal_spacing=0.005,
        subplot_titles=[COL_TITLE[m] for m in COLS] + [""])

    trace_case = []          # case index that each trace belongs to
    for ci, case in enumerate(cases):
        for col, m in enumerate(COLS, start=1):
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
                showscale=False, hoverinfo="skip",
                visible=(ci == 0)), row=1, col=col)
            trace_case.append(ci)
        header, cols = case_table(case, load_loo(res, case), recon)
        fig.add_trace(go.Table(
            header=dict(values=header, fill_color="#e8e8e8",
                        align="left", font=dict(size=12)),
            cells=dict(values=cols, align="left", font=dict(size=12),
                       height=24), visible=(ci == 0)), row=2, col=1)
        trace_case.append(ci)

    ntr = len(trace_case)
    buttons = []
    for ci, case in enumerate(cases):
        vis = [trace_case[t] == ci for t in range(ntr)]
        buttons.append(dict(label=case, method="update",
                            args=[{"visible": vis},
                                  {"title.text": f"Tumor {case}  —  ground truth vs "
                                                 f"linear / spline / SDF"}]))

    fig.update_scenes(aspectmode="data", xaxis_visible=False,
                      yaxis_visible=False, zaxis_visible=False)
    fig.update_layout(
        height=860, margin=dict(l=0, r=0, t=110, b=0),
        title=dict(text=f"Tumor {cases[0]}  —  ground truth vs linear / spline / SDF",
                   x=0.5, xanchor="center", y=0.99),
        updatemenus=[dict(
            buttons=buttons, direction="down", showactive=True,
            x=0.5, xanchor="center", y=1.10, yanchor="top",
            pad=dict(t=4, b=4), bgcolor="#f4f4f4")],
        annotations=list(fig.layout.annotations) + [dict(
            text="Select tumor:", x=0.5, xanchor="right", y=1.135,
            yref="paper", xref="paper", showarrow=False,
            font=dict(size=13), xshift=-120)])

    fig.write_html(args.out, include_plotlyjs=True, full_html=True)
    size = os.path.getsize(args.out) / 1e6
    print(f"wrote {args.out}  ({len(cases)} cases, {size:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
