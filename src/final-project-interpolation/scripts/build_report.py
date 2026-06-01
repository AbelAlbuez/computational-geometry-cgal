#!/usr/bin/env python3
# =============================================================================
# Build a single self-contained HTML REPORT from a results folder:
#   - description + math of the three interpolation methods (MathJax)
#   - interactive BAR plots (per-method Dice/IoU, Hausdorff, % Dice>0.8)
#   - interactive LINE plots (Dice CDF per method; per-case mean Dice)
#   - a tumor SELECTOR dropdown: GT vs linear/spline/sdf 3D + per-case table
#   - aggregate + reconstruction tables
# plotly.js is embedded once; MathJax is loaded from CDN.
#
#   python3 build_report.py --results data/results --out report.html
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
from plotly.offline import get_plotlyjs

COLS = ["gt", "linear", "spline", "sdf"]
COL_TITLE = {"gt": "ground truth", "linear": "linear", "spline": "spline", "sdf": "sdf"}
COL_COLOR = {"gt": "#8fce8f", "linear": "#9ec5fe", "spline": "#f5c2a0", "sdf": "#c5a3e0"}
M_COLOR = {"linear": "#1f77b4", "spline": "#ff7f0e", "sdf": "#9467bd"}
METHODS = ["linear", "spline", "sdf"]


# ---------- data loading -----------------------------------------------------
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


def load_agg(res):
    p = res / "loo_aggregate.csv"
    return list(csv.DictReader(open(p))) if p.exists() else []


def load_raw(res):
    rows = []
    for fp in glob.glob(str(res / "loo_*" / "results_raw.csv")):
        case = os.path.basename(os.path.dirname(fp)).replace("loo_", "")
        for r in csv.DictReader(open(fp)):
            r["case"] = case; rows.append(r)
    return rows


def load_percase(res):
    out = {}
    for fp in glob.glob(str(res / "loo_*" / "results.csv")):
        case = os.path.basename(os.path.dirname(fp)).replace("loo_", "")
        out[case] = {r["method"]: r for r in csv.DictReader(open(fp))}
    return out


def load_recon(res):
    out = {}
    p = res / "reconstruction_summary.csv"
    if p.exists():
        for r in csv.DictReader(open(p)):
            out[(r["case"], r["method"])] = r
    return out


# ---------- plotly figures ---------------------------------------------------
def fig_bars(agg):
    fig = make_subplots(rows=1, cols=3, subplot_titles=(
        "Overlap (higher = better)", "Hausdorff mean (mm, lower = better)",
        "% slices Dice > 0.8"))
    methods = [a["method"] for a in agg]
    fig.add_trace(go.Bar(name="Dice mean", x=methods,
                         y=[float(a["dice_mean"]) for a in agg], marker_color="#1f77b4"), 1, 1)
    fig.add_trace(go.Bar(name="Dice median", x=methods,
                         y=[float(a["dice_median"]) for a in agg], marker_color="#7fb3e0"), 1, 1)
    fig.add_trace(go.Bar(name="IoU mean", x=methods,
                         y=[float(a["iou_mean"]) for a in agg], marker_color="#aec7e8"), 1, 1)
    fig.add_trace(go.Bar(showlegend=False, x=methods,
                         y=[float(a["hausdorff_mean"]) for a in agg],
                         marker_color="#ff7f0e"), 1, 2)
    fig.add_trace(go.Bar(showlegend=False, x=methods,
                         y=[float(a["dice_pct_gt80"]) for a in agg],
                         marker_color="#2ca02c"), 1, 3)
    fig.update_yaxes(range=[0, 1], row=1, col=1)
    fig.update_yaxes(range=[0, 100], row=1, col=3)
    fig.update_layout(barmode="group", height=380, margin=dict(l=10, r=10, t=40, b=10),
                      legend=dict(orientation="h", y=1.18, x=0.0))
    return fig


def fig_cdf(raw):
    fig = go.Figure()
    for m in METHODS:
        d = np.sort([float(r["dice"]) for r in raw if r["method"] == m])
        if len(d) == 0:
            continue
        y = np.arange(1, len(d) + 1) / len(d)
        fig.add_trace(go.Scatter(x=d, y=y, mode="lines", name=m,
                                 line=dict(color=M_COLOR[m], width=2)))
    fig.add_vline(x=0.8, line_dash="dash", line_color="gray",
                  annotation_text="Dice = 0.8")
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=30, b=10),
                      xaxis_title="Dice", yaxis_title="fraction of held-out slices ≤ Dice",
                      xaxis_range=[0, 1], yaxis_range=[0, 1],
                      legend=dict(x=0.02, y=0.98))
    return fig


def fig_percase(percase):
    fig = go.Figure()
    cases = list(percase.keys())
    base = sorted(cases, key=lambda c: -float(percase[c].get("linear", {}).get("dice", 0)))
    x = list(range(1, len(base) + 1))
    for m in METHODS:
        y = [float(percase[c].get(m, {}).get("dice", "nan")) for c in base]
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines+markers", name=m,
                                 line=dict(color=M_COLOR[m], width=1.5),
                                 marker=dict(size=4)))
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=30, b=10),
                      xaxis_title="case (sorted by linear Dice)",
                      yaxis_title="per-case mean Dice", yaxis_range=[0, 1],
                      legend=dict(x=0.02, y=0.05))
    return fig


def case_table_cols(case, percase, recon):
    def vol(m):
        return recon.get((case, m), {}).get("mesh_volume_mm3", "").replace(" mm^3", "").strip()
    gt = ""
    for m in METHODS:
        if recon.get((case, m), {}).get("voxel_gt_volume_mm3"):
            gt = recon[(case, m)]["voxel_gt_volume_mm3"].strip(); break
    def f(m, k, nd=3):
        try:
            return f"{float(percase.get(case, {}).get(m, {}).get(k, '')):.{nd}f}"
        except Exception:
            return "—"
    def ratio(m):
        try:
            return f"{float(vol(m)) / float(gt):.2f}"
        except Exception:
            return "—"
    header = ["method", "Dice", "IoU", "Hausdorff (mm)", "area err",
              "mesh vol (mm³)", "GT vol (mm³)", "ratio"]
    cols = [METHODS,
            [f(m, "dice") for m in METHODS], [f(m, "iou") for m in METHODS],
            [f(m, "hausdorff", 2) for m in METHODS], [f(m, "rel_area_err") for m in METHODS],
            [vol(m) or "—" for m in METHODS], [gt or "—"] * 3, [ratio(m) for m in METHODS]]
    return header, cols


def fig_selector(res, percase, recon):
    cases = sorted(os.path.basename(p)[5:-7] for p in glob.glob(str(res / "mesh_*_gt.off")))
    if not cases:
        return None, 0
    fig = make_subplots(rows=2, cols=4,
        specs=[[{"type": "scene"}] * 4, [{"type": "table", "colspan": 4}, None, None, None]],
        row_heights=[0.74, 0.26], vertical_spacing=0.06, horizontal_spacing=0.005,
        subplot_titles=[COL_TITLE[m] for m in COLS] + [""])
    trace_case = []
    for ci, case in enumerate(cases):
        for col, m in enumerate(COLS, start=1):
            off = res / f"mesh_{case}_{m}.off"
            if not off.exists():
                continue
            v, f = read_off(str(off))
            if len(f) == 0:
                continue
            fig.add_trace(go.Mesh3d(x=v[:, 0], y=v[:, 1], z=v[:, 2],
                          i=f[:, 0], j=f[:, 1], k=f[:, 2], color=COL_COLOR[m],
                          lighting=dict(ambient=0.5, diffuse=0.85), showscale=False,
                          hoverinfo="skip", visible=(ci == 0)), row=1, col=col)
            trace_case.append(ci)
        header, cols = case_table_cols(case, percase, recon)
        fig.add_trace(go.Table(header=dict(values=header, fill_color="#e8e8e8",
                      align="left", font=dict(size=12)),
                      cells=dict(values=cols, align="left", font=dict(size=12), height=24),
                      visible=(ci == 0)), row=2, col=1)
        trace_case.append(ci)
    ntr = len(trace_case)
    buttons = [dict(label=c, method="update",
                    args=[{"visible": [trace_case[t] == ci for t in range(ntr)]},
                          {"title.text": f"Tumor {c}  —  ground truth vs linear / spline / sdf"}])
               for ci, c in enumerate(cases)]
    fig.update_scenes(aspectmode="data", xaxis_visible=False, yaxis_visible=False,
                      zaxis_visible=False)
    fig.update_layout(height=820, margin=dict(l=0, r=0, t=110, b=0),
                      title=dict(text=f"Tumor {cases[0]}  —  ground truth vs linear / spline / sdf",
                                 x=0.5, xanchor="center", y=0.99),
                      updatemenus=[dict(buttons=buttons, direction="down", showactive=True,
                                        x=0.5, xanchor="center", y=1.10, yanchor="top",
                                        bgcolor="#f4f4f4")],
                      annotations=list(fig.layout.annotations) + [dict(
                          text="Select tumor:", x=0.5, xanchor="right", y=1.135,
                          yref="paper", xref="paper", showarrow=False,
                          font=dict(size=13), xshift=-120)])
    return fig, len(cases)


# ---------- HTML assembly ----------------------------------------------------
MATH = r"""
<h2>Methods &amp; mathematics</h2>
<p>A contour at axial slice <i>z</i> is a closed polyline
$C=(p_0,\dots,p_{n-1})$, $p_i=(x_i,y_i)$. Three quantities are reused: the
<b>signed (shoelace) area</b>
$$A(C)=\tfrac12\sum_{i=0}^{n-1}\bigl(x_i\,y_{i+1}-x_{i+1}\,y_i\bigr),$$
its sign giving orientation (we normalise every contour to counter-clockwise),
the cumulative <b>arc length</b> $s_i=\sum_{k<i}\lVert p_{k+1}-p_k\rVert$, and the
turning-angle <b>curvature</b> $\kappa_i=\dfrac{2\theta_i}{\lVert p_{i-1}-p_i\rVert+\lVert p_i-p_{i+1}\rVert}$.</p>

<h3>M1 — Linear (vertex-correspondence LERP)</h3>
<p>Resample $A,B$ to a common $N$ points equally spaced in arc length, find the
cyclic shift that best matches them,
$$\sigma^\star=\arg\min_{\sigma}\sum_{i=0}^{N-1}\bigl\lVert A_i-B_{(i+\sigma)\bmod N}\bigr\rVert^2,$$
then blend vertex-by-vertex at parameter $t\in[0,1]$:
$$C^{(i)}_{t}=(1-t)\,A_i+t\,B_{(i+\sigma^\star)\bmod N}.$$
Residual self-intersections are removed with a Bentley–Ottmann sweep.</p>

<h3>M2 — Polynomial / natural cubic spline (along z)</h3>
<p>Chaining the correspondence across slices makes each vertex trace a trajectory
$\gamma_i(z)=(x_i(z),y_i(z))$. Each coordinate is fitted independently. The
<b>natural cubic spline</b> is piecewise cubic with $C^2$ continuity; its second
derivatives $M_k$ solve the tridiagonal system
$$h_{k-1}M_{k-1}+2(h_{k-1}+h_k)M_k+h_kM_{k+1}
 =6\!\left(\frac{y_{k+1}-y_k}{h_k}-\frac{y_k-y_{k-1}}{h_{k-1}}\right),\quad M_0=M_{m}=0,$$
with $h_k=z_{k+1}-z_k$, and is evaluated as a cubic Hermite on each interval. The
<b>polynomial</b> variant is a least-squares fit of degree $d=\min(3,m-1)$, kept as a
teaching contrast (it under-fits / oscillates).</p>

<h3>M3 — Shape-based (signed distance field)</h3>
<p>Each contour becomes a signed distance field, $\Phi_C(\mathbf q)=\pm\min_{\mathbf c\in\partial C}\lVert\mathbf q-\mathbf c\rVert$
(positive inside; distance from a CGAL AABB-tree of the edges, sign from
$\texttt{Polygon\_2::bounded\_side}$). The map $\mathbf q\mapsto\min_{\mathbf c}\lVert\mathbf q-\mathbf c\rVert$
is the distance to the contour, whose nearest-edge partition is the <b>Voronoi diagram</b>
of the edges and whose ridge set is the <b>medial axis</b>. The fields are blended and the
new contour is the <b>zero level set</b>, extracted with marching squares:
$$\Phi_t=(1-t)\,\Phi_A+t\,\Phi_B,\qquad C_t=\{\mathbf q:\Phi_t(\mathbf q)=0\}.$$</p>

<h3>Reconstruction &amp; metrics</h3>
<p>Each densified contour stack becomes an oriented 3-D point cloud
$(x,y,k\,\Delta z)$ and is closed by <b>Poisson</b> surface reconstruction (which builds a
3-D Delaunay triangulation). Volume uses the divergence theorem,
$V=\tfrac16\sum_{\triangle}\mathbf a\cdot(\mathbf b\times\mathbf c)$. Interpolation accuracy
(leave-one-slice-out) uses <b>Dice</b> $=\dfrac{2\lvert X\cap Y\rvert}{\lvert X\rvert+\lvert Y\rvert}$,
<b>IoU</b> $=\dfrac{\lvert X\cap Y\rvert}{\lvert X\cup Y\rvert}$, the symmetric <b>Hausdorff</b>
distance, and the relative area error $\lvert A(\hat C)-A(C)\rvert/A(C)$.</p>
"""


def tbl(rows, headers, title, note=""):
    th = "".join(f"<th>{h}</th>" for h in headers)
    trs = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    n = f"<p class='note'>{note}</p>" if note else ""
    return f"<h3>{title}</h3>{n}<table><thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table>"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    res = Path(args.results)

    agg = load_agg(res)
    raw = load_raw(res)
    percase = load_percase(res)
    recon = load_recon(res)
    ncases = len(percase)
    nslices = len([r for r in raw if r["method"] == "linear"])

    bar = fig_bars(agg).to_html(full_html=False, include_plotlyjs=False)
    cdf = fig_cdf(raw).to_html(full_html=False, include_plotlyjs=False)
    line = fig_percase(percase).to_html(full_html=False, include_plotlyjs=False)
    sel_fig, n_vis = fig_selector(res, percase, recon)
    sel = sel_fig.to_html(full_html=False, include_plotlyjs=False) if sel_fig else ""

    # tables
    agg_rows = [[a["method"], a["n"], f"{float(a['dice_mean']):.3f}",
                 f"{float(a['dice_median']):.3f}", f"{float(a['dice_pct_gt80']):.1f}%",
                 f"{float(a['iou_mean']):.3f}", f"{float(a['hausdorff_mean']):.2f}",
                 f"{float(a['area_err_median']):.3f}"] for a in agg]
    agg_tbl = tbl(agg_rows, ["method", "n", "Dice mean", "Dice median", "% Dice>0.8",
                             "IoU mean", "Hausdorff mean (mm)", "area-err median"],
                  "Leave-one-slice-out accuracy (aggregate)")
    rec_rows = []
    for (case, m), r in sorted(recon.items()):
        mv = r["mesh_volume_mm3"].replace(" mm^3", "").strip()
        gt = r["voxel_gt_volume_mm3"].strip()
        try:
            ratio = f"{float(mv)/float(gt):.2f}"
        except Exception:
            ratio = "—"
        rec_rows.append([case, m, r.get("status", ""), mv or "—", gt, ratio])
    rec_tbl = tbl(rec_rows, ["case", "method", "status", "mesh vol (mm³)",
                             "voxel GT (mm³)", "ratio"],
                  "Reconstruction volume vs voxel ground truth (3D cases)")

    style = """<style>
      body{font-family:system-ui,Arial,sans-serif;margin:24px;color:#222;max-width:1180px}
      h1{font-size:24px} h2{margin-top:34px;border-bottom:1px solid #ddd;padding-bottom:4px}
      h3{margin-top:22px} p{line-height:1.5;max-width:900px}
      table{border-collapse:collapse;margin:8px 0 18px}
      th,td{border:1px solid #ccc;padding:4px 10px;text-align:right;font-size:13px}
      th{background:#f0f0f0} td:first-child,th:first-child{text-align:left}
      .note{color:#555;font-size:13px}
      .fig{margin:6px 0 26px}
    </style>"""
    mathjax = ("<script>MathJax={tex:{inlineMath:[['$','$'],['\\\\(','\\\\)']],"
               "displayMath:[['$$','$$'],['\\\\[','\\\\]']]}};</script>"
               "<script async src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>")
    head = (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>BraTS interpolation report</title>{style}{mathjax}"
            f"<script>{get_plotlyjs()}</script></head><body>")
    intro = (f"<h1>Tumor contour interpolation, 3-D reconstruction &amp; comparison</h1>"
             f"<p>Real-data evaluation on <b>BraTS 2024 GLI</b> (enhancing tumor). "
             f"Leave-one-slice-out over <b>{ncases} cases</b> "
             f"(<b>{nslices} held-out slices</b> per method); "
             f"<b>{n_vis}</b> tumors reconstructed in 3-D below.</p>")
    body = (intro + MATH +
            "<h2>Metrics — bar plots</h2><div class='fig'>" + bar + "</div>" +
            "<p class='note'>Left: region overlap (Dice/IoU, higher is better). "
            "Centre: mean symmetric Hausdorff distance (mm, lower is better). "
            "Right: fraction of held-out slices reproduced with Dice &gt; 0.8.</p>" +
            "<h2>Metrics — line plots</h2>"
            "<h3>Dice distribution (empirical CDF)</h3><div class='fig'>" + cdf + "</div>"
            "<h3>Per-case mean Dice</h3><div class='fig'>" + line + "</div>" +
            "<h2>3-D reconstructions (select a tumor)</h2><div class='fig'>" + sel + "</div>" +
            "<h2>Tables</h2>" + agg_tbl + rec_tbl +
            "</body></html>")
    Path(args.out).write_text(head + body, encoding="utf-8")
    print(f"wrote {args.out}  ({os.path.getsize(args.out)/1e6:.1f} MB, "
          f"{ncases} cases, {n_vis} in 3D)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
