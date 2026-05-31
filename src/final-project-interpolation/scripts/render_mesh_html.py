#!/usr/bin/env python3
# =============================================================================
# Export a reconstructed surface (.off) to a self-contained interactive 3D HTML
# (Plotly Mesh3d). Open the .html in any browser and rotate / zoom / pan — no
# server, no app, works offline (plotly.js is embedded).
#
#   python3 render_mesh_html.py mesh.off out.html [--title "..."]
# =============================================================================
from __future__ import annotations

import argparse
import sys

import numpy as np


def read_off(path: str):
    with open(path) as f:
        toks = f.read().split()
    if not toks or not toks[0].upper().startswith("OFF"):
        raise ValueError(f"{path}: not an OFF file")
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("out")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    import plotly.graph_objects as go
    v, f = read_off(args.mesh)
    if len(f) == 0:
        print(f"{args.mesh}: no faces", file=sys.stderr); return 1
    x, y, z = v.T
    i, j, k = f.T
    fig = go.Figure(data=[go.Mesh3d(
        x=x, y=y, z=z, i=i, j=j, k=k,
        color="lightsteelblue", opacity=1.0, flatshading=False,
        lighting=dict(ambient=0.45, diffuse=0.85, specular=0.15),
        lightposition=dict(x=200, y=200, z=400))])
    fig.update_layout(
        title=f"{args.title}  (V={len(v)} F={len(f)})",
        scene=dict(aspectmode="data", xaxis_title="x (mm)",
                   yaxis_title="y (mm)", zaxis_title="z (mm)"),
        margin=dict(l=0, r=0, t=40, b=0))
    fig.write_html(args.out, include_plotlyjs=True, full_html=True)
    print(f"wrote {args.out}  ({len(v)} verts, {len(f)} tris)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
