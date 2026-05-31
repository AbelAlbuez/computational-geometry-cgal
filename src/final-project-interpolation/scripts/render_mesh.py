#!/usr/bin/env python3
# =============================================================================
# Render a reconstructed surface mesh (.off) to PNG snapshots — a headless 3D
# preview of the Poisson output (no OpenGL needed; matplotlib software render).
# Two viewing angles per figure, simple Lambert shading from face normals.
#
#   python3 render_mesh.py mesh.off out.png [--title "linear"]
# =============================================================================
from __future__ import annotations

import argparse
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def read_off(path: str):
    """Minimal ASCII OFF reader -> (vertices Nx3, triangles Mx3)."""
    with open(path) as f:
        toks = f.read().split()
    if not toks or not toks[0].upper().startswith("OFF"):
        raise ValueError(f"{path}: not an OFF file")
    # The header keyword may be glued to the counts ("OFF") or standalone.
    p = 1 if toks[0].upper() == "OFF" else 0
    if p == 0:                       # e.g. "OFF123 ..." is not expected; be safe
        toks = toks[0][3:] and ([toks[0][3:]] + toks[1:]) or toks[1:]
        p = 0
    nv, nf = int(toks[p]), int(toks[p + 1])
    vals = toks[p + 3:]
    verts = np.array(vals[:nv * 3], dtype=float).reshape(nv, 3)
    rest = vals[nv * 3:]
    faces = []
    q = 0
    for _ in range(nf):
        k = int(rest[q]); ids = list(map(int, rest[q + 1:q + 1 + k])); q += 1 + k
        for t in range(1, k - 1):    # fan-triangulate polygons (meshes are tris)
            faces.append((ids[0], ids[t], ids[t + 1]))
    return verts, np.asarray(faces, dtype=int)


def render(verts, faces, out: str, title: str):
    tris = verts[faces]
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    ln = np.linalg.norm(n, axis=1); ln[ln == 0] = 1.0
    n /= ln[:, None]
    light = np.array([0.3, 0.3, 0.9]); light /= np.linalg.norm(light)
    shade = 0.45 + 0.55 * np.clip(n @ light, 0, 1)
    colors = np.zeros((len(faces), 4))
    colors[:, 0] = 0.40 * shade; colors[:, 1] = 0.55 * shade
    colors[:, 2] = 0.85 * shade; colors[:, 3] = 1.0

    mn, mx = verts.min(0), verts.max(0)
    c = (mn + mx) / 2.0; r = (mx - mn).max() / 2.0 + 1e-6

    fig = plt.figure(figsize=(11, 5))
    for i, (elev, azim) in enumerate([(20, 35), (18, 125)]):
        ax = fig.add_subplot(1, 2, i + 1, projection="3d")
        coll = Poly3DCollection(tris, linewidths=0.05)
        coll.set_facecolor(colors)
        coll.set_edgecolor((0, 0, 0, 0.08))
        ax.add_collection3d(coll)
        ax.set_xlim(c[0] - r, c[0] + r)
        ax.set_ylim(c[1] - r, c[1] + r)
        ax.set_zlim(c[2] - r, c[2] + r)
        try:
            ax.set_box_aspect((1, 1, 1))
        except Exception:
            pass
        ax.view_init(elev, azim)
        ax.set_axis_off()
    fig.suptitle(f"{title}   (V={len(verts)}  F={len(faces)})")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("out")
    ap.add_argument("--title", default="")
    args = ap.parse_args()
    verts, faces = read_off(args.mesh)
    if len(faces) == 0:
        print(f"{args.mesh}: no faces", file=sys.stderr); return 1
    render(verts, faces, args.out, args.title or args.mesh)
    print(f"wrote {args.out}  ({len(verts)} verts, {len(faces)} tris)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
