#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


def parse_obj(path: Path):
    vertices = []
    edges = set()
    if not path.exists() or path.stat().st_size == 0:
        return None, None

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            kind = parts[0]

            if kind == "v" and len(parts) >= 3:
                try:
                    x = float(parts[1])
                    y = float(parts[2])
                    vertices.append((x, y))
                except ValueError:
                    continue
            elif kind == "f" and len(parts) >= 4:
                idx = []
                for token in parts[1:]:
                    item = token.split("/")[0]
                    if not item:
                        continue
                    try:
                        idx.append(int(item))
                    except ValueError:
                        continue
                if len(idx) < 2:
                    continue
                for i in range(len(idx)):
                    a = idx[i]
                    b = idx[(i + 1) % len(idx)]
                    if a == b:
                        continue
                    edges.add(tuple(sorted((a, b))))
            elif kind == "l" and len(parts) >= 3:
                idx = []
                for token in parts[1:]:
                    item = token.split("/")[0]
                    if not item:
                        continue
                    try:
                        idx.append(int(item))
                    except ValueError:
                        continue
                for i in range(len(idx) - 1):
                    a, b = idx[i], idx[i + 1]
                    if a == b:
                        continue
                    edges.add(tuple(sorted((a, b))))

    if not vertices:
        return None, None

    return np.asarray(vertices, dtype=float), edges


def mesh_style_for_density(n_verts: int):
    if n_verts < 1_000:
        return 0.6, 3.0
    if n_verts < 10_000:
        return 0.3, 1.0
    if n_verts < 50_000:
        return 0.15, 0.4
    return 0.08, 0.2


def draw_mesh(ax, vertices, edges, title, force_square_box=False):
    ax.set_axis_off()
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)

    if vertices is None:
        ax.text(0.5, 0.5, "archivo no disponible", ha="center", va="center", transform=ax.transAxes)
        return

    lw, ms = mesh_style_for_density(len(vertices))
    segments = []
    for a, b in edges:
        ia = a - 1
        ib = b - 1
        if ia < 0 or ib < 0 or ia >= len(vertices) or ib >= len(vertices):
            continue
        pa = vertices[ia]
        pb = vertices[ib]
        segments.append(((pa[0], pa[1]), (pb[0], pb[1])))

    if segments:
        lc = LineCollection(segments, colors="#666666", linewidths=lw)
        ax.add_collection(lc)

    x_min = float(np.min(vertices[:, 0]))
    x_max = float(np.max(vertices[:, 0]))
    y_min = float(np.min(vertices[:, 1]))
    y_max = float(np.max(vertices[:, 1]))
    x_span = x_max - x_min
    y_span = y_max - y_min
    x_pad = 0.02 * (x_span if x_span > 0 else 1.0)
    y_pad = 0.02 * (y_span if y_span > 0 else 1.0)
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    if force_square_box:
        ax.set_box_aspect(1)

    ax.scatter(vertices[:, 0], vertices[:, 1], s=ms, c="black")


def draw_heightmap(ax, input_png: Path, title: str, fill_subplot=False):
    ax.set_axis_off()
    ax.set_title(title, fontsize=10)
    try:
        img = mpimg.imread(str(input_png))
        if img.ndim == 3:
            img = img[..., :3].mean(axis=2)
        if fill_subplot:
            ax.imshow(img, cmap="gray", aspect="auto")
            ax.set_box_aspect(1)
        else:
            ax.imshow(img, cmap="gray")
    except Exception:
        ax.text(0.5, 0.5, "archivo no disponible", ha="center", va="center", transform=ax.transAxes)


def save_single_heightmap(input_png: Path, output_file: Path, name: str):
    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    draw_heightmap(ax, input_png, f"{name} - Heightmap")
    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)


def save_single_mesh(obj_path: Path, output_file: Path, title: str):
    vertices, edges = parse_obj(obj_path)
    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    draw_mesh(ax, vertices, edges, title)
    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)
    return 0 if vertices is None else len(vertices)


def save_comparative(
    input_png: Path,
    original_obj: Path,
    simplified_obj: Path,
    output_file: Path,
    name: str,
    order=None,
    gamma=None,
):
    v_original, e_original = parse_obj(original_obj)
    v_simplified, e_simplified = parse_obj(simplified_obj)

    n = 0 if v_original is None else len(v_original)
    m = 0 if v_simplified is None else len(v_simplified)
    reduction = 0.0
    if n > 0:
        reduction = 100.0 * (1.0 - (m / n))

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18, 6),
        dpi=150,
        gridspec_kw={"width_ratios": [1, 1, 1]},
    )
    if order is not None and gamma is not None:
        suptitle = f"{name}  -  order={order}, gamma={gamma}"
    else:
        suptitle = name
    fig.suptitle(suptitle, fontsize=14, y=0.98)

    draw_heightmap(axes[0], input_png, f"{name} - Heightmap", fill_subplot=True)
    draw_mesh(axes[1], v_original, e_original, f"{name} - Original ({n} vertices)", force_square_box=True)
    draw_mesh(
        axes[2],
        v_simplified,
        e_simplified,
        f"{name} - Simplificado ({m} vertices, {reduction:.2f}% reduccion)",
        force_square_box=True,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_file)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-png", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--order", type=int, default=None)
    parser.add_argument("--gamma", type=float, default=None)
    args = parser.parse_args()

    input_png = Path(args.input_png)
    output_dir = Path(args.output_dir)
    name = args.name

    output_dir.mkdir(parents=True, exist_ok=True)

    original_obj = output_dir / "original.obj"
    simplified_obj = output_dir / "simplificado.obj"

    save_single_heightmap(input_png, output_dir / "00_heightmap.png", name)

    v_original, e_original = parse_obj(original_obj)
    v_simplified, e_simplified = parse_obj(simplified_obj)

    n_original = 0 if v_original is None else len(v_original)
    n_simplified = 0 if v_simplified is None else len(v_simplified)
    reduction = 100.0 * (1.0 - (n_simplified / n_original)) if n_original > 0 else 0.0

    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    draw_mesh(ax, v_original, e_original if e_original is not None else set(), f"{name} - Original ({n_original} vertices)")
    fig.tight_layout()
    fig.savefig(output_dir / "01_original.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    draw_mesh(
        ax,
        v_simplified,
        e_simplified if e_simplified is not None else set(),
        f"{name} - Simplificado ({n_simplified} vertices, {reduction:.2f}% reduccion)",
    )
    fig.tight_layout()
    fig.savefig(output_dir / "02_simplificado.png")
    plt.close(fig)

    save_comparative(
        input_png,
        original_obj,
        simplified_obj,
        output_dir / "03_comparativo.png",
        name,
        order=args.order,
        gamma=args.gamma,
    )


if __name__ == "__main__":
    main()
