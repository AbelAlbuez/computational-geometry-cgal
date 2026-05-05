#!/usr/bin/env python3
"""
visualizer.py — Taller 3: Geometria Computacional
Lee output/original.obj y output/simplificado.obj,
renderiza ambas mallas como imagenes y genera output/resultado.gif
"""

import os

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from PIL import Image


def load_obj(path):
    """Carga vertices y caras de un .obj simple."""
    verts, faces = [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] == "v":
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif parts[0] == "f":
                # indices 1-based en .obj -> 0-based
                face = [int(p.split("/")[0]) - 1 for p in parts[1:]]
                if len(face) == 3:
                    faces.append(face)
    return np.array(verts), np.array(faces)


def render_mesh(verts, faces, title, output_path):
    """Renderiza la malla como PNG vista desde arriba (x,y coloreado por z)."""
    fig, ax = plt.subplots(figsize=(6, 6), facecolor="#1a1a1a")
    ax.set_facecolor("#1a1a1a")

    if len(faces) > 0:
        x, y, z = verts[:, 0], verts[:, 1], verts[:, 2]
        triang = mtri.Triangulation(x, y, faces)
        ax.tripcolor(triang, z, cmap="gray", shading="gouraud", alpha=0.9)
        ax.triplot(triang, color="#378ADD", linewidth=0.3, alpha=0.6)

    ax.set_title(title, color="white", fontsize=13, pad=10)
    ax.set_xlabel(
        f"{len(verts)} vertices  |  {len(faces)} triangulos",
        color="#888780",
        fontsize=10,
    )
    ax.axis("equal")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(output_path, dpi=120, bbox_inches="tight", facecolor="#1a1a1a")
    plt.close()


def main():
    os.makedirs("output", exist_ok=True)

    original_path = "output/original.obj"
    simplificado_path = "output/simplificado.obj"

    if not os.path.exists(original_path) or not os.path.exists(simplificado_path):
        print("Error: no se encontraron los archivos .obj en output/")
        return

    print("Leyendo original...")
    verts_orig, faces_orig = load_obj(original_path)
    print("Leyendo simplificado...")
    verts_simp, faces_simp = load_obj(simplificado_path)

    print("Renderizando frames...")
    render_mesh(
        verts_orig,
        faces_orig,
        f"Original ({len(verts_orig)} vertices)",
        "output/frame_original.png",
    )
    render_mesh(
        verts_simp,
        faces_simp,
        f"Simplificado ({len(verts_simp)} vertices)",
        "output/frame_simplificado.png",
    )

    print("Generando GIF...")
    frames = [
        Image.open("output/frame_original.png"),
        Image.open("output/frame_simplificado.png"),
    ]
    frames[0].save(
        "output/resultado.gif",
        save_all=True,
        append_images=frames[1:],
        duration=1200,
        loop=0,
    )
    print("GIF guardado en output/resultado.gif")


if __name__ == "__main__":
    main()
