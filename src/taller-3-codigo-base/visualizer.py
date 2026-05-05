#!/usr/bin/env python3
"""
visualizer.py - Taller 3: Geometria Computacional
Genera:
  - output/resultado.gif    : antes vs despues (2 frames)
  - output/paso_a_paso.gif  : 7 pasos del algoritmo
"""

import argparse
import os

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from PIL import Image

PASOS = [
    (
        "Paso 1 - Malla original",
        "Un vertice por pixel. Triangulacion de Delaunay completa.\nMiles de triangulos, muchos redundantes.",
    ),
    (
        "Paso 2 - Calcular error de planitud",
        "error(v) = promedio |v.z - qi.z| para cada vecino qi.\nVertices en zonas planas tienen error bajo (verde).",
    ),
    (
        "Paso 3 - Min-heap: sacar vertice mas plano",
        "El priority_queue entrega el vertice con menor error.\nDijkstra: siempre procesar el mas prescindible primero.",
    ),
    (
        "Paso 4 - Recorrer estrella con Face_circulator",
        "DCEL half-edge: girar alrededor de p con twin->next.\nRecolectar vecinos directos (anillo orden 1).",
    ),
    (
        "Paso 5 - Expandir vecindad a orden k",
        "Repetir circulator desde cada vecino qi.\nAcumular triangulos sin duplicados (std::set).",
    ),
    (
        "Paso 6 - Eliminar vertice redundante",
        "T.remove(p) -> CGAL borra estrella y re-triangula.\nSolo si error(p) < e. Bordes y crestas se conservan.",
    ),
    (
        "Paso 7 - Resultado final",
        "Loop hasta heap vacio o error_min >= e.\nZonas planas simplificadas. Relieve preservado.",
    ),
]


def load_obj(path):
    verts, faces = [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] == "v":
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif parts[0] == "f":
                face = [int(p.split("/")[0]) - 1 for p in parts[1:]]
                if len(face) == 3:
                    faces.append(face)
    return np.array(verts), np.array(faces)


def render_mesh(verts, faces, title, subtitle, output_path, cmap="terrain"):
    fig, ax = plt.subplots(figsize=(6, 6.5), facecolor="#1a1a1a")
    ax.set_facecolor("#1a1a1a")
    if len(faces) > 0:
        x, y, z = verts[:, 0], verts[:, 1], verts[:, 2]
        triang = mtri.Triangulation(x, y, faces)
        ax.tripcolor(triang, z, cmap=cmap, shading="gouraud", alpha=0.92)
        ax.triplot(triang, color="#378ADD", linewidth=0.25, alpha=0.5)
    ax.set_title(title, color="white", fontsize=13, fontweight="bold", pad=8)
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


def render_paso(verts, faces, paso_num, titulo, descripcion, output_path, cmap="terrain"):
    fig = plt.figure(figsize=(6, 7.5), facecolor="#1a1a1a")
    ax_mesh = fig.add_axes([0.05, 0.28, 0.9, 0.65])
    ax_mesh.set_facecolor("#1a1a1a")
    if len(faces) > 0:
        x, y, z = verts[:, 0], verts[:, 1], verts[:, 2]
        triang = mtri.Triangulation(x, y, faces)
        ax_mesh.tripcolor(triang, z, cmap=cmap, shading="gouraud", alpha=0.92)
        ax_mesh.triplot(triang, color="#378ADD", linewidth=0.25, alpha=0.5)
    ax_mesh.set_title(titulo, color="white", fontsize=12, fontweight="bold", pad=6)
    ax_mesh.axis("equal")
    ax_mesh.axis("off")

    ax_text = fig.add_axes([0.05, 0.02, 0.9, 0.24])
    ax_text.set_facecolor("#111111")
    ax_text.axis("off")
    ax_text.text(
        0.5,
        0.5,
        descripcion,
        color="#D3D1C7",
        fontsize=9.5,
        ha="center",
        va="center",
        transform=ax_text.transAxes,
        wrap=True,
        multialignment="center",
        fontfamily="monospace",
    )

    fig.text(0.5, 0.96, f"Paso {paso_num} de 7", color="#378ADD", fontsize=10, ha="center", va="top")
    plt.savefig(output_path, dpi=110, bbox_inches="tight", facecolor="#1a1a1a")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="output")
    args = parser.parse_args()
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    orig_path = os.path.join(output_dir, "original.obj")
    simp_path = os.path.join(output_dir, "simplificado.obj")
    if not os.path.exists(orig_path) or not os.path.exists(simp_path):
        print("Error: no se encontraron los .obj en", output_dir)
        return

    print("Leyendo mallas...")
    vo, fo = load_obj(orig_path)
    vs, fs = load_obj(simp_path)

    print("Renderizando resultado...")
    p1 = os.path.join(output_dir, "frame_original.png")
    p2 = os.path.join(output_dir, "frame_simplificado.png")
    render_mesh(vo, fo, f"Original ({len(vo)} vertices)", "", p1, cmap="terrain")
    render_mesh(vs, fs, f"Simplificado ({len(vs)} vertices)", "", p2, cmap="terrain")
    frames = [Image.open(p1), Image.open(p2)]
    gif1 = os.path.join(output_dir, "resultado.gif")
    frames[0].save(gif1, save_all=True, append_images=frames[1:], duration=1500, loop=0)
    print(f"GIF guardado: {gif1}")

    print("Renderizando paso a paso...")
    mallas = [(vo, fo), (vs, fs), (vs, fs), (vs, fs), (vs, fs), (vs, fs), (vs, fs)]
    paso_paths = []
    for i, (titulo, desc) in enumerate(PASOS):
        p = os.path.join(output_dir, f"paso_{i + 1:02d}.png")
        v, f = mallas[i]
        render_paso(v, f, i + 1, titulo, desc, p, cmap="terrain")
        paso_paths.append(p)
        print(f"  Frame {i + 1}/7 OK")

    frames2 = [Image.open(p) for p in paso_paths]
    gif2 = os.path.join(output_dir, "paso_a_paso.gif")
    frames2[0].save(gif2, save_all=True, append_images=frames2[1:], duration=2000, loop=0)
    print(f"GIF guardado: {gif2}")


if __name__ == "__main__":
    main()
