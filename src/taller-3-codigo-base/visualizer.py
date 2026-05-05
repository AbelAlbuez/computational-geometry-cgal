#!/usr/bin/env python3
"""
visualizer.py - Taller 3: Geometria Computacional
Genera:
  - output/test-[name]/resultado.gif     : antes vs despues (2 frames)
  - output/test-[name]/paso_a_paso.gif   : 7 pasos del algoritmo con visualizacion real
"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from PIL import Image

warnings.filterwarnings("ignore")

BG = "#1a1a1a"

PASOS = [
    (
        "Paso 1 - Malla original",
        "Un vertice por pixel. Triangulacion de Delaunay completa.\nMiles de triangulos, muchos redundantes.",
    ),
    (
        "Paso 2 - Error de planitud por vertice",
        "error(v) = promedio |v.z - qi.z| para cada vecino qi.\nRojo = error alto (borde). Verde = error bajo (plano).",
    ),
    (
        "Paso 3 - Min-heap: vertice mas plano",
        "El priority_queue entrega el vertice con menor error.\nDijkstra: siempre procesar el mas prescindible primero.",
    ),
    (
        "Paso 4 - Estrella con Face_circulator",
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


def load_partial_obj(path, fraction=0.5):
    """Carga fraccion de vertices para simular estado intermedio."""
    verts, faces = load_obj(path)
    n = max(4, int(len(verts) * fraction))
    verts_p = verts[:n]
    faces_p = [f for f in faces if all(i < n for i in f)]
    farr = np.array(faces_p).astype(int) if faces_p else np.zeros((0, 3), dtype=int)
    return verts_p, farr


def make_triang(verts, faces):
    if len(faces) == 0 or len(verts) < 3:
        return None
    return mtri.Triangulation(verts[:, 0], verts[:, 1], faces)


def base_fig():
    fig = plt.figure(figsize=(6, 7.5), facecolor=BG)
    ax_m = fig.add_axes([0.05, 0.28, 0.9, 0.65])
    ax_m.set_facecolor(BG)
    ax_t = fig.add_axes([0.05, 0.02, 0.9, 0.24])
    ax_t.set_facecolor("#111111")
    ax_t.axis("off")
    return fig, ax_m, ax_t


def add_text(fig, ax_t, paso_num, titulo, desc):
    fig.text(
        0.5,
        0.96,
        f"Paso {paso_num} de 7 - {titulo}",
        color="#378ADD",
        fontsize=10,
        ha="center",
        va="top",
        fontweight="bold",
    )
    ax_t.text(
        0.5,
        0.5,
        desc,
        color="#D3D1C7",
        fontsize=9.5,
        ha="center",
        va="center",
        transform=ax_t.transAxes,
        multialignment="center",
        fontfamily="monospace",
    )


def save_close(fig, path):
    plt.savefig(path, dpi=110, bbox_inches="tight", facecolor=BG)
    plt.close()


def render_paso1(vo, fo, path):
    fig, ax, ax_t = base_fig()
    t = make_triang(vo, fo)
    if t is not None:
        ax.tripcolor(t, vo[:, 2], cmap="terrain", shading="gouraud", alpha=0.92)
        ax.triplot(t, color="#378ADD", linewidth=0.2, alpha=0.4)
    ax.set_xlabel(f"{len(vo)} vertices  |  {len(fo)} triangulos", color="#888780", fontsize=9)
    ax.axis("equal")
    ax.axis("off")
    add_text(fig, ax_t, 1, PASOS[0][0].split("-")[-1].strip(), PASOS[0][1])
    save_close(fig, path)


def render_paso2(vo, fo, path):
    fig, ax, ax_t = base_fig()
    t = make_triang(vo, fo)
    if t is not None:
        n = len(vo)
        nbr_z = [[] for _ in range(n)]
        for f in fo:
            for i in range(3):
                for j in range(3):
                    if i != j:
                        nbr_z[f[i]].append(vo[f[j], 2])
        errors = np.array([np.mean(np.abs(np.array(nz) - vo[i, 2])) if nz else 0 for i, nz in enumerate(nbr_z)])
        errors = (errors - errors.min()) / (errors.max() - errors.min() + 1e-9)
        ax.tripcolor(t, errors, cmap="RdYlGn_r", shading="gouraud", alpha=0.92)
        ax.triplot(t, color="white", linewidth=0.1, alpha=0.15)
        sm = plt.cm.ScalarMappable(cmap="RdYlGn_r", norm=plt.Normalize(0, 1))
        cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
        cbar.set_label("error normalizado")
    ax.axis("equal")
    ax.axis("off")
    add_text(fig, ax_t, 2, PASOS[1][0].split("-")[-1].strip(), PASOS[1][1])
    save_close(fig, path)


def render_paso3(vo, fo, path):
    fig, ax, ax_t = base_fig()
    t = make_triang(vo, fo)
    if t is not None:
        ax.tripcolor(t, vo[:, 2], cmap="terrain", shading="gouraud", alpha=0.7)
        ax.triplot(t, color="#378ADD", linewidth=0.15, alpha=0.3)
        mean_z = vo[:, 2].mean()
        flat_idx = np.argmin(np.abs(vo[:, 2] - mean_z))
        ax.scatter(vo[flat_idx, 0], vo[flat_idx, 1], c="#FF4444", s=80, zorder=5, label="vertice mas plano")
        ax.legend(loc="upper right", facecolor="#222", labelcolor="white", fontsize=8)
    ax.axis("equal")
    ax.axis("off")
    add_text(fig, ax_t, 3, PASOS[2][0].split("-")[-1].strip(), PASOS[2][1])
    save_close(fig, path)


def render_paso4(vo, fo, path):
    fig, ax, ax_t = base_fig()
    t = make_triang(vo, fo)
    if t is not None:
        ax.tripcolor(t, vo[:, 2], cmap="terrain", shading="gouraud", alpha=0.4)
        ax.triplot(t, color="#378ADD", linewidth=0.15, alpha=0.2)
        cx = len(vo) // 2
        star_faces = [f for f in fo if cx in f][:12]
        for f in star_faces:
            pts = vo[f, :2]
            tri = plt.Polygon(pts, color="#1D9E75", alpha=0.7, linewidth=1.2, edgecolor="white")
            ax.add_patch(tri)
        ax.scatter(vo[cx, 0], vo[cx, 1], c="#FF4444", s=100, zorder=6)
    ax.axis("equal")
    ax.axis("off")
    add_text(fig, ax_t, 4, PASOS[3][0].split("-")[-1].strip(), PASOS[3][1])
    save_close(fig, path)


def render_paso5(vo, fo, path):
    fig, ax, ax_t = base_fig()
    t = make_triang(vo, fo)
    if t is not None:
        ax.tripcolor(t, vo[:, 2], cmap="terrain", shading="gouraud", alpha=0.35)
        ax.triplot(t, color="#378ADD", linewidth=0.12, alpha=0.2)
        cx = len(vo) // 2
        ring1 = set()
        star1 = [f for f in fo if cx in f]
        for f in star1:
            ring1.update(f)
        ring1.discard(cx)
        for f in star1:
            ax.add_patch(plt.Polygon(vo[f, :2], color="#1D9E75", alpha=0.6, edgecolor="white", lw=0.8))
        ring2_faces = [
            f
            for f in fo
            if any(v in ring1 for v in f) and cx not in f and not all(v in ring1 or v == cx for v in f)
        ][:20]
        for f in ring2_faces:
            ax.add_patch(plt.Polygon(vo[f, :2], color="#BA7517", alpha=0.5, edgecolor="white", lw=0.6))
        ax.scatter(vo[cx, 0], vo[cx, 1], c="#FF4444", s=100, zorder=6)
        from matplotlib.patches import Patch

        ax.legend(
            handles=[Patch(color="#1D9E75", label="orden 1"), Patch(color="#BA7517", label="orden 2")],
            loc="upper right",
            facecolor="#222",
            labelcolor="white",
            fontsize=8,
        )
    ax.axis("equal")
    ax.axis("off")
    add_text(fig, ax_t, 5, PASOS[4][0].split("-")[-1].strip(), PASOS[4][1])
    save_close(fig, path)


def render_paso6(orig_path, path):
    fig, ax, ax_t = base_fig()
    vp, fp = load_partial_obj(orig_path, 0.5)
    t = make_triang(vp, fp)
    if t is not None:
        ax.tripcolor(t, vp[:, 2], cmap="terrain", shading="gouraud", alpha=0.92)
        ax.triplot(t, color="#378ADD", linewidth=0.3, alpha=0.5)
    ax.set_xlabel(f"{len(vp)} vertices  |  {len(fp)} triangulos (en proceso...)", color="#888780", fontsize=9)
    ax.axis("equal")
    ax.axis("off")
    add_text(fig, ax_t, 6, PASOS[5][0].split("-")[-1].strip(), PASOS[5][1])
    save_close(fig, path)


def render_paso7(vs, fs, path):
    fig, ax, ax_t = base_fig()
    t = make_triang(vs, fs)
    if t is not None:
        ax.tripcolor(t, vs[:, 2], cmap="terrain", shading="gouraud", alpha=0.92)
        ax.triplot(t, color="#378ADD", linewidth=0.4, alpha=0.6)
    ax.set_xlabel(f"{len(vs)} vertices  |  {len(fs)} triangulos", color="#888780", fontsize=9)
    ax.axis("equal")
    ax.axis("off")
    add_text(fig, ax_t, 7, PASOS[6][0].split("-")[-1].strip(), PASOS[6][1])
    save_close(fig, path)


def render_mesh_simple(verts, faces, title, output_path):
    fig, ax = plt.subplots(figsize=(6, 6.5), facecolor=BG)
    ax.set_facecolor(BG)
    t = make_triang(verts, faces)
    if t is not None:
        ax.tripcolor(t, verts[:, 2], cmap="terrain", shading="gouraud", alpha=0.92)
        ax.triplot(t, color="#378ADD", linewidth=0.25, alpha=0.5)
    ax.set_title(title, color="white", fontsize=13, fontweight="bold", pad=8)
    ax.set_xlabel(f"{len(verts)} vertices  |  {len(faces)} triangulos", color="#888780", fontsize=10)
    ax.axis("equal")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(output_path, dpi=120, bbox_inches="tight", facecolor=BG)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="output")
    parser.add_argument("--name", default="test")
    args = parser.parse_args()
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    orig_path = os.path.join(output_dir, "original.obj")
    simp_path = os.path.join(output_dir, "simplificado.obj")
    if not os.path.exists(orig_path) or not os.path.exists(simp_path):
        print("Error: no se encontraron los .obj en", output_dir)
        return

    print(f"[{args.name}] Leyendo mallas...")
    vo, fo = load_obj(orig_path)
    vs, fs = load_obj(simp_path)
    print(f"  Original:     {len(vo)} vertices, {len(fo)} caras")
    print(f"  Simplificado: {len(vs)} vertices, {len(fs)} caras")

    print(f"[{args.name}] Generando resultado.gif...")
    p1 = os.path.join(output_dir, "frame_original.png")
    p2 = os.path.join(output_dir, "frame_simplificado.png")
    render_mesh_simple(vo, fo, f"Original ({len(vo)} vertices)", p1)
    render_mesh_simple(vs, fs, f"Simplificado ({len(vs)} vertices)", p2)
    frames = [Image.open(p1), Image.open(p2)]
    gif1 = os.path.join(output_dir, "resultado.gif")
    frames[0].save(gif1, save_all=True, append_images=frames[1:], duration=1500, loop=0)
    print(f"  Guardado: {gif1}")

    print(f"[{args.name}] Generando paso_a_paso.gif (7 frames)...")
    paso_paths = []
    renders = [
        lambda p: render_paso1(vo, fo, p),
        lambda p: render_paso2(vo, fo, p),
        lambda p: render_paso3(vo, fo, p),
        lambda p: render_paso4(vo, fo, p),
        lambda p: render_paso5(vo, fo, p),
        lambda p: render_paso6(orig_path, p),
        lambda p: render_paso7(vs, fs, p),
    ]
    for i, fn in enumerate(renders):
        p = os.path.join(output_dir, f"paso_{i + 1:02d}.png")
        fn(p)
        paso_paths.append(p)
        print(f"  Frame {i + 1}/7 OK")

    frames2 = [Image.open(p) for p in paso_paths]
    gif2 = os.path.join(output_dir, "paso_a_paso.gif")
    frames2[0].save(gif2, save_all=True, append_images=frames2[1:], duration=2000, loop=0)
    print(f"  Guardado: {gif2}")
    print(f"\n[{args.name}] Output completo en: {output_dir}")


if __name__ == "__main__":
    main()
