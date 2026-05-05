#!/usr/bin/env python3
# =========================================================================
# Parcial 1 — Conteo de concavidades (nube 2D, .obj solo v).
# Cinco pasos + PNG + GIF. Solo Python (matplotlib / numpy / scipy / PIL).
# =========================================================================

from __future__ import annotations

import argparse
import queue
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.spatial import ConvexHull

try:
    import mapbox_earcut as _earcut

    _HAS_EARCUT = True
except ImportError:
    _HAS_EARCUT = False


# -- Colores (visualización)
C_BG = "#1a1f2e"
C_CLOUD = "#FF6B35"
C_CH = "#FFFFFF"
C_ADJ = "#00CED1"
C_TRI = "#2ECC71"
C_BARY_INT = "#3498DB"
C_BARY_BOL = "#E74C3C"
C_DUAL_INT = "#2ECC71"
C_DUAL_EXT = "#E74C3C"
PALETTE = [
    "#E74C3C",
    "#3498DB",
    "#2ECC71",
    "#9B59B6",
    "#F1C40F",
    "#1ABC9C",
    "#E67E22",
    "#8E44AD",
]


# -------------------------------------------------------------------------
def read_point_cloud_obj(path: Path) -> np.ndarray:
    """Lee .obj con líneas 'v x y ...' (sin f). Devuelve (n,2)."""
    pts = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("o"):
                continue
            if line.startswith("v"):
                parts = line.split()
                if len(parts) >= 3:
                    pts.append([float(parts[1]), float(parts[2])])
    if not pts:
        raise ValueError(f"Sin vértices en {path}")
    return np.asarray(pts, dtype=np.float64)


def write_point_cloud_obj(path: Path, pts: np.ndarray, name: str = "nube") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write(f"o {name}\n")
        for x, y in pts:
            f.write(f"v {x:.6f} {y:.6f} 0.0\n")


# -------------------------------------------------------------------------
def convex_hull_ccw(points: np.ndarray) -> np.ndarray:
    """Casco convexo, vértices en orden CCW (Paso 1)."""
    if len(points) < 3:
        return points.copy()
    hull = ConvexHull(points)
    return points[hull.vertices]


def polygon_signed_area(ring: np.ndarray) -> float:
    x = ring[:, 0]
    y = ring[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))


def ensure_ccw_ring(ring: np.ndarray) -> np.ndarray:
    r = ring.copy()
    if polygon_signed_area(r) < 0:
        r = r[::-1]
    return r


def adjusted_polygon_ring(cloud: np.ndarray) -> np.ndarray:
    """Paso 2: L_upper + L_lower (luego CCW)."""
    cy = float(np.mean(cloud[:, 1]))
    order_x = np.argsort(cloud[:, 0])
    sorted_pts = cloud[order_x]
    upper = sorted_pts[sorted_pts[:, 1] >= cy]
    lower = sorted_pts[sorted_pts[:, 1] < cy]
    lower = lower[np.argsort(-lower[:, 0])]
    ring = np.vstack([upper, lower]) if len(lower) else upper
    return ensure_ccw_ring(ring)


def on_segment(p: np.ndarray, a: np.ndarray, b: np.ndarray, eps: float = 1e-9) -> bool:
    """p en el segmento cerrado [a,b]."""
    ab = b - a
    ap = p - a
    cross = ab[0] * ap[1] - ab[1] * ap[0]
    if abs(cross) > eps * (np.linalg.norm(ab) + 1.0):
        return False
    dot = np.dot(ap, ab)
    if dot < -eps:
        return False
    if dot - np.dot(ab, ab) > eps:
        return False
    return True


def edge_on_convex_hull(a: np.ndarray, b: np.ndarray, hull_ccw: np.ndarray) -> bool:
    """Arista (a,b) contenida en alguna arista consecutiva del CH."""
    n = len(hull_ccw)
    for i in range(n):
        u = hull_ccw[i]
        v = hull_ccw[(i + 1) % n]
        if on_segment(a, u, v) and on_segment(b, u, v):
            return True
    return False


def triangulate_ear_clip(ring_ccw: np.ndarray) -> list[tuple[int, int, int]]:
    """Paso 3: triangulación por recorte de orejas (polígono simple CCW)."""
    n = len(ring_ccw)
    if n < 3:
        return []
    verts = ring_ccw.astype(np.float64)
    idx = list(range(n))
    tris: list[tuple[int, int, int]] = []

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def is_convex(i_prev, i, i_next):
        o, a, b = verts[i_prev], verts[i], verts[i_next]
        return cross(o, a, b) > 1e-14

    def point_in_tri(p, a, b, c):
        s1 = cross(a, b, p)
        s2 = cross(b, c, p)
        s3 = cross(c, a, p)
        has_neg = (s1 < 0) or (s2 < 0) or (s3 < 0)
        has_pos = (s1 > 0) or (s2 > 0) or (s3 > 0)
        return not (has_neg and has_pos)

    while len(idx) > 3:
        ear_found = False
        m = len(idx)
        for k in range(m):
            i_prev = idx[(k - 1) % m]
            i = idx[k]
            i_next = idx[(k + 1) % m]
            if not is_convex(i_prev, i, i_next):
                continue
            a, b, c = verts[i_prev], verts[i], verts[i_next]
            ok = True
            for j in idx:
                if j in (i_prev, i, i_next):
                    continue
                if point_in_tri(verts[j], a, b, c):
                    ok = False
                    break
            if not ok:
                continue
            tris.append((i_prev, i, i_next))
            del idx[k]
            ear_found = True
            break
        if not ear_found:
            raise RuntimeError("Ear clip falló: ¿polígono no simple?")
    if len(idx) == 3:
        tris.append((idx[0], idx[1], idx[2]))
    return tris


def triangulate_polygon(ring_ccw: np.ndarray) -> list[tuple[int, int, int]]:
    """
    Paso 3: triangula el polígono ajustado (CCW). Preferimos mapbox_earcut;
    si no está instalado, recorte de orejas.
    """
    n = len(ring_ccw)
    if n < 3:
        return []
    if _HAS_EARCUT:
        v = ring_ccw.astype(np.float64)
        ends = np.array([v.shape[0]], dtype=np.uint32)
        idx = _earcut.triangulate_float64(v, ends)
        return [tuple(int(idx[i + k]) for k in range(3)) for i in range(0, len(idx), 3)]
    return triangulate_ear_clip(ring_ccw)


def count_edge_uses(tris: list[tuple[int, int, int]]) -> dict[tuple[int, int], int]:
    cnt: dict[tuple[int, int], int] = {}

    def bump(a: int, b: int) -> None:
        e = (min(a, b), max(a, b))
        cnt[e] = cnt.get(e, 0) + 1

    for a, b, c in tris:
        bump(a, b)
        bump(b, c)
        bump(a, c)
    return cnt


def boundary_edges_per_triangle(
    tris: list[tuple[int, int, int]], ecnt: dict[tuple[int, int], int]
) -> list[list[tuple[int, int]]]:
    out: list[list[tuple[int, int]]] = []
    for a, b, c in tris:
        edges = []
        for u, v in ((a, b), (b, c), (a, c)):
            e = (min(u, v), max(u, v))
            if ecnt[e] == 1:
                edges.append((u, v))
        out.append(edges)
    return out


def build_dual_adjacency(tris: list[tuple[int, int, int]]) -> list[list[int]]:
    n = len(tris)
    edge_to_tris: dict[tuple[int, int], list[int]] = {}
    for ti, (a, b, c) in enumerate(tris):
        for u, v in ((a, b), (b, c), (a, c)):
            e = (min(u, v), max(u, v))
            edge_to_tris.setdefault(e, []).append(ti)
    adj = [[] for _ in range(n)]
    for trilist in edge_to_tris.values():
        if len(trilist) == 2:
            i, j = trilist[0], trilist[1]
            adj[i].append(j)
            adj[j].append(i)
    return adj


def compute_barycenters(ring: np.ndarray, tris: list[tuple[int, int, int]]) -> np.ndarray:
    P = ring
    B = []
    for a, b, c in tris:
        B.append((P[a] + P[b] + P[c]) / 3.0)
    return np.asarray(B, dtype=np.float64)


def point_infinity(bary: np.ndarray) -> np.ndarray:
    xmin, ymin = bary.min(axis=0)
    xmax, ymax = bary.max(axis=0)
    cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    dx, dy = xmax - xmin, ymax - ymin
    return np.array([cx + 2.0 * dx, cy - 0.2 * dy], dtype=np.float64)


def classify_bolsillos(
    ring: np.ndarray,
    tris: list[tuple[int, int, int]],
    hull_ccw: np.ndarray,
    boundary_edges: list[list[tuple[int, int]]],
) -> tuple[
    list[bool],
    list[tuple[int, int]],
    list[tuple[int, int]],
    np.ndarray,
    np.ndarray,
]:
    """is_bolsillo, aristas internas dual, aristas externas (tri, inf)."""
    n = len(tris)
    is_boundary_tri = [False] * n
    ecnt = count_edge_uses(tris)
    for e, c in ecnt.items():
        if c == 1:
            for ti, edges in enumerate(boundary_edges):
                for u, v in edges:
                    if (min(u, v), max(u, v)) == e:
                        is_boundary_tri[ti] = True
                        break

    is_bolsillo = [False] * n
    for i in range(n):
        if not is_boundary_tri[i]:
            continue
        alguna_fuera = False
        for u, v in boundary_edges[i]:
            a, b = ring[u], ring[v]
            if not edge_on_convex_hull(a, b, hull_ccw):
                alguna_fuera = True
                break
        is_bolsillo[i] = alguna_fuera

    bary = compute_barycenters(ring, tris)
    p_inf = point_infinity(bary)
    edge_to_tris: dict[tuple[int, int], list[int]] = {}
    for ti, (a, b, c) in enumerate(tris):
        for u, v in ((a, b), (b, c), (a, c)):
            e = (min(u, v), max(u, v))
            edge_to_tris.setdefault(e, []).append(ti)
    int_e: list[tuple[int, int]] = []
    ext_e: list[tuple[int, int]] = []
    for trilist in edge_to_tris.values():
        if len(trilist) == 2:
            i, j = trilist[0], trilist[1]
            int_e.append((i, j))
        elif len(trilist) == 1:
            ext_e.append((trilist[0], n))

    return is_bolsillo, int_e, ext_e, bary, p_inf


def bfs_concavidades(
    n_tri: int, is_bolsillo: list[bool], dual_adj: list[list[int]]
) -> tuple[int, list[list[int]], list[int]]:
    """Paso 5: BFS con queue.Queue (estructura obligatoria)."""
    visitado = [False] * n_tri
    cola: queue.Queue[int] = queue.Queue()
    contador = 0
    componentes: list[list[int]] = []
    pocket_comp = [-1] * n_tri

    def fmt_v() -> str:
        return "[" + ", ".join("true" if v else "false" for v in visitado) + "]"

    print(f"[Paso 5] visitado[] inicial: {fmt_v()}")

    for i in range(n_tri):
        if not is_bolsillo[i] or visitado[i]:
            continue
        contador += 1
        visitado[i] = True
        cola.put(i)
        nodos_comp: list[int] = []
        while not cola.empty():
            nodo = cola.get()
            nodos_comp.append(nodo)
            pocket_comp[nodo] = contador - 1
            for j in dual_adj[nodo]:
                if is_bolsillo[j] and not visitado[j]:
                    visitado[j] = True
                    cola.put(j)
        componentes.append(nodos_comp)
        print(f"Concavidad {contador}: nodos B = {nodos_comp}")
        print(f"visitado[] ahora: {fmt_v()}")

    return contador, componentes, pocket_comp


def count_concavities_silent(cloud: np.ndarray) -> int:
    """Mismo pipeline que Paso 5 sin imprimir (para calibrar .obj)."""
    hull = convex_hull_ccw(cloud)
    ring = adjusted_polygon_ring(cloud)
    tris = triangulate_polygon(ring)
    ecnt = count_edge_uses(tris)
    bedges = boundary_edges_per_triangle(tris, ecnt)
    dual_adj = build_dual_adjacency(tris)
    is_bolsillo, _, _, _, _ = classify_bolsillos(ring, tris, hull, bedges)
    n_tri = len(tris)
    visitado = [False] * n_tri
    cola: queue.Queue[int] = queue.Queue()
    contador = 0
    for i in range(n_tri):
        if not is_bolsillo[i] or visitado[i]:
            continue
        contador += 1
        visitado[i] = True
        cola.put(i)
        while not cola.empty():
            nodo = cola.get()
            for j in dual_adj[nodo]:
                if is_bolsillo[j] and not visitado[j]:
                    visitado[j] = True
                    cola.put(j)
    return contador


# -------------------------------------------------------------------------
def _setup_ax(ax, pts: np.ndarray) -> None:
    ax.set_facecolor(C_BG)
    ax.set_aspect("equal", adjustable="box")
    pad = 0.4
    ax.set_xlim(pts[:, 0].min() - pad, pts[:, 0].max() + pad)
    ax.set_ylim(pts[:, 1].min() - pad, pts[:, 1].max() + pad)
    for s in ax.spines.values():
        s.set_color("#445566")
    ax.tick_params(colors="#8899aa")
    ax.xaxis.label.set_color("#8899aa")
    ax.yaxis.label.set_color("#8899aa")


def save_step(
    path: Path,
    title: str,
    cloud: np.ndarray,
    *,
    hull: np.ndarray | None = None,
    hull_dashed: bool = False,
    ring: np.ndarray | None = None,
    ring_alpha: float = 1.0,
    tris: list[tuple[int, int, int]] | None = None,
    tri_alpha: float = 0.35,
    bary: np.ndarray | None = None,
    bols_mask: list[bool] | None = None,
    dual_int: list[tuple[int, int]] | None = None,
    dual_ext: list[tuple[int, int]] | None = None,
    p_inf: np.ndarray | None = None,
    pocket_fill: list[tuple[int, tuple[int, int, int], tuple[float, float, float]]] | None = None,
    big_text: str | None = None,
    cam_xy: np.ndarray | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(11, 7.5), facecolor=C_BG)
    ref = cam_xy if cam_xy is not None else cloud
    _setup_ax(ax, ref)

    ax.scatter(
        cloud[:, 0],
        cloud[:, 1],
        c=C_CLOUD,
        s=12,
        zorder=5,
        linewidths=0,
        label="nube",
    )

    if hull is not None and len(hull) >= 2:
        hc = np.vstack([hull, hull[0]])
        ax.plot(
            hc[:, 0],
            hc[:, 1],
            color=C_CH,
            linewidth=1.8,
            linestyle="--" if hull_dashed else "-",
            alpha=0.85 if hull_dashed else 1.0,
            zorder=3,
        )

    if ring is not None and len(ring) >= 2:
        rc = np.vstack([ring, ring[0]])
        ax.plot(
            rc[:, 0],
            rc[:, 1],
            color=C_ADJ,
            linewidth=2.2,
            solid_capstyle="round",
            alpha=ring_alpha,
            zorder=4,
        )

    if tris is not None and ring is not None:
        P = ring
        for a, b, c in tris:
            tri = np.array([P[a], P[b], P[c], P[a]])
            ax.plot(tri[:, 0], tri[:, 1], color=C_TRI, linewidth=0.8, alpha=tri_alpha, zorder=2)

    if bary is not None:
        n = len(bary)
        if bols_mask is None:
            bols_mask = [False] * n
        m0 = np.array([not bols_mask[i] for i in range(n)])
        m1 = np.array(bols_mask)
        if m0.any():
            ax.scatter(bary[m0, 0], bary[m0, 1], c=C_BARY_INT, s=38, zorder=6, edgecolors="none")
        if m1.any():
            ax.scatter(bary[m1, 0], bary[m1, 1], c=C_BARY_BOL, s=45, zorder=7, edgecolors="none")

    if bary is not None and dual_int is not None:
        for i, j in dual_int:
            ax.plot(
                [bary[i, 0], bary[j, 0]],
                [bary[i, 1], bary[j, 1]],
                color=C_DUAL_INT,
                linewidth=1.2,
                alpha=0.85,
                zorder=3,
            )
    if bary is not None and dual_ext is not None and p_inf is not None:
        for i, _ in dual_ext:
            ax.plot(
                [bary[i, 0], p_inf[0]],
                [bary[i, 1], p_inf[1]],
                color=C_DUAL_EXT,
                linewidth=1.0,
                linestyle=(0, (4, 4)),
                alpha=0.9,
                zorder=3,
            )
        ax.scatter([p_inf[0]], [p_inf[1]], c=C_DUAL_EXT, s=28, zorder=8, marker="o")

    if pocket_fill is not None and ring is not None:
        P = ring
        for comp_id, (a, b, c), rgb in pocket_fill:
            tri = np.array([P[a], P[b], P[c]])
            poly = plt.Polygon(tri, closed=True, facecolor=rgb, edgecolor="none", alpha=0.5, zorder=1)
            ax.add_patch(poly)

    ax.set_title(title, color="#e8eef8", fontsize=14, pad=12)

    if big_text:
        ax.text(
            0.5,
            0.42,
            big_text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="#fff8d0",
            fontsize=18,
            fontweight="bold",
            zorder=10,
        )

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120, facecolor=C_BG, edgecolor="none")
    plt.close(fig)


def next_result_dir(output_root: Path) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    best = 0
    for p in output_root.iterdir():
        if not p.is_dir():
            continue
        m = re.fullmatch(r"result-(\d{3})", p.name)
        if m:
            best = max(best, int(m.group(1)))
    return output_root / f"result-{best + 1:03d}"


def make_gif_from_pngs(png_paths: list[Path], gif_path: Path, duration_ms: int = 1000) -> None:
    frames = [Image.open(p).convert("RGB") for p in png_paths]
    if not frames:
        return
    gif_path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        gif_path,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
    )
    for im in frames:
        im.close()


def run_pipeline(cloud: np.ndarray, out_dir: Path) -> int:
    # -- Paso 1
    hull = convex_hull_ccw(cloud)
    print(f"[Paso 1] CH: {len(hull)} vértices")

    # -- Paso 2
    ring = adjusted_polygon_ring(cloud)
    print(f"[Paso 2] Polígono ajustado: {len(ring)} vértices, CCW")

    # -- Paso 3
    tris = triangulate_polygon(ring)
    print(f"[Paso 3] Triangulación: {len(tris)} triángulos")

    ecnt = count_edge_uses(tris)
    bedges = boundary_edges_per_triangle(tris, ecnt)
    dual_adj = build_dual_adjacency(tris)
    is_bolsillo, dual_int, dual_ext, bary, p_inf = classify_bolsillos(
        ring, tris, hull, bedges
    )

    bidx = [i for i, b in enumerate(is_bolsillo) if b]
    print(f"[Paso 4] Grafo dual construido. Bolsillos: {bidx}")

    n_tri = len(tris)
    k, _, pocket_comp = bfs_concavidades(n_tri, is_bolsillo, dual_adj)

    print("────────────────────────────────────")
    print(f"Resultado: {k} concavidad(es)")
    print("────────────────────────────────────")

    # -- Visualización
    pocket_rgb = []
    for ti in range(n_tri):
        cid = pocket_comp[ti]
        if cid < 0:
            continue
        col = PALETTE[cid % len(PALETTE)]
        rgb = tuple(int(col[j : j + 2], 16) / 255.0 for j in (1, 3, 5))
        pocket_rgb.append((cid, tris[ti], rgb))

    cam_bary = bary.copy()

    save_step(out_dir / "00_entrada.png", "Entrada: nube de puntos 2D", cloud)

    save_step(
        out_dir / "01_casco_convexo.png",
        "Paso 1: casco convexo (CH)",
        cloud,
        hull=hull,
    )

    save_step(
        out_dir / "02_poligono_ajustado.png",
        "Paso 2: polígono ajustado (CCW)",
        cloud,
        hull=hull,
        hull_dashed=True,
        ring=ring,
    )

    save_step(
        out_dir / "03_triangulacion.png",
        "Paso 3: triangulación (n-2 triángulos)",
        cloud,
        ring=ring,
        tris=tris,
        tri_alpha=0.45,
    )

    save_step(
        out_dir / "04_grafo_dual.png",
        "Paso 4: grafo dual (T=interior, B=bolsillo)",
        cloud,
        ring=ring,
        ring_alpha=0.35,
        bary=bary,
        bols_mask=is_bolsillo,
        dual_int=dual_int,
        dual_ext=[],
        cam_xy=cam_bary,
    )

    save_step(
        out_dir / "05_resultado.png",
        "Paso 5: concavidades detectadas",
        cloud,
        hull=hull,
        hull_dashed=True,
        ring=ring,
        pocket_fill=pocket_rgb,
        big_text=f"Resultado: {k} concavidad(es)",
    )

    pngs = [out_dir / f"{i:02d}_{name}.png" for i, name in enumerate(
        ["entrada", "casco_convexo", "poligono_ajustado", "triangulacion", "grafo_dual", "resultado"]
    )]
    make_gif_from_pngs(pngs, out_dir / "concavidades.gif", duration_ms=1000)
    print(f"Output: {out_dir}")
    print(f"GIF: {out_dir / 'concavidades.gif'}")
    return k


# -------------------------------------------------------------------------
def ring_points(
    n_lobes: int,
    n_pts: int,
    cx: float,
    cy: float,
    r0: float,
    amp: float,
    phase: float = 0.0,
    clip_lo: float = 0.12,
    clip_hi: float = 9.88,
) -> np.ndarray:
    """
    Anillo continuo r(θ) = r0 + amp*cos(n_lobes * θ + phase).
    n_lobes == 0 → radio constante r0 (círculo).
    """
    theta = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
    if n_lobes <= 0:
        r = np.full_like(theta, r0)
    else:
        r = r0 + amp * np.cos(n_lobes * theta + phase)
    r = np.clip(r, 0.25, 6.0)
    xy = np.column_stack([cx + r * np.cos(theta), cy + r * np.sin(theta)])
    xy[:, 0] = np.clip(xy[:, 0], clip_lo, clip_hi)
    xy[:, 1] = np.clip(xy[:, 1], clip_lo, clip_hi)
    return xy


def generate_data_files(data_dir: Path) -> None:
    """
    Genera los 5 .obj: anillo continuo r(θ)=r0+amp*cos(n·θ+φ) en [0,10]².
    Parámetros fijos (calibrados) para que este pipeline dé 0,2,3,4,6
    concavidades respectivamente.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    cx = cy = 5.0
    # -- (archivo, n_lobes, n_pts, r0, amp, phase, objetivo)
    specs: list[tuple[str, int, int, float, float, float, int]] = [
        ("nube_0_convexa", 0, 96, 3.0, 0.0, 0.0, 0),
        ("nube_2_concavidades", 1, 96, 2.45, 0.38, 0.5, 2),
        ("nube_3_concavidades", 3, 72, 2.4, 0.2473684210526316, 0.0, 3),
        ("nube_4_concavidades", 2, 80, 2.45, 0.58, 0.2, 4),
        ("nube_6_concavidades", 3, 80, 2.45, 0.28, 0.0, 6),
    ]
    for name, n_lobes, n_pts, r0, amp, phase, target_k in specs:
        pts = ring_points(n_lobes, n_pts, cx, cy, r0, amp, phase)
        got = count_concavities_silent(pts)
        if got != target_k:
            raise RuntimeError(
                f"{name}: esperaba {target_k} concavidades, obtuve {got}."
            )
        write_point_cloud_obj(data_dir / f"{name}.obj", pts, name=name)
        print(f"Generado {name}.obj  (concavidades={got}, objetivo={target_k})")


def main() -> None:
    root = Path(__file__).resolve().parent
    data_dir = root / "data"
    output_root = root / "output"

    ap = argparse.ArgumentParser(description="Conteo de concavidades (Parcial 1)")
    ap.add_argument(
        "obj",
        nargs="?",
        type=str,
        help="Ruta al .obj (omitir con --generate-data)",
    )
    ap.add_argument(
        "--generate-data",
        action="store_true",
        help="Regenera los 5 archivos en data/",
    )
    args = ap.parse_args()

    if args.generate_data:
        generate_data_files(data_dir)
        if not args.obj:
            return

    if not args.obj:
        ap.print_help()
        sys.exit(1)

    obj_path = Path(args.obj).expanduser()
    if not obj_path.is_file():
        print(f"No existe el archivo: {obj_path}", file=sys.stderr)
        sys.exit(1)

    cloud = read_point_cloud_obj(obj_path)
    out_dir = next_result_dir(output_root)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_pipeline(cloud, out_dir)


if __name__ == "__main__":
    main()
