#!/usr/bin/env python3
"""
generar_diagramas.py
Taller 3 - Geometria Computacional - Pontificia Universidad Javeriana
Genera diagramas explicativos del algoritmo de simplificacion de heightmap
y un GIF animado con los 7 pasos del proceso.

Dependencias: matplotlib numpy Pillow
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image

# --- Configuracion global -------------------------------------------------
OUT = 'diagramas'
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
})

C = {
    'azul':        '#378ADD',
    'azul_osc':    '#1F4E79',
    'azul_clar':   '#D5E8F0',
    'verde':       '#1D9E75',
    'verde_osc':   '#0F6E56',
    'verde_bg':    '#EAF3DE',
    'rojo':        '#D85A30',
    'rojo_osc':    '#993C1D',
    'rojo_bg':     '#FAECE7',
    'naranja':     '#BA7517',
    'naranja_bg':  '#FAEEDA',
    'gris':        '#888780',
    'gris_clar':   '#B4B2A9',
    'negro':       '#2C2C2A',
}

DPI = 150


def guardar(nombre, fig=None):
    path = os.path.join(OUT, nombre)
    if fig is None:
        fig = plt.gcf()
    fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  OK: {nombre}')
    return path


# --- Helpers de diagramas ------------------------------------------------
def estrella(ax, center, radio, n=6, colores=None, etiquetas=None,
             fill_colors=None, stroke_colors=None, titulo_p=None,
             mostrar_p=True, punto_extra=None, punto_extra_color='#D85A30',
             alpha_fill=0.85):
    """Dibuja una estrella de n triangulos alrededor de un punto central."""
    cx, cy = center
    angulos = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, n + 1)
    pts = [(cx + radio * np.cos(a), cy + radio * np.sin(a)) for a in angulos]

    if fill_colors is None:
        fill_colors = [C['azul_clar']] * n
    if stroke_colors is None:
        stroke_colors = [C['azul']] * n

    for i in range(n):
        tri = plt.Polygon(
            [center, pts[i], pts[i + 1]],
            facecolor=fill_colors[i % len(fill_colors)],
            edgecolor=stroke_colors[i % len(stroke_colors)],
            linewidth=1.2, alpha=alpha_fill
        )
        ax.add_patch(tri)

    if mostrar_p:
        ax.scatter(*center, color=C['negro'], s=320, zorder=12)
        ax.text(cx, cy, 'p', color='white', ha='center', va='center',
                fontsize=12, fontweight='bold', zorder=13)
        if titulo_p:
            ax.text(cx, cy - radio * 0.22, titulo_p, ha='center',
                    fontsize=9, color=C['negro'], zorder=14)

    if punto_extra:
        ax.scatter(*punto_extra, color=punto_extra_color, s=200,
                   zorder=11, edgecolors='white', linewidths=1.5)

    vecinos = []
    for i in range(n):
        col = colores[i] if colores else C['azul']
        ax.scatter(pts[i][0], pts[i][1], color=col, s=160, zorder=10,
                   edgecolors='white', linewidths=1.5)
        if etiquetas:
            off = 1.22
            ax.text(pts[i][0] * off + cx * (1 - off),
                    pts[i][1] * off + cy * (1 - off),
                    etiquetas[i], ha='center', va='center',
                    fontsize=10, fontweight='bold', color=col)
        vecinos.append(pts[i])

    return vecinos, pts


# ========================================================================
# DIAGRAMA 1 - Criterio de planitud
# ========================================================================
def d1_criterio():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Criterio de Planitud: icuando sobra un vertice?',
                 fontsize=15, fontweight='bold', color=C['azul_osc'])

    center = (0, 0)
    radio = 1.0
    ang = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    px = [radio * np.cos(a) for a in ang]
    py = [radio * np.sin(a) for a in ang]

    # - Caso plano (izquierda) -
    for i in range(5):
        tri = plt.Polygon(
            [center, (px[i], py[i]), (px[(i + 1) % 5], py[(i + 1) % 5])],
            facecolor=C['verde_bg'], edgecolor=C['verde'], lw=1.3, alpha=0.9
        )
        ax1.add_patch(tri)
    ax1.scatter(0, 0, color=C['verde_osc'], s=320, zorder=10)
    ax1.text(0, 0, 'p', color='white', ha='center', va='center',
             fontsize=12, fontweight='bold', zorder=11)
    ax1.text(0, -0.22, 'z = 0.51', ha='center', fontsize=9,
             color=C['verde_osc'])
    vals_p = [0.50, 0.52, 0.51, 0.50, 0.52]
    for i in range(5):
        ax1.scatter(px[i], py[i], color=C['verde'], s=140, zorder=9,
                    edgecolors='white', lw=1.5)
        ox = 1.26 * px[i]
        oy = 1.26 * py[i]
        ax1.text(ox, oy, f'z={vals_p[i]}', ha='center', fontsize=9,
                 color=C['verde_osc'])
    caja1 = FancyBboxPatch((-1.15, -1.65), 2.3, 0.50,
                           boxstyle='round,pad=0.06',
                           facecolor=C['verde_bg'],
                           edgecolor=C['verde'], lw=1.8)
    ax1.add_patch(caja1)
    ax1.text(0, -1.40, 'error(p) = 0.008 < e (0.05)', ha='center',
             fontsize=11, fontweight='bold', color=C['verde_osc'])
    ax1.text(0, -1.58, '->  ELIMINAR', ha='center', fontsize=12,
             fontweight='bold', color=C['verde_osc'])
    ax1.set_title('Zona plana', fontsize=13, fontweight='bold',
                  color=C['verde'], pad=10)
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.8, 1.4)
    ax1.axis('off')

    # - Caso borde (derecha) -
    vals_b = [0.10, 0.90, 0.85, 0.12, 0.88]
    for i in range(5):
        tri = plt.Polygon(
            [center, (px[i], py[i]), (px[(i + 1) % 5], py[(i + 1) % 5])],
            facecolor=C['rojo_bg'], edgecolor=C['rojo'], lw=1.3, alpha=0.9
        )
        ax2.add_patch(tri)
    ax2.scatter(0, 0, color='#4A1B0C', s=320, zorder=10)
    ax2.text(0, 0, 'p', color='white', ha='center', va='center',
             fontsize=12, fontweight='bold', zorder=11)
    ax2.text(0, -0.22, 'z = 0.50', ha='center', fontsize=9, color='#4A1B0C')
    for i in range(5):
        ax2.scatter(px[i], py[i], color=C['rojo'], s=140, zorder=9,
                    edgecolors='white', lw=1.5)
        ox = 1.26 * px[i]
        oy = 1.26 * py[i]
        ax2.text(ox, oy, f'z={vals_b[i]}', ha='center', fontsize=9,
                 color=C['rojo_osc'])
    caja2 = FancyBboxPatch((-1.15, -1.65), 2.3, 0.50,
                           boxstyle='round,pad=0.06',
                           facecolor=C['rojo_bg'],
                           edgecolor=C['rojo'], lw=1.8)
    ax2.add_patch(caja2)
    ax2.text(0, -1.40, 'error(p) = 0.38 >= e (0.05)', ha='center',
             fontsize=11, fontweight='bold', color=C['rojo_osc'])
    ax2.text(0, -1.58, '->  CONSERVAR', ha='center', fontsize=12,
             fontweight='bold', color=C['rojo_osc'])
    ax2.set_title('Zona de borde', fontsize=13, fontweight='bold',
                  color=C['rojo'], pad=10)
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(-1.8, 1.4)
    ax2.axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    return guardar('01_criterio_planitud.png', fig)


# ========================================================================
# DIAGRAMA 2 - Face_circulator (DCEL)
# ========================================================================
def d2_circulator():
    fig, ax = plt.subplots(figsize=(9, 9))
    fig.suptitle('Face_circulator: Recorrido de la Estrella de p\n'
                 '(DCEL half-edge twin->next)',
                 fontsize=14, fontweight='bold', color=C['azul_osc'])

    fills = [C['azul_clar'], C['verde_bg'], C['naranja_bg'],
             '#EEEDFE', C['rojo_bg'], '#E1F5EE']
    strks = [C['azul'], C['verde'], C['naranja'],
             '#7F77DD', C['rojo'], C['verde_osc']]
    cols_q = strks[:]
    etiqs = [f'q{i}' for i in range(6)]

    vecinos, pts = estrella(
        ax, center=(0, 0), radio=1.1, n=6,
        colores=cols_q, etiquetas=etiqs,
        fill_colors=fills, stroke_colors=strks,
        titulo_p='(center)', alpha_fill=0.82
    )

    # Etiquetas de triangulos
    ang6 = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, 7)
    for i in range(6):
        mid_a = (ang6[i] + ang6[i + 1]) / 2
        ax.text(0.48 * np.cos(mid_a), 0.48 * np.sin(mid_a),
                f'T{i}', ha='center', va='center',
                fontsize=10, fontweight='bold', color=strks[i])

    # Half-edge inicial (p -> q0)
    ax.annotate('', xy=(pts[0][0] * 0.82, pts[0][1] * 0.82),
                xytext=(0.1, 0.05),
                arrowprops=dict(arrowstyle='->', color=C['azul'],
                                lw=2.5, mutation_scale=18))
    ax.text(0.52, 0.72, 'half-edge\ninicial', fontsize=9,
            color=C['azul_osc'], ha='center')

    # Arco de rotacion
    theta = np.linspace(0.3, 5.5, 80)
    rx, ry = 0.38 * np.cos(theta), 0.38 * np.sin(theta)
    ax.plot(rx, ry, color=C['negro'], lw=1.8, linestyle='--')
    ax.annotate('', xy=(rx[-1], ry[-1]),
                xytext=(rx[-2], ry[-2]),
                arrowprops=dict(arrowstyle='->', color=C['negro'], lw=1.5))
    ax.text(-0.52, 0.22, 'girar <->', fontsize=14, color=C['negro'])

    ax.text(0, -1.52,
            'Face_circulator recorre T0->T1->...->T5->T0',
            ha='center', fontsize=10, style='italic', color=C['gris'])
    ax.text(0, -1.72,
            'Cada cara tiene 3 vertices: p + 2 vecinos -> std::set elimina duplicados',
            ha='center', fontsize=10, color='#5F5E5A')

    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.85, 1.55)
    ax.axis('off')
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    return guardar('02_face_circulator.png', fig)


# ========================================================================
# DIAGRAMA 3 - Min-heap (Dijkstra)
# ========================================================================
def d3_heap():
    fig, ax = plt.subplots(figsize=(12, 7))
    fig.suptitle('priority_queue - Menor error de planitud sale primero\n'
                 '(Dijkstra: siempre procesar el mas prescindible)',
                 fontsize=13, fontweight='bold', color=C['azul_osc'])

    nodos = [
        ('r',   5.0, 4.2, 'v42',  'err=0.002', '#085041', C['verde_bg'],  '* sale primero'),
        ('l1',  2.5, 2.8, 'v17',  'err=0.018', '#27500A', C['verde_bg'],  ''),
        ('r1',  7.5, 2.8, 'v99',  'err=0.031', '#3B6D11', C['verde_bg'],  ''),
        ('l2a', 1.0, 1.2, 'v55',  'err=0.045', '#633806', C['naranja_bg'], ''),
        ('l2b', 4.0, 1.2, 'v23',  'err=0.052', '#854F0B', C['naranja_bg'], ''),
        ('r2a', 6.0, 1.2, 'v3',   'err=0.071', '#993C1D', C['rojo_bg'],   ''),
        ('r2b', 9.0, 1.2, 'v8',   'err=0.089', '#4A1B0C', C['rojo_bg'],   ''),
    ]
    pos = {n[0]: (n[1], n[2]) for n in nodos}
    edges = [('r', 'l1'), ('r', 'r1'), ('l1', 'l2a'), ('l1', 'l2b'),
             ('r1', 'r2a'), ('r1', 'r2b')]

    for u, v in edges:
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=C['gris_clar'], lw=1.5, zorder=1)

    for key, x, y, lbl, err, tc, fc, extra in nodos:
        caja = FancyBboxPatch((x - 0.9, y - 0.5), 1.8, 1.0,
                              boxstyle='round,pad=0.12',
                              facecolor=fc, edgecolor=tc, lw=1.8, zorder=2)
        ax.add_patch(caja)
        ax.text(x, y + 0.15, lbl, ha='center', va='center',
                fontsize=11, fontweight='bold', color=tc, zorder=3)
        ax.text(x, y - 0.18, err, ha='center', va='center',
                fontsize=9.5, color=tc, zorder=3)
        if extra:
            ax.text(x, y + 0.72, extra, ha='center', fontsize=10,
                    fontweight='bold', color=C['verde_osc'],
                    bbox=dict(facecolor='white', edgecolor=C['verde_osc'],
                              boxstyle='round,pad=0.25', lw=1.2))

    # Flecha "pop"
    ax.annotate('heap.pop() -> procesar v42 primero',
                xy=(4.1, 4.2), xytext=(0.8, 5.1),
                fontsize=11, fontweight='bold', color=C['verde_osc'],
                arrowprops=dict(arrowstyle='->', color=C['verde_osc'], lw=2.2))

    ax.set_xlim(0, 10.5)
    ax.set_ylim(0.3, 6.0)
    ax.axis('off')
    plt.tight_layout(rect=[0, 0, 1, 0.91])
    return guardar('03_min_heap.png', fig)


# ========================================================================
# DIAGRAMA 4 - T.remove(p) y re-triangulacion
# ========================================================================
def d4_remove():
    fig, (ax1, axm, ax2) = plt.subplots(
        1, 3, figsize=(16, 7),
        gridspec_kw={'width_ratios': [5, 1, 5]}
    )
    fig.suptitle('Efecto de T.remove(p) en la Triangulacion',
                 fontsize=14, fontweight='bold', color=C['azul_osc'])

    ang6 = np.linspace(np.pi / 6, np.pi / 6 + 2 * np.pi, 7)
    pts6 = [(np.cos(a), np.sin(a)) for a in ang6]
    fills6 = [C['azul_clar'], C['verde_bg'], C['naranja_bg'],
              '#EEEDFE', C['rojo_bg'], '#E1F5EE']
    strks6 = [C['azul'], C['verde'], C['naranja'],
              '#7F77DD', C['rojo'], C['verde_osc']]
    cols6 = strks6[:]
    etiqs6 = [f'q{i}' for i in range(6)]

    # ANTES
    for i in range(6):
        ax1.add_patch(plt.Polygon(
            [(0, 0), pts6[i], pts6[i + 1]],
            facecolor=fills6[i], edgecolor=strks6[i], lw=1.3, alpha=0.88
        ))
    ax1.scatter(0, 0, color=C['negro'], s=320, zorder=10)
    ax1.text(0, 0, 'p', color='white', ha='center', va='center',
             fontsize=13, fontweight='bold', zorder=11)
    for i in range(6):
        ax1.scatter(pts6[i][0], pts6[i][1], color=cols6[i], s=160,
                    zorder=10, edgecolors='white', lw=1.5)
        ax1.text(pts6[i][0] * 1.22, pts6[i][1] * 1.22, etiqs6[i],
                 ha='center', fontsize=10, fontweight='bold', color=cols6[i])
    ax1.set_title('ANTES\n6 triangulos incidentes a p',
                  fontsize=12, fontweight='bold', color=C['negro'], pad=8)
    ax1.set_xlim(-1.4, 1.4)
    ax1.set_ylim(-1.5, 1.45)
    ax1.axis('off')

    # Flecha central
    axm.text(0.5, 0.55, 'T.remove(p)\n->', ha='center', va='center',
             fontsize=13, fontweight='bold', color=C['azul_osc'],
             transform=axm.transAxes)
    axm.axis('off')

    # DESPUES
    new_tris = [[0, 1, 2], [0, 2, 3], [0, 3, 5], [3, 4, 5]]
    fill4 = [C['azul_clar'], C['verde_bg'], C['naranja_bg'], '#EEEDFE']
    strk4 = [C['azul'], C['verde'], C['naranja'], '#7F77DD']
    for tri, fc, sc in zip(new_tris, fill4, strk4):
        ax2.add_patch(plt.Polygon(
            [pts6[j] for j in tri],
            facecolor=fc, edgecolor=sc, lw=1.3, alpha=0.88
        ))
    for i in range(6):
        ax2.scatter(pts6[i][0], pts6[i][1], color=cols6[i], s=160,
                    zorder=10, edgecolors='white', lw=1.5)
        ax2.text(pts6[i][0] * 1.22, pts6[i][1] * 1.22, etiqs6[i],
                 ha='center', fontsize=10, fontweight='bold', color=cols6[i])
    ax2.text(0, 0, 'x', color=C['rojo'], ha='center', va='center',
             fontsize=36, fontweight='bold', zorder=10)
    ax2.set_title('DESPUES\nCGAL re-triangula el hueco (Delaunay)',
                  fontsize=12, fontweight='bold', color=C['verde_osc'], pad=8)
    ax2.set_xlim(-1.4, 1.4)
    ax2.set_ylim(-1.5, 1.45)
    ax2.axis('off')

    fig.text(0.5, 0.01,
             'La malla sigue siendo Delaunay - CGAL garantiza validez automaticamente',
             ha='center', fontsize=10, style='italic', color=C['gris'])
    plt.tight_layout(rect=[0, 0.05, 1, 0.92])
    return guardar('04_remove_retrigulacion.png', fig)


# ========================================================================
# DIAGRAMA 5 - Lazy deletion
# ========================================================================
def d5_lazy():
    fig, ax = plt.subplots(figsize=(14, 6))
    fig.suptitle('Lazy Deletion: deteccion de entradas obsoletas en el heap',
                 fontsize=14, fontweight='bold', color=C['azul_osc'])

    pasos = [
        ('Momento 1\nHeap inicial',
         [('qi  err=0.020', C['verde_bg'], C['verde_osc'], 'entrada valida')],
         ''),
        ('Momento 2\nDespues de T.remove(p)',
         [('qi  err=0.020', C['rojo_bg'], C['rojo_osc'], 'OBSOLETA'),
          ('qi  err=0.041', C['verde_bg'], C['verde_osc'], 'NUEVA (recalculada)')],
         'El entorno de qi cambio\nsu error real es 0.041'),
        ('Momento 3\nSale del heap',
         [('qi  err=0.020', C['rojo_bg'], C['rojo_osc'], 'sale primero\n(menor error)')],
         'iqi en visitados? -> skip\nierr != error_map[qi]? -> skip\n(entrada obsoleta)'),
    ]

    xs = [1.2, 5.2, 9.5]
    for j, (titulo, entradas, nota) in enumerate(pasos):
        x = xs[j]
        ax.text(x, 5.4, titulo, ha='center', fontsize=11,
                fontweight='bold', color=C['azul_osc'])
        rect_bg = FancyBboxPatch((x - 1.1, 1.2), 2.2, 3.8,
                                 boxstyle='round,pad=0.1',
                                 facecolor='#F8F8F8',
                                 edgecolor=C['gris_clar'], lw=1)
        ax.add_patch(rect_bg)
        for k, (texto, fc, tc, sublbl) in enumerate(entradas):
            y = 4.0 - k * 1.6
            caja = FancyBboxPatch((x - 0.95, y - 0.42), 1.9, 0.84,
                                  boxstyle='round,pad=0.08',
                                  facecolor=fc, edgecolor=tc, lw=1.5)
            ax.add_patch(caja)
            ax.text(x, y + 0.10, texto, ha='center', va='center',
                    fontsize=10, fontweight='bold', color=tc)
            ax.text(x, y - 0.18, sublbl, ha='center', va='center',
                    fontsize=8.5, color=tc, style='italic')
        if nota:
            ax.text(x, 1.0, nota, ha='center', fontsize=9,
                    color=C['naranja'], style='italic', va='top')

    # Flechas entre momentos
    for xa, xb in [(2.3, 4.2), (6.3, 8.4)]:
        ax.annotate('', xy=(xb, 3.2), xytext=(xa, 3.2),
                    arrowprops=dict(arrowstyle='->', color=C['gris'],
                                    lw=2, mutation_scale=18))

    ax.set_xlim(0, 12)
    ax.set_ylim(0.5, 6.0)
    ax.axis('off')
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    return guardar('05_lazy_deletion.png', fig)


# ========================================================================
# DIAGRAMA 6 - Vecindad de orden k
# ========================================================================
def d6_vecindad():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    fig.suptitle('Vecindad de Orden k - Expansion por Anillos',
                 fontsize=14, fontweight='bold', color=C['azul_osc'])

    fills1 = [C['azul_clar']] * 6
    strks1 = [C['azul']] * 6
    cols1 = [C['azul']] * 6
    etiqs1 = [f'q{i}' for i in range(6)]

    # Orden 1
    estrella(ax1, (0, 0), 1.0, 6, cols1, etiqs1, fills1, strks1)
    ax1.set_title('Orden 1 - Estrella directa de p\n'
                  'Face_circulator sobre p',
                  fontsize=12, fontweight='bold', color=C['azul'], pad=8)
    ax1.text(0, -1.35, 'std::set recolecta 6 vecinos unicos', ha='center',
             fontsize=10, color=C['gris'])
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.55, 1.45)
    ax1.axis('off')

    # Orden 2 - anillo exterior
    ang_ext = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, 13)
    px_ext = [1.9 * np.cos(a) for a in ang_ext]
    py_ext = [1.9 * np.sin(a) for a in ang_ext]

    ang6 = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, 7)
    px6 = [1.0 * np.cos(a) for a in ang6]
    py6 = [1.0 * np.sin(a) for a in ang6]

    # Triangulos anillo 2 (gris)
    for i in range(12):
        tri = plt.Polygon(
            [(px_ext[i], py_ext[i]),
             (px_ext[i + 1], py_ext[i + 1]),
             (px6[i // 2], py6[i // 2])],
            facecolor='#F1EFE8', edgecolor=C['gris_clar'], lw=0.8, alpha=0.7
        )
        ax2.add_patch(tri)

    # Triangulos anillo 1 (azul)
    for i in range(6):
        tri = plt.Polygon(
            [(0, 0), (px6[i], py6[i]), (px6[i + 1], py6[i + 1])],
            facecolor=C['azul_clar'], edgecolor=C['azul'], lw=1.2, alpha=0.88
        )
        ax2.add_patch(tri)

    # Vertices
    ax2.scatter(0, 0, color=C['negro'], s=320, zorder=12)
    ax2.text(0, 0, 'p', color='white', ha='center', va='center',
             fontsize=12, fontweight='bold', zorder=13)
    for i in range(6):
        ax2.scatter(px6[i], py6[i], color=C['azul'], s=160, zorder=11,
                    edgecolors='white', lw=1.5)
        ax2.text(px6[i] * 1.22, py6[i] * 1.22, f'q{i}', ha='center',
                 fontsize=9, fontweight='bold', color=C['azul'])
    for i in range(0, 12, 2):
        ax2.scatter(px_ext[i], py_ext[i], color=C['gris'], s=100, zorder=10,
                    edgecolors='white', lw=1.2)

    # Leyenda
    leyenda = [
        mpatches.Patch(facecolor=C['azul_clar'], edgecolor=C['azul'], label='Orden 1 (estrella de p)'),
        mpatches.Patch(facecolor='#F1EFE8', edgecolor=C['gris_clar'], label='Orden 2 (estrella de vecinos)'),
    ]
    ax2.legend(handles=leyenda, loc='upper right', fontsize=9,
               framealpha=0.9)
    ax2.set_title('Orden 2 - Segundo anillo\n'
                  'Face_circulator sobre cada qi',
                  fontsize=12, fontweight='bold', color=C['gris'], pad=8)
    ax2.text(0, -2.3, 'El mismo std::set acumula ambos anillos sin duplicados',
             ha='center', fontsize=10, color=C['gris'])
    ax2.set_xlim(-2.4, 2.4)
    ax2.set_ylim(-2.5, 2.4)
    ax2.axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    return guardar('06_vecindad_orden_k.png', fig)


# ========================================================================
# DIAGRAMA 7 - Resultados por imagen
# ========================================================================
def d7_resultados():
    nombres = ['SRTM_256', 'SRTM_512', 'NormalMap\n-second',
               'NormalMap', 'SRTM_1024', 'SRTM_2048']
    antes = [32768, 131072, 196608, 262144, 524288, 2097152]
    despues = [4840, 14190, 22352, 1647, 37188, 81758]
    perc = [85.2, 89.2, 88.6, 99.4, 92.9, 96.1]

    fig, ax = plt.subplots(figsize=(13, 7.5))
    fig.suptitle('Vertices antes vs despues por imagen\n'
                 '(e = 0.05, orden_k = 2)',
                 fontsize=14, fontweight='bold', color=C['azul_osc'])

    y = np.arange(len(nombres))
    h = 0.38

    bars_a = ax.barh(y + h / 2, np.log10(antes), h,
                     color=C['gris_clar'], label='Antes',
                     edgecolor='white', linewidth=1.2)
    bars_d = ax.barh(y - h / 2, np.log10(despues), h,
                     color=C['azul'], label='Despues',
                     edgecolor='white', linewidth=1.2)

    for i, p in enumerate(perc):
        color_p = C['verde_osc'] if p >= 95 else C['verde']
        ax.text(np.log10(antes[i]) + 0.07, y[i] + h / 2,
                f'{p}%', va='center', fontsize=11,
                fontweight='bold', color=color_p)

    ax.set_yticks(y)
    ax.set_yticklabels(nombres, fontsize=11)
    ax.set_xlabel('Numero de vertices (escala log10)', fontsize=11)
    ax.set_xticks([3, 4, 5, 6])
    ax.set_xticklabels(['10^3', '10^4', '10^5', '10^6'], fontsize=11)
    ax.legend(fontsize=11, loc='lower right')
    ax.grid(True, axis='x', linestyle='--', alpha=0.35)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.91])
    return guardar('07_resultados.png', fig)


# ========================================================================
# GIF - Algoritmo completo animado
# ========================================================================
def generar_gif(paths_png):
    """
    Construye un GIF animado a partir de los 7 PNGs.
    Cada frame dura 2 segundos. Loop infinito.
    """
    pasos = [
        'Paso 1: Criterio de planitud\nicuando sobra un vertice?',
        'Paso 2: Face_circulator\nDCEL twin->next recorre la estrella de p',
        'Paso 3: Min-heap (Dijkstra)\nsiempre procesar el mas plano primero',
        'Paso 4: T.remove(p)\nCGAL elimina estrella y re-triangula',
        'Paso 5: Lazy deletion\ndetectar y descartar entradas obsoletas',
        'Paso 6: Vecindad de orden k\nsegundo anillo para actualizacion mas amplia',
        'Paso 7: Resultados\nsimplificacion 85-99% segun el terreno',
    ]

    frames = []
    for i, (path, paso) in enumerate(zip(paths_png, pasos)):
        img = Image.open(path).convert('RGBA')
        w, h = img.size

        # Agregar banda inferior con el numero de paso
        banda_h = 60
        frame = Image.new('RGBA', (w, h + banda_h), (255, 255, 255, 255))
        frame.paste(img, (0, 0))

        # Dibujar texto de paso en matplotlib y pegar
        fig_t, ax_t = plt.subplots(figsize=(w / 100, banda_h / 100))
        fig_t.patch.set_facecolor('#1F4E79')
        ax_t.set_facecolor('#1F4E79')
        ax_t.text(0.5, 0.5, f'({i + 1}/7)  {paso}',
                  ha='center', va='center', color='white',
                  fontsize=max(8, int(w / 120)),
                  fontweight='bold',
                  transform=ax_t.transAxes)
        ax_t.axis('off')
        fig_t.tight_layout(pad=0)

        # Guardar banda como imagen temporal
        tmp = os.path.join(OUT, f'_tmp_banda_{i}.png')
        fig_t.savefig(tmp, dpi=100, bbox_inches='tight',
                      facecolor='#1F4E79')
        plt.close(fig_t)

        banda = Image.open(tmp).convert('RGBA').resize((w, banda_h))
        frame.paste(banda, (0, h))
        os.remove(tmp)

        frames.append(frame.convert('RGB'))

    # Agregar frame final de resumen
    fig_r, ax_r = plt.subplots(figsize=(frames[0].size[0] / 100,
                                        frames[0].size[1] / 100))
    fig_r.patch.set_facecolor('white')
    ax_r.set_facecolor('white')
    ax_r.text(0.5, 0.65,
              'Algoritmo de Simplificacion de Heightmap\npor Vecindad de Orden k',
              ha='center', va='center', fontsize=18, fontweight='bold',
              color='#1F4E79', transform=ax_r.transAxes)
    ax_r.text(0.5, 0.42,
              'Taller 3 - Geometria Computacional\nPontificia Universidad Javeriana',
              ha='center', va='center', fontsize=13,
              color='#888780', transform=ax_r.transAxes)
    ax_r.text(0.5, 0.22,
              'DCEL - Delaunay - Dijkstra - Lazy deletion',
              ha='center', va='center', fontsize=12, style='italic',
              color='#2E75B6', transform=ax_r.transAxes)
    ax_r.axis('off')
    tmp_r = os.path.join(OUT, '_tmp_resumen.png')
    fig_r.savefig(tmp_r, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig_r)
    frames.append(Image.open(tmp_r).convert('RGB').resize(frames[0].size))
    os.remove(tmp_r)

    gif_path = os.path.join(OUT, 'algoritmo_completo.gif')
    frames[0].save(
        gif_path,
        save_all=True,
        append_images=frames[1:],
        duration=2200,
        loop=0,
        optimize=False,
    )
    print(f'  OK: algoritmo_completo.gif  ({len(frames)} frames)')
    return gif_path


# ========================================================================
# MAIN
# ========================================================================
if __name__ == '__main__':
    print('\n=== Generando diagramas del Taller 3 ===\n')
    paths = []
    paths.append(d1_criterio())
    paths.append(d2_circulator())
    paths.append(d3_heap())
    paths.append(d4_remove())
    paths.append(d5_lazy())
    paths.append(d6_vecindad())
    paths.append(d7_resultados())

    print('\n=== Generando GIF animado ===\n')
    generar_gif(paths)

    print(f'\n\N{CHECK MARK} Listo. Archivos en: ./{OUT}/')
    print('\nArchivos generados:')
    for f in sorted(os.listdir(OUT)):
        size = os.path.getsize(os.path.join(OUT, f))
        print(f'  {f:<45}  {size//1024} KB')