#!/usr/bin/env python3
# =============================================================================
# Genera un INFORME HTML autocontenido (en español) desde la carpeta de
# resultados: resumen ejecutivo, métodos y matemáticas (MathJax), gráficas de
# barras y de líneas (interactivas), distribución (box plots), selector 3-D de
# tumores (verdad de terreno vs lineal/spline/sdf), tablas y un análisis
# post-mortem. plotly.js va embebido una sola vez; MathJax desde CDN.
#
#   python3 build_report.py --results data/results --out report.html
# =============================================================================
from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics as st
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs

COLS = ["gt", "linear", "spline", "sdf"]
COL_ES = {"gt": "verdad de terreno", "linear": "lineal", "spline": "spline", "sdf": "sdf"}
COL_COLOR = {"gt": "#8fce8f", "linear": "#9ec5fe", "spline": "#f5c2a0", "sdf": "#c5a3e0"}
M_COLOR = {"linear": "#1f77b4", "spline": "#ff7f0e", "sdf": "#9467bd"}
METHODS = ["linear", "spline", "sdf"]
M_ES = {"linear": "lineal", "spline": "spline", "sdf": "sdf"}


# ---------- carga de datos ---------------------------------------------------
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


def pct(v, q):
    s = sorted(v)
    if not s:
        return 0.0
    i = min(len(s) - 1, max(0, int(q * len(s))))
    return s[i]


# ---------- figuras ----------------------------------------------------------
def _resp(fig):
    return fig.to_html(full_html=False, include_plotlyjs=False,
                       config={"responsive": True, "displaylogo": False})


def fig_bars(agg):
    fig = make_subplots(rows=1, cols=3, subplot_titles=(
        "Solape (Dice / IoU, mayor = mejor)",
        "Hausdorff medio (mm, menor = mejor)", "% de cortes con Dice > 0.8"))
    m = [M_ES[a["method"]] for a in agg]
    fig.add_trace(go.Bar(name="Dice media", x=m,
                  y=[float(a["dice_mean"]) for a in agg], marker_color="#1f77b4"), 1, 1)
    fig.add_trace(go.Bar(name="Dice mediana", x=m,
                  y=[float(a["dice_median"]) for a in agg], marker_color="#7fb3e0"), 1, 1)
    fig.add_trace(go.Bar(name="IoU media", x=m,
                  y=[float(a["iou_mean"]) for a in agg], marker_color="#aec7e8"), 1, 1)
    fig.add_trace(go.Bar(showlegend=False, x=m,
                  y=[float(a["hausdorff_mean"]) for a in agg], marker_color="#ff7f0e"), 1, 2)
    fig.add_trace(go.Bar(showlegend=False, x=m,
                  y=[float(a["dice_pct_gt80"]) for a in agg], marker_color="#2ca02c"), 1, 3)
    fig.update_yaxes(range=[0, 1], row=1, col=1)
    fig.update_yaxes(range=[0, 100], row=1, col=3)
    fig.update_layout(barmode="group", height=380, margin=dict(l=10, r=10, t=46, b=10),
                      legend=dict(orientation="h", y=1.2, x=0.0),
                      paper_bgcolor="white", plot_bgcolor="#fafbfd")
    return _resp(fig)


def fig_cdf(raw):
    fig = go.Figure()
    for m in METHODS:
        d = np.sort([float(r["dice"]) for r in raw if r["method"] == m])
        if len(d) == 0:
            continue
        y = np.arange(1, len(d) + 1) / len(d)
        fig.add_trace(go.Scatter(x=d, y=y, mode="lines", name=M_ES[m],
                                 line=dict(color=M_COLOR[m], width=2.5)))
    fig.add_vline(x=0.8, line_dash="dash", line_color="gray",
                  annotation_text="Dice = 0.8")
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=20, b=10),
                      xaxis_title="Dice", yaxis_title="fracción de cortes ≤ Dice",
                      xaxis_range=[0, 1], yaxis_range=[0, 1],
                      legend=dict(x=0.02, y=0.98), plot_bgcolor="#fafbfd")
    return _resp(fig)


def fig_percase(percase):
    fig = go.Figure()
    base = sorted(percase, key=lambda c: -float(percase[c].get("linear", {}).get("dice", 0)))
    x = list(range(1, len(base) + 1))
    for m in METHODS:
        y = [float(percase[c].get(m, {}).get("dice", "nan")) for c in base]
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=M_ES[m],
                      line=dict(color=M_COLOR[m], width=1.6)))
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=20, b=10),
                      xaxis_title="caso (ordenado por Dice lineal)",
                      yaxis_title="Dice medio por caso", yaxis_range=[0, 1],
                      legend=dict(x=0.02, y=0.05), plot_bgcolor="#fafbfd")
    return _resp(fig)


def fig_box(raw, key, title, color_lower=True):
    fig = go.Figure()
    for m in METHODS:
        v = [float(r[key]) for r in raw if r["method"] == m]
        fig.add_trace(go.Box(y=v, name=M_ES[m], marker_color=M_COLOR[m],
                             boxmean=True, boxpoints=False))
    fig.update_layout(height=340, margin=dict(l=10, r=10, t=36, b=10),
                      title=title, showlegend=False, plot_bgcolor="#fafbfd")
    return _resp(fig)


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
    header = ["método", "Dice", "IoU", "Hausdorff (mm)", "error área",
              "vol malla (mm³)", "vol VT (mm³)", "razón"]
    cols = [[M_ES[m] for m in METHODS],
            [f(m, "dice") for m in METHODS], [f(m, "iou") for m in METHODS],
            [f(m, "hausdorff", 2) for m in METHODS], [f(m, "rel_area_err") for m in METHODS],
            [vol(m) or "—" for m in METHODS], [gt or "—"] * 3, [ratio(m) for m in METHODS]]
    return header, cols


def fig_selector(res, percase, recon):
    cases = sorted(os.path.basename(p)[5:-7] for p in glob.glob(str(res / "mesh_*_gt.off")))
    if not cases:
        return "", 0
    fig = make_subplots(rows=2, cols=4,
        specs=[[{"type": "scene"}] * 4, [{"type": "table", "colspan": 4}, None, None, None]],
        row_heights=[0.74, 0.26], vertical_spacing=0.06, horizontal_spacing=0.005,
        subplot_titles=[COL_ES[m] for m in COLS] + [""])
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
        fig.add_trace(go.Table(header=dict(values=header, fill_color="#dfe8f3",
                      align="left", font=dict(size=12)),
                      cells=dict(values=cols, align="left", font=dict(size=12), height=24),
                      visible=(ci == 0)), row=2, col=1)
        trace_case.append(ci)
    ntr = len(trace_case)
    buttons = [dict(label=c, method="update",
               args=[{"visible": [trace_case[t] == ci for t in range(ntr)]},
                     {"title.text": f"Tumor {c}  —  verdad de terreno vs lineal / spline / sdf"}])
               for ci, c in enumerate(cases)]
    fig.update_scenes(aspectmode="data", xaxis_visible=False, yaxis_visible=False,
                      zaxis_visible=False)
    fig.update_layout(height=820, margin=dict(l=0, r=0, t=110, b=0),
                      title=dict(text=f"Tumor {cases[0]}  —  verdad de terreno vs lineal / spline / sdf",
                                 x=0.5, xanchor="center", y=0.99),
                      updatemenus=[dict(buttons=buttons, direction="down", showactive=True,
                                        x=0.5, xanchor="center", y=1.10, yanchor="top",
                                        bgcolor="#eef2f8")],
                      annotations=list(fig.layout.annotations) + [dict(
                          text="Seleccione un tumor:", x=0.5, xanchor="right", y=1.135,
                          yref="paper", xref="paper", showarrow=False, font=dict(size=13),
                          xshift=-130)])
    return _resp(fig), len(cases)


# ---------- bloques HTML -----------------------------------------------------
def tbl(rows, headers, cls=""):
    th = "".join(f"<th>{h}</th>" for h in headers)
    trs = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table class='{cls}'><thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table>"


MATE = r"""
<div class="grid2">
  <div class="card method linear">
    <h3>M1 — Lineal (correspondencia de vértices)</h3>
    <p>Se remuestrean $A,B$ a $N$ puntos equiespaciados en longitud de arco, se
    busca el desfase cíclico que mejor los alinea,
    $$\sigma^\star=\arg\min_{\sigma}\sum_{i=0}^{N-1}\bigl\lVert A_i-B_{(i+\sigma)\bmod N}\bigr\rVert^2,$$
    y se interpola vértice a vértice con $t\in[0,1]$:
    $$C^{(i)}_{t}=(1-t)\,A_i+t\,B_{(i+\sigma^\star)\bmod N}.$$
    Las auto-intersecciones residuales se eliminan con un barrido de
    Bentley–Ottmann.</p>
  </div>
  <div class="card method spline">
    <h3>M2 — Polinomial / spline cúbico natural (en z)</h3>
    <p>Encadenando la correspondencia entre cortes, cada vértice traza una
    trayectoria $\gamma_i(z)$. El <b>spline cúbico natural</b> es $C^2$ y sus
    segundas derivadas $M_k$ resuelven el sistema tridiagonal
    $$h_{k-1}M_{k-1}+2(h_{k-1}+h_k)M_k+h_kM_{k+1}
     =6\!\left(\frac{y_{k+1}-y_k}{h_k}-\frac{y_k-y_{k-1}}{h_{k-1}}\right),$$
    con $M_0=M_m=0$. La variante <b>polinomial</b> es un ajuste por mínimos
    cuadrados de grado $d=\min(3,m-1)$ (contraste didáctico: sub-ajusta).</p>
  </div>
  <div class="card method sdf">
    <h3>M3 — Campo de distancia con signo (SDF)</h3>
    <p>Cada contorno se vuelve un campo de distancia con signo,
    $\Phi_C(\mathbf q)=\pm\min_{\mathbf c\in\partial C}\lVert\mathbf q-\mathbf c\rVert$
    (positivo dentro; distancia con un árbol AABB de CGAL, signo con
    <code>Polygon_2::bounded_side</code>). La partición por arista más cercana es
    el <b>diagrama de Voronoi</b> de las aristas y su cresta es el <b>eje medio</b>.
    El nuevo contorno es el conjunto de nivel cero:
    $$\Phi_t=(1-t)\,\Phi_A+t\,\Phi_B,\qquad C_t=\{\mathbf q:\Phi_t(\mathbf q)=0\}.$$</p>
  </div>
  <div class="card method recon">
    <h3>Reconstrucción 3-D y métricas</h3>
    <p>La pila densificada se convierte en una nube orientada $(x,y,k\,\Delta z)$ y
    se cierra con reconstrucción de superficie de <b>Poisson</b> (que construye una
    triangulación de Delaunay 3-D). El volumen usa el teorema de la divergencia,
    $V=\tfrac16\sum_{\triangle}\mathbf a\cdot(\mathbf b\times\mathbf c)$. La exactitud
    (dejar-un-corte-fuera) usa <b>Dice</b> $=\tfrac{2\lvert X\cap Y\rvert}{\lvert X\rvert+\lvert Y\rvert}$,
    <b>IoU</b>, la <b>Hausdorff</b> simétrica y el error relativo de área.</p>
  </div>
</div>"""


def postmortem():
    items = [
        ("Entorno CGAL 5.6 (no 6.x)",
         "El proyecto se construye en WSL con CGAL 5.6, no 6.x. Hubo que usar "
         "<code>CGAL::AABB_traits</code> (no <code>AABB_traits_3</code>) y, como 5.6 no "
         "incluye <code>interpolated_corrected_curvatures</code>, la curvatura media se "
         "implementó a mano (Laplaciano cotangente)."),
        ("Eigen ausente, sin sudo",
         "Poisson y <code>jet_estimate_normals</code> requieren Eigen, que no estaba "
         "instalado y sin contraseña de sudo. Solución: Eigen es solo cabeceras → se "
         "descargó el .deb como usuario y se apuntó CMake con <code>-DEIGEN3_INCLUDE_DIR</code>."),
        ("La librería compartida no compilaba",
         "<code>lib/pujCGAL/SegmentsIntersection.hxx</code> usaba <code>std::get_if</code> "
         "sobre lo que en CGAL 5.6 es un <code>boost::variant</code>. Se hizo portable: "
         "<code>boost::get</code> en 5.x y <code>std::get_if</code> en ≥ 6.0 según "
         "<code>CGAL_VERSION_NR</code>."),
        ("Convención de etiquetas BraTS",
         "El tumor activo es la etiqueta 3, pero los primeros casos están "
         "post-tratamiento (solo etiquetas 2 y 4) y casi no tienen ET, por lo que la "
         "extracción no producía contornos. El estudio selecciona los casos con más "
         "cortes de ET reales."),
        ("Colapso de la reconstrucción",
         "En pilas de anillos irregulares, <code>mst_orient_normals</code> descartaba la "
         "mayoría de los puntos y la malla colapsaba (un caso: 387 vs 29015 mm³). Solución: "
         "no descartar; usar como respaldo la normal exterior en el plano del contorno "
         "(fiable en la superficie lateral). Tras el arreglo: 32957 mm³."),
        ("El spline sobre-oscila en datos reales",
         "En sintético el spline ganaba; en datos reales es el peor. Los cortes reales "
         "cambian de forma de manera irregular en z, y el ajuste suave del spline "
         "sobre-dispara. Lineal y SDF (dos vecinos) son más robustos."),
        ("Tumores muy pequeños aún fallan",
         "Tumores de pocos miles de mm³ (p. ej. 02093-100/101) tienen contornos demasiado "
         "escasos para una superficie de Poisson estable y se reconstruyen incompletos. "
         "Se dejan en el selector como modo de fallo honesto."),
        ("Huecos en z en el dejar-uno-fuera",
         "Las regiones de ET tienen huecos (cortes con poco tumor). La validación corrige "
         "$t$ con las alturas reales de los cortes y descarta tríos con hueco &gt; 3."),
        ("Escalabilidad",
         "La validación lanza el binario una vez por corte: los 1350 casos serían horas. "
         "Se ejecutó en segundo plano (110 casos); el dataset completo requeriría una "
         "validación en proceso (sin subproceso por corte)."),
        ("Cuotas WSL/PowerShell e integración",
         "El paso de comandos PowerShell→WSL→bash rompía el entrecomillado; se resolvió "
         "escribiendo scripts y ejecutándolos por ruta con espacios escapados. Además se "
         "integró la interpolación lineal del repositorio compartido "
         "(<code>LinearInterpolator</code> + resolución de auto-intersecciones) con el resto "
         "del pipeline."),
    ]
    cards = "".join(
        f"<div class='card pm'><h4>{i+1}. {t}</h4><p>{d}</p></div>"
        for i, (t, d) in enumerate(items))
    return f"<div class='grid2'>{cards}</div>"


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
    aggm = {a["method"]: a for a in agg}

    # figuras
    bar = fig_bars(agg)
    cdf = fig_cdf(raw)
    line = fig_percase(percase)
    box_d = fig_box(raw, "dice", "Distribución de Dice por método")
    box_h = fig_box(raw, "hausdorff", "Distribución de Hausdorff (mm) por método")
    sel, n_vis = fig_selector(res, percase, recon)

    # estadísticos detallados
    stat_rows = []
    for m in METHODS:
        d = [float(r["dice"]) for r in raw if r["method"] == m]
        h = [float(r["hausdorff"]) for r in raw if r["method"] == m]
        a = [float(r["rel_area_err"]) for r in raw if r["method"] == m]
        stat_rows.append([M_ES[m], len(d), f"{st.mean(d):.3f}", f"{st.median(d):.3f}",
                          f"{st.pstdev(d):.3f}", f"{pct(d,0.1):.3f}", f"{pct(d,0.9):.3f}",
                          f"{st.mean(h):.2f}", f"{st.median(a):.3f}"])
    stat_tbl = tbl(stat_rows, ["método", "n", "Dice media", "Dice mediana", "Dice σ",
                               "Dice p10", "Dice p90", "Hausd. media (mm)", "err.área med."],
                   "data")

    # victorias por caso (qué método logra el mayor Dice en cada caso)
    wins = {m: 0 for m in METHODS}
    for c, d in percase.items():
        best, bv = None, -1
        for m in METHODS:
            try:
                v = float(d.get(m, {}).get("dice", -1))
            except Exception:
                v = -1
            if v > bv:
                bv, best = v, m
        if best:
            wins[best] += 1
    win_rows = [[M_ES[m], wins[m], f"{100*wins[m]/max(1,ncases):.0f}%"] for m in METHODS]
    win_tbl = tbl(win_rows, ["método", "casos ganados", "porcentaje"], "data")

    # tabla reconstrucción
    rec_rows = []
    for (case, m), r in sorted(recon.items()):
        mv = r["mesh_volume_mm3"].replace(" mm^3", "").strip()
        gt = r["voxel_gt_volume_mm3"].strip()
        try:
            ratio = f"{float(mv)/float(gt):.2f}"
        except Exception:
            ratio = "—"
        rec_rows.append([case, M_ES.get(m, m), r.get("status", ""), mv or "—", gt, ratio])
    rec_tbl = tbl(rec_rows, ["caso", "método", "estado", "vol malla (mm³)",
                             "vol VT (mm³)", "razón"], "data")

    def g(m, k):
        return float(aggm[m][k])
    style = """
<style>
:root{--azul:#1f77b4;--naranja:#ff7f0e;--morado:#9467bd;--verde:#2ca02c;
  --tinta:#1f2a37;--suave:#5b6b7c;--linea:#e3e8ef;--fondo:#f6f8fb}
*{box-sizing:border-box}
body{font-family:'Segoe UI',system-ui,Arial,sans-serif;margin:0;color:var(--tinta);
  background:var(--fondo);line-height:1.55}
.wrap{max-width:1200px;margin:0 auto;padding:0 22px 60px}
header.hero{background:linear-gradient(120deg,#1e3c72,#2a5298 60%,#3a7bd5);
  color:#fff;padding:38px 22px 30px;margin-bottom:8px}
header.hero .wrap{padding-bottom:0}
header.hero h1{margin:0 0 6px;font-size:27px;font-weight:700}
header.hero p{margin:0;color:#dce6f7;max-width:820px}
.kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:20px 0 6px}
.kpi{background:#fff;border-radius:12px;padding:14px 16px;box-shadow:0 1px 3px rgba(20,40,80,.08);
  border-top:4px solid var(--azul)}
.kpi .v{font-size:24px;font-weight:700} .kpi .l{color:var(--suave);font-size:13px}
.kpi.l2{border-color:var(--morado)} .kpi.l3{border-color:var(--verde)} .kpi.l4{border-color:var(--naranja)}
h2{margin:38px 0 6px;font-size:21px;border-left:5px solid var(--azul);padding-left:12px}
h3{margin:18px 0 6px;font-size:16px} h4{margin:0 0 6px;font-size:15px}
p{max-width:980px}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:14px 0}
.grid2.plots{align-items:start}
.card{background:#fff;border-radius:12px;padding:16px 18px;box-shadow:0 1px 3px rgba(20,40,80,.07);
  border:1px solid var(--linea)}
.card.method{border-top:5px solid var(--azul)}
.card.method.spline{border-top-color:var(--naranja)}
.card.method.sdf{border-top-color:var(--morado)}
.card.method.recon{border-top-color:var(--verde)}
.card.pm{border-left:4px solid var(--naranja);background:#fff8f1}
.callout{background:#eef5ff;border:1px solid #cfe0fb;border-left:5px solid var(--azul);
  border-radius:10px;padding:14px 18px;margin:14px 0}
.callout b{color:#15396b}
.fig{background:#fff;border:1px solid var(--linea);border-radius:12px;padding:8px;margin:12px 0}
table{border-collapse:collapse;margin:10px 0 6px;width:100%;background:#fff;border-radius:10px;
  overflow:hidden;box-shadow:0 1px 3px rgba(20,40,80,.06)}
th,td{border-bottom:1px solid var(--linea);padding:7px 12px;text-align:right;font-size:13.5px}
th{background:#eef2f8;color:#2a3f5f} td:first-child,th:first-child{text-align:left}
tbody tr:nth-child(even){background:#fafbfd}
.note{color:var(--suave);font-size:13px}
.pill{display:inline-block;padding:1px 9px;border-radius:999px;color:#fff;font-size:12px;font-weight:600}
.pill.lin{background:var(--azul)} .pill.spl{background:var(--naranja)} .pill.sdf{background:var(--morado)}
code{background:#eef1f6;padding:1px 5px;border-radius:5px;font-size:12.5px}
@media(max-width:860px){.grid2,.kpis{grid-template-columns:1fr}}
</style>"""
    mathjax = ("<script>MathJax={tex:{inlineMath:[['$','$'],['\\\\(','\\\\)']],"
               "displayMath:[['$$','$$'],['\\\\[','\\\\]']]}};</script>"
               "<script async src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>")
    head = (f"<!doctype html><html lang='es'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>Informe — interpolación de contornos tumorales</title>"
            f"{style}{mathjax}<script>{get_plotlyjs()}</script></head><body>")

    hero = (f"<header class='hero'><div class='wrap'>"
            f"<h1>Interpolación de contornos tumorales, reconstrucción 3-D y comparación</h1>"
            f"<p>Evaluación con datos reales de <b>BraTS 2024 GLI</b> (tumor realzante). "
            f"Validación <i>dejar-un-corte-fuera</i> sobre <b>{ncases} casos</b> "
            f"(<b>{nslices} cortes</b> ocultados por método) y <b>{n_vis}</b> tumores "
            f"reconstruidos en 3-D.</p></div></header>")

    kpis = (f"<div class='kpis'>"
            f"<div class='kpi'><div class='v'>{g('linear','dice_median'):.3f}</div>"
            f"<div class='l'>Dice mediana — lineal</div></div>"
            f"<div class='kpi l2'><div class='v'>{g('sdf','dice_median'):.3f}</div>"
            f"<div class='l'>Dice mediana — sdf</div></div>"
            f"<div class='kpi l3'><div class='v'>{g('linear','dice_pct_gt80'):.0f}%</div>"
            f"<div class='l'>cortes con Dice &gt; 0.8 (lineal)</div></div>"
            f"<div class='kpi l4'><div class='v'>{ncases}</div>"
            f"<div class='l'>casos evaluados</div></div></div>")

    resumen = (
        "<h2>1. Resumen ejecutivo</h2>"
        + kpis +
        "<div class='callout'><b>Conclusión principal.</b> La interpolación es muy buena en "
        "el corte típico (Dice mediana ≈ 0.92): cerca del 80 % de los cortes ocultados se "
        f"reproducen con Dice &gt; 0.8. El método <span class='pill lin'>lineal</span> y el "
        f"<span class='pill sdf'>sdf</span> lideran; el <span class='pill spl'>spline</span> "
        "queda por detrás en datos reales. La <b>media</b> es menor que la <b>mediana</b> "
        "porque una minoría de cortes difíciles (tumor diminuto o que se divide) puntúa casi "
        "0 y arrastra el promedio.</div>")

    metodos = ("<h2>2. Métodos y matemáticas</h2>"
               "<p>Un contorno en el corte axial $z$ es una poligonal cerrada "
               "$C=(p_0,\\dots,p_{n-1})$. Se normaliza a sentido antihorario mediante el "
               "área con signo (fórmula del cordón) "
               "$A(C)=\\tfrac12\\sum_i (x_i y_{i+1}-x_{i+1} y_i)$.</p>" + MATE)

    barras = ("<h2>3. Métricas — gráficas de barras</h2><div class='fig'>" + bar + "</div>"
              "<p class='note'>Izquierda: solape de región (Dice/IoU, mayor es mejor). "
              "Centro: Hausdorff simétrica media (mm, menor es mejor). Derecha: fracción de "
              "cortes con Dice &gt; 0.8.</p>")

    lineas = ("<h2>4. Métricas — gráficas de líneas y distribución</h2>"
              "<div class='grid2 plots'>"
              "<div><h3>Distribución de Dice (CDF empírica)</h3><div class='fig'>" + cdf +
              "</div></div>"
              "<div><h3>Dice medio por caso</h3><div class='fig'>" + line + "</div></div>"
              "</div>"
              "<div class='grid2 plots'><div class='fig'>" + box_d + "</div>"
              "<div class='fig'>" + box_h + "</div></div>")

    tablas = ("<h2>5. Tablas de métricas</h2>"
              "<h3>Estadísticos por método (110 casos)</h3>" + stat_tbl +
              "<h3>¿Qué método gana en cada caso? (mayor Dice por caso)</h3>" + win_tbl)

    tres_d = ("<h2>6. Reconstrucciones 3-D (seleccione un tumor)</h2>"
              "<p class='note'>Cada panel es interactivo: arrastre para rotar, rueda para "
              "acercar. La verdad de terreno es la superficie por <i>marching cubes</i> de los "
              "vóxeles; las demás son reconstrucciones de Poisson de la pila interpolada.</p>"
              "<div class='fig'>" + sel + "</div>"
              "<h3>Volumen reconstruido vs. verdad de terreno (vóxeles)</h3>" + rec_tbl)

    disc = (
        "<h2>7. Discusión detallada</h2>"
        "<div class='grid2'>"
        "<div class='card'><h4>Lineal y SDF por encima del spline</h4>"
        f"<p>Lineal alcanza Dice media {g('linear','dice_mean'):.3f} (mediana "
        f"{g('linear','dice_median'):.3f}) y SDF {g('sdf','dice_mean'):.3f} "
        f"({g('sdf','dice_median'):.3f}); el spline se queda en {g('spline','dice_mean'):.3f} "
        f"({g('spline','dice_median'):.3f}). En datos reales la forma cambia de manera "
        "irregular entre cortes y el ajuste suave del spline sobre-dispara; los métodos de "
        "dos vecinos son más estables. Es lo <i>contrario</i> al caso sintético suave, donde "
        "el spline ganaba — buena señal de que la métrica es sensible a la irregularidad real.</p></div>"
        "<div class='card'><h4>Mediana ≫ media</h4>"
        f"<p>La diferencia entre mediana ({g('linear','dice_median'):.2f}) y media "
        f"({g('linear','dice_mean'):.2f}) indica una distribución sesgada: el corte típico es "
        "excelente, pero una cola de cortes difíciles (extremos del tumor, regiones que se "
        "dividen o se funden) hunde la media. Por eso reportamos también <b>mediana</b> y "
        "<b>% &gt; 0.8</b> como cifras más representativas.</p></div>"
        "<div class='card'><h4>Volumen vs. verdad de terreno</h4>"
        "<p>Para tumores sólidos la razón volumen-malla / volumen-vóxel ≈ 1.0–1.1. Cuando es "
        "mayor, Poisson envuelve el contorno exterior y rellena núcleos necróticos/anulares "
        "que la máscara de vóxeles deja huecos. Es una sobre-estimación interpretable, no un "
        "error de cálculo.</p></div>"
        "<div class='card'><h4>Hausdorff y error de área</h4>"
        f"<p>La Hausdorff media (lineal {g('linear','hausdorff_mean'):.1f} mm) está dominada por "
        "los mismos cortes atípicos; la mediana es mucho menor. El error de área mediano "
        f"(~{g('sdf','area_err_median')*100:.0f}–{g('spline','area_err_median')*100:.0f} %) "
        "confirma que el tamaño del corte interpolado se conserva bien salvo en los casos "
        "difíciles.</p></div>"
        "</div>")

    pm = ("<h2>8. Análisis post-mortem — problemas y soluciones</h2>"
          "<p>El camino hasta estos resultados tuvo varios obstáculos técnicos. Se documentan "
          "para reproducibilidad y como aprendizaje.</p>" + postmortem())

    concl = (
        "<h2>9. Conclusiones</h2>"
        "<div class='callout'><p>Se construyó un pipeline completo en C++/CGAL: tres métodos "
        "de interpolación de contornos, reconstrucción 3-D por Poisson y un conjunto de "
        "métricas, validado con datos reales de BraTS. La interpolación <b>lineal</b> es la "
        "opción más fiable en datos reales, con <b>SDF</b> muy cerca y ventaja potencial en "
        "cambios de topología; el <b>spline</b> es útil pero sobre-dispara. Las "
        "reconstrucciones son cerradas y su volumen coincide con la verdad de terreno salvo en "
        "tumores muy pequeños. Próximos pasos: validación en proceso para todo el dataset y "
        "manejo de multifocalidad (varios contornos por corte).</p></div>")

    repro = ("<h2>10. Reproducibilidad</h2>"
             "<p class='note'>Dentro de WSL (ver <code>RUN_LOCAL.md</code>):</p>"
             "<div class='card'><code>python3 scripts/run_real_study.py --dataset "
             "/ruta/BraTS/training_data1_v2 --bin ~/builds/contour/contour_interpolator "
             "--extract 250 --topk 110 --n-vis 10</code><br><br>"
             "<code>python3 scripts/build_report.py --results data/results --out report.html</code></div>")

    body = ("<div class='wrap'>" + resumen + metodos + barras + lineas + tablas +
            tres_d + disc + pm + concl + repro + "</div></body></html>")

    Path(args.out).write_text(head + hero + body, encoding="utf-8")
    print(f"escrito {args.out}  ({os.path.getsize(args.out)/1e6:.1f} MB, "
          f"{ncases} casos, {n_vis} en 3D)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
