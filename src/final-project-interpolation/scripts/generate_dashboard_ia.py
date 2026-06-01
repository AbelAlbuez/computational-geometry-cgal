"""
Generate dashboard_ia.html: interactive dashboard with metrics
(Dice, IoU, Hausdorff, areaErr) y malla 3D generada por la IA (VFI).
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MplPath
from scipy.spatial.distance import directed_hausdorff
import plotly.graph_objects as go

ROOT       = Path(__file__).resolve().parent.parent
BINARY     = ROOT / "build" / "interpolation_lineal"
BINARY_IA  = ROOT / "build" / "interpolador_ia"
CONTOURS   = ROOT / "data" / "contours"
OUT_HTML   = ROOT / "src" / "interpolation-lineal" / "output" / "dashboard_ia.html"

T_STEPS  = [round(i / 10.0, 1) for i in range(11)]
T_INDEX_HALF = 5

# Casos de prueba reales
PAIRS_SPEC = [
    {"case": "BraTS-GLI-00008-100", "slice_a": "slice_0070", "slice_gt": "slice_0071", "slice_b": "slice_0072", "label": "Caso_00008_100", "descripcion": "Slices 70 y 72 (GT: 71)"},
    {"case": "BraTS-GLI-00008-101", "slice_a": "slice_0070", "slice_gt": "slice_0071", "slice_b": "slice_0072", "label": "Caso_00008_101", "descripcion": "Slices 70 y 72 (GT: 71)"},
    {"case": "BraTS-GLI-00008-103", "slice_a": "slice_0052", "slice_gt": "slice_0053", "slice_b": "slice_0054", "label": "Caso_00008_103", "descripcion": "Slices 52 y 54 (GT: 53)"},
    {"case": "BraTS-GLI-00009-100", "slice_a": "slice_0088", "slice_gt": "slice_0089", "slice_b": "slice_0090", "label": "Caso_00009_100", "descripcion": "Slices 88 y 90 (GT: 89)"},
    {"case": "BraTS-GLI-00009-101", "slice_a": "slice_0085", "slice_gt": "slice_0086", "slice_b": "slice_0087", "label": "Caso_00009_101", "descripcion": "Slices 85 y 87 (GT: 86)"},
    {"case": "BraTS-GLI-00020-100", "slice_a": "slice_0074", "slice_gt": "slice_0075", "slice_b": "slice_0076", "label": "Caso_00020_100", "descripcion": "Slices 74 y 76 (GT: 75)"},
    {"case": "BraTS-GLI-00020-101", "slice_a": "slice_0073", "slice_gt": "slice_0074", "slice_b": "slice_0075", "label": "Caso_00020_101", "descripcion": "Slices 73 y 75 (GT: 74)"},
    {"case": "BraTS-GLI-00027-100", "slice_a": "slice_0095", "slice_gt": "slice_0096", "slice_b": "slice_0097", "label": "Caso_00027_100", "descripcion": "Slices 95 y 97 (GT: 96)"},
    {"case": "BraTS-GLI-00027-101", "slice_a": "slice_0095", "slice_gt": "slice_0096", "slice_b": "slice_0097", "label": "Caso_00027_101", "descripcion": "Slices 95 y 97 (GT: 96)"},
    {"case": "BraTS-GLI-00046-100", "slice_a": "slice_0126", "slice_gt": "slice_0127", "slice_b": "slice_0128", "label": "Caso_00046_100", "descripcion": "Slices 126 y 128 (GT: 127)"},
]

# -----------------------------------------------------------------------------
# I/O helpers.
# -----------------------------------------------------------------------------

def read_obj_xy(path: Path):
    verts = []
    if not path.exists():
        return verts
    with open(path, "r") as fh:
        for line in fh:
            if not line or line[0] != "v":
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                verts.append((float(parts[1]), float(parts[2])))
            except ValueError:
                continue
    return verts

def resolve_pair(spec):
    base = CONTOURS / spec["case"]
    a  = base / f"{spec['slice_a']}.obj"
    b  = base / f"{spec['slice_b']}.obj"
    gt = base / f"{spec['slice_gt']}.obj"
    if a.exists() and b.exists() and gt.exists():
        return {**spec, "obj_a": a, "obj_b": b, "obj_gt": gt}
    raise FileNotFoundError(f"No triple found for {spec['case']}")

# -----------------------------------------------------------------------------
# Geometry / metrics.
# -----------------------------------------------------------------------------

def polygon_area(contour):
    if len(contour) < 3: return 0.0
    n = len(contour)
    s = 0.0
    for i in range(n):
        x1, y1 = contour[i]
        x2, y2 = contour[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2.0

def _bbox(contours, pad=2.0):
    pts = np.vstack([np.asarray(c) for c in contours if len(c)])
    return (pts[:, 0].min() - pad, pts[:, 0].max() + pad,
            pts[:, 1].min() - pad, pts[:, 1].max() + pad)

def _rasterise(contour, xmin, xmax, ymin, ymax, n=512):
    if len(contour) < 3: return np.zeros((n, n), dtype=bool)
    xs = np.linspace(xmin, xmax, n)
    ys = np.linspace(ymin, ymax, n)
    X, Y = np.meshgrid(xs, ys)
    pts = np.column_stack([X.ravel(), Y.ravel()])
    poly = MplPath(np.asarray(contour))
    return poly.contains_points(pts).reshape(n, n)

def dice_iou(contour_a, contour_b):
    if len(contour_a) < 3 or len(contour_b) < 3: return 0.0, 0.0
    xmin, xmax, ymin, ymax = _bbox([contour_a, contour_b])
    ma = _rasterise(contour_a, xmin, xmax, ymin, ymax)
    mb = _rasterise(contour_b, xmin, xmax, ymin, ymax)
    inter = np.logical_and(ma, mb).sum()
    union = np.logical_or(ma, mb).sum()
    sum_  = ma.sum() + mb.sum()
    return (0.0 if sum_ == 0 else 2.0 * inter / sum_), (0.0 if union == 0 else inter / union)

def hausdorff(contour_a, contour_b):
    a, b = np.asarray(contour_a), np.asarray(contour_b)
    if a.size == 0 or b.size == 0: return float("nan")
    return float(max(directed_hausdorff(a, b)[0], directed_hausdorff(b, a)[0]))

def area_err(contour, contour_gt):
    ag = polygon_area(contour_gt)
    return float("nan") if ag == 0 else abs(polygon_area(contour) - ag) / ag

# -----------------------------------------------------------------------------
# Binary invocations.
# -----------------------------------------------------------------------------

def run_binary(obj_a: Path, obj_b: Path, t: float, out_obj: Path):
    proc = subprocess.run([str(BINARY), str(obj_a), str(obj_b), str(out_obj), f"{t:.3f}"], cwd=str(ROOT), capture_output=True, text=True)
    self_int = False
    for line in proc.stdout.splitlines():
        if line.startswith("Self-intersections detected:"):
            self_int = line.split(":")[1].strip().lower() == "yes"
    return self_int

def run_ia_binary(obj_a: Path, obj_b: Path, out_obj: Path):
    proc = subprocess.run([str(BINARY_IA), str(obj_a), str(obj_b), str(out_obj)], cwd=str(ROOT), capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"[!] Advertencia: La IA falló para {obj_a.name}. Detalle: {proc.stderr}")
        return False
    return True

# -----------------------------------------------------------------------------
# Per-pair processing.
# -----------------------------------------------------------------------------

def process_pair(spec):
    print(f"=== {spec['label']}  ({spec['descripcion']}) ===")
    a  = read_obj_xy(spec["obj_a"])
    b  = read_obj_xy(spec["obj_b"])
    gt = read_obj_xy(spec["obj_gt"])
    
    interpolated = {}
    self_int_05  = False
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        
        # 1. C++ Lineal
        for t in T_STEPS:
            out_obj = tmp_dir / f"t_{int(round(t * 100)):03d}.obj"
            si = run_binary(spec["obj_a"], spec["obj_b"], t, out_obj)
            interpolated[t] = read_obj_xy(out_obj)
            if abs(t - 0.5) < 1e-9: self_int_05 = si

        # 2. Pipeline IA (VFI)
        ia_obj_path = tmp_dir / "ia_interp_mid.obj"
        ia_metrics = None
        ia_contour = []
        if run_ia_binary(spec["obj_a"], spec["obj_b"], ia_obj_path) and ia_obj_path.exists():
            ia_contour = read_obj_xy(ia_obj_path)
            if len(ia_contour) > 2:
                i_dice, i_iou = dice_iou(ia_contour, gt)
                i_hd = hausdorff(ia_contour, gt)
                i_aerr = area_err(ia_contour, gt)
                ia_metrics = {"dice": i_dice, "iou": i_iou, "hausdorff": i_hd, "area_err": i_aerr}
                print(f"  [IA-VFI] Dice={i_dice:.4f}  IoU={i_iou:.4f}  Hausdorff={i_hd:.3f}")

    interp_05 = interpolated[0.5]
    dice, iou = dice_iou(interp_05, gt)
    
    return {
        "label": spec["label"], "descripcion": spec["descripcion"],
        "case": spec["case"], "slice_a": spec["slice_a"], "slice_b": spec["slice_b"], "slice_gt": spec["slice_gt"],
        "contour_a": a, "contour_b": b, "contour_gt": gt, "contour_ia": ia_contour, 
        "interpolated": {f"{t:.1f}": interpolated[t] for t in T_STEPS},
        "metrics": {"dice": dice, "iou": iou, "hausdorff": hausdorff(interp_05, gt), "area_err": area_err(interp_05, gt), "self_int": self_int_05},
        "metrics_ia": ia_metrics
    }

# -----------------------------------------------------------------------------
# IA 3D Mesh Generation (Correspondencia Polar & Arc-Length en Python)
# -----------------------------------------------------------------------------

def align_and_resample(contour, n=150):
    """ Remuestreo por longitud de arco y alineación por ángulo polar (como en C++) """
    if not contour or len(contour) < 3:
        return [[0.0, 0.0]] * (n + 1)
    
    pts = np.array(contour)
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack((pts, pts[0]))
        
    dp = np.diff(pts, axis=0)
    l = np.cumsum(np.sqrt(np.sum(dp**2, axis=1)))
    l = np.insert(l, 0, 0)
    if l[-1] == 0: return [[pts[0][0], pts[0][1]]] * (n + 1)
        
    l_norm = l / l[-1]
    t_target = np.linspace(0, 1, n, endpoint=False)
    x = np.interp(t_target, l_norm, pts[:,0])
    y = np.interp(t_target, l_norm, pts[:,1])
    resampled = np.column_stack((x, y))
    
    # Alineación por ángulo polar (Centroide)
    cx, cy = np.mean(resampled, axis=0)
    angles = np.arctan2(resampled[:,1] - cy, resampled[:,0] - cx)
    start_idx = np.argmin(angles)
    
    aligned = np.roll(resampled, -start_idx, axis=0)
    aligned = np.vstack((aligned, aligned[0])) # Cerrar el loop para malla
    return aligned.tolist()

def build_3d_figure(results):
    fig = go.Figure()
    buttons = []
    
    for idx, r in enumerate(results):
        label = r["label"]
        ca = r.get("contour_a", [])
        cia = r.get("contour_ia", [])
        cb = r.get("contour_b", [])
        
        if not cia: continue
            
        # Generar las 3 capas remuestreadas uniformemente
        n_pts = 150
        layer_a = align_and_resample(ca, n_pts)
        layer_ia = align_and_resample(cia, n_pts)
        layer_b = align_and_resample(cb, n_pts)
        
        layers = [layer_a, layer_ia, layer_b]
        z_vals = [0.0, 5.0, 10.0]
        
        cx = sum(p[0] for p in layer_ia) / len(layer_ia)
        cy = sum(p[1] for p in layer_ia) / len(layer_ia)
        
        xs, ys, zs = [], [], []
        for z, pts in zip(z_vals, layers):
            for p in pts:
                xs.append(p[0] - cx)
                ys.append(p[1] - cy)
                zs.append(z)
                
        ii, jj, kk = [], [], []
        N = n_pts + 1
        for layer in range(2):
            base = layer * N
            for v in range(n_pts):
                nxt = v + 1
                ii.append(base + v);   jj.append(base + nxt); kk.append(base + N + v)
                ii.append(base + nxt); jj.append(base + N + nxt); kk.append(base + N + v)
                
        fig.add_trace(go.Mesh3d(
            x=xs, y=ys, z=zs, i=ii, j=jj, k=kk,
            intensity=zs, colorscale="Purples", opacity=0.85, name=label, showscale=True,
            colorbar=dict(title="Cortes (Z)"), flatshading=False,
            visible=(idx == 0) # Solo mostramos el primer caso por defecto
        ))
        
        # Opciones para el Dropdown
        vis_array = [False] * len(results)
        vis_array[idx] = True
        buttons.append(dict(
            label=label, method="update",
            args=[{"visible": vis_array}]
        ))

    fig.update_layout(
        paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e", font=dict(color="#ddd"),
        updatemenus=[dict(
            active=0, buttons=buttons,
            x=0.0, y=1.15, xanchor="left", yanchor="top",
            bgcolor="#333", font=dict(color="#fff")
        )],
        scene=dict(
            xaxis=dict(title="x", backgroundcolor="#1a1a2e", gridcolor="#444466", zerolinecolor="#666688"),
            yaxis=dict(title="y", backgroundcolor="#1a1a2e", gridcolor="#444466", zerolinecolor="#666688"),
            zaxis=dict(title="z", backgroundcolor="#1a1a2e", gridcolor="#444466", zerolinecolor="#666688", range=[-1, 11]),
            bgcolor="#1a1a2e", aspectmode="auto" # <--- ARREGLA LA ALTURA APLASTADA
        ),
        margin=dict(l=0, r=0, t=60, b=0), height=550,
    )
    return fig

# -----------------------------------------------------------------------------
# HTML rendering.
# -----------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<title>BraTS interpolation híbrida — dashboard</title>
<style>
  body { font-family: system-ui, Arial, sans-serif; margin: 18px; color: #222; background: #fafafa; }
  h1   { font-size: 22px; margin-bottom: 4px; }
  h2   { font-size: 16px; margin-top: 28px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }
  .note { color: #555; font-size: 13px; max-width: 980px; }
  table { border-collapse: collapse; margin: 12px 0 24px; background: #fff; width: 100%; max-width: 980px; }
  th, td { border: 1px solid #ccc; padding: 6px 12px; text-align: right; font-size: 13px; }
  th { background: #f0f0f0; }
  td:first-child, th:first-child, td.txt, th.txt { text-align: left; }
  .good { background: #d4edda; }
  .mid  { background: #fff3cd; }
  .bad  { background: #f8d7da; }
  .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }
  .panel { background: #fff; border: 1px solid #ddd; padding: 12px; border-radius: 6px; }
  .panel h3 { margin: 0 0 6px; font-size: 14px; }
  .slider-row { display: flex; align-items: center; gap: 10px; margin: 6px 0 10px; }
  .slider-row input[type=range] { flex: 1; }
  .t-label { font-family: ui-monospace, Menlo, monospace; font-size: 13px; min-width: 64px; }
  #viz-grid svg { width: 100%; height: 420px; background: #fff; }
  .legend { font-size: 12px; color: #444; margin-top: 4px; }
  .legend span { display: inline-block; margin-right: 12px; }
  .swatch { display: inline-block; width: 12px; height: 3px; vertical-align: middle; margin-right: 4px; }
  .summary { background: #fff; border: 1px solid #ddd; padding: 12px 16px; border-radius: 6px; max-width: 980px; }
</style>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
</head>
<body>

<h1>BraTS glioma — Comparativa C++ Geométrica vs Red Neuronal (VFI)</h1>
<p class="note">
  Dashboard interactivo con métricas por par comparando el enfoque paramétrico lineal en C++ 
  contra el modelo de flujo óptico neuronal (Baseline IA). 
</p>

<h2>1. Métricas por par</h2>
<table id="metrics-table">
  <thead>
    <tr>
      <th class="txt">Caso</th>
      <th class="txt">Slices A / GT / B</th>
      <th class="txt">Método</th>
      <th>Dice</th>
      <th>IoU</th>
      <th>Hausdorff</th>
      <th>areaErr</th>
      <th class="txt">Self-int</th>
    </tr>
  </thead>
  <tbody></tbody>
</table>

<h2>2. Comparativa por métrica</h2>
<div class="grid">
  <div class="panel"><h3>Dice (mayor es mejor)</h3><canvas id="chart-dice"></canvas></div>
  <div class="panel"><h3>IoU (mayor es mejor)</h3><canvas id="chart-iou"></canvas></div>
  <div class="panel"><h3>Hausdorff (menor es mejor)</h3><canvas id="chart-hd"></canvas></div>
  <div class="panel"><h3>areaErr (menor es mejor)</h3><canvas id="chart-aerr"></canvas></div>
</div>

<h2>3. Visualizador interactivo — slider de <em>t</em> y comparativa IA</h2>
<div id="viz-grid" class="grid"></div>

<h2>4. Visualización 3D — Reconstrucción IA (VFI)</h2>
<p class="note">
  Se muestra la reconstrucción de la malla tridimensional basada <strong>exclusivamente en la inferencia de la Inteligencia Artificial</strong>. 
  Se mapean los vértices mediante correspondencia de ángulo polar nativa desde el slice A (z=0) hacia la predicción de la red neuronal (z=5) y finalmente al slice B (z=10). Selecciona el caso en el menú desplegable.
</p>
__SECTION4_HTML__

<h2>5. Resumen General</h2>
<div class="summary" id="summary"></div>

<script>
const PAIRS = __PAIRS_JSON__;

function cls(v, good, mid, lowerBetter=false) {
  if (Number.isNaN(v) || v === null || v === undefined) return "";
  if (lowerBetter) {
    if (v <= good) return "good";
    if (v <= mid)  return "mid";
    return "bad";
  } else {
    if (v >= good) return "good";
    if (v >= mid)  return "mid";
    return "bad";
  }
}

const tbody = document.querySelector("#metrics-table tbody");
for (const p of PAIRS) {
  const m = p.metrics;
  const hasIA = !!p.metrics_ia;
  
  const tr = document.createElement("tr");
  tr.innerHTML = `
    <td class="txt" rowspan="${hasIA ? 2 : 1}">${p.label}<br><span style="color:#777;font-size:11px">${p.case}</span></td>
    <td class="txt" rowspan="${hasIA ? 2 : 1}">${p.slice_a} / <b>${p.slice_gt}</b> / ${p.slice_b}</td>
    <td class="txt"><b>Lineal (C++)</b></td>
    <td class="${cls(m.dice, 0.90, 0.75)}">${m.dice.toFixed(4)}</td>
    <td class="${cls(m.iou,  0.80, 0.60)}">${m.iou.toFixed(4)}</td>
    <td class="${cls(m.hausdorff, 3.0, 6.0, true)}">${m.hausdorff.toFixed(3)}</td>
    <td class="${cls(m.area_err, 0.05, 0.15, true)}">${m.area_err.toFixed(4)}</td>
    <td class="txt">${m.self_int ? "yes" : "no"}</td>
  `;
  tbody.appendChild(tr);

  if (hasIA) {
    const mia = p.metrics_ia;
    const tria = document.createElement("tr");
    tria.innerHTML = `
      <td class="txt" style="background:#f4ebff"><b>IA (VFI)</b></td>
      <td class="${cls(mia.dice, 0.90, 0.75)}"><b>${mia.dice.toFixed(4)}</b></td>
      <td class="${cls(mia.iou,  0.80, 0.60)}"><b>${mia.iou.toFixed(4)}</b></td>
      <td class="${cls(mia.hausdorff, 3.0, 6.0, true)}"><b>${mia.hausdorff.toFixed(3)}</b></td>
      <td class="${cls(mia.area_err, 0.05, 0.15, true)}"><b>${mia.area_err.toFixed(4)}</b></td>
      <td class="txt" style="background:#f4ebff">N/A</td>
    `;
    tbody.appendChild(tria);
  }
}

const labels = PAIRS.map(p => p.label);
function bar(id, key, title) {
  new Chart(document.getElementById(id), {
    type: "bar",
    data: {
      labels: labels,
      datasets: [
        {
          label: "Lineal (C++)",
          data:  PAIRS.map(p => p.metrics[key]),
          backgroundColor: "rgba(33,150,243,0.7)",
          borderColor:     "rgba(33,150,243,1)",
          borderWidth: 1,
        },
        {
          label: "IA (VFI)",
          data:  PAIRS.map(p => p.metrics_ia ? p.metrics_ia[key] : 0),
          backgroundColor: "rgba(156,39,176,0.7)",
          borderColor:     "rgba(156,39,176,1)",
          borderWidth: 1,
        }
      ],
    },
    options: {
      responsive: true,
      plugins: { legend: { display: true } },
      scales:  { y: { beginAtZero: true } },
    },
  });
}
bar("chart-dice", "dice",      "Dice");
bar("chart-iou",  "iou",       "IoU");
bar("chart-hd",   "hausdorff", "Hausdorff");
bar("chart-aerr", "area_err",  "areaErr");

function bboxOf(...contours) {
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (const c of contours) for (const [x, y] of c) {
    if (x < xmin) xmin = x; if (x > xmax) xmax = x;
    if (y < ymin) ymin = y; if (y > ymax) ymax = y;
  }
  const pad = 0.05 * Math.max(xmax - xmin, ymax - ymin, 1);
  return { xmin: xmin - pad, xmax: xmax + pad, ymin: ymin - pad, ymax: ymax + pad };
}

function pathD(contour, bb, w, h) {
  if (!contour || !contour.length) return "";
  const sx = x => ((x - bb.xmin) / (bb.xmax - bb.xmin)) * w;
  const sy = y => h - ((y - bb.ymin) / (bb.ymax - bb.ymin)) * h;
  let d = `M ${sx(contour[0][0])} ${sy(contour[0][1])}`;
  for (let i = 1; i < contour.length; ++i)
    d += ` L ${sx(contour[i][0])} ${sy(contour[i][1])}`;
  d += " Z";
  return d;
}

const vizGrid = document.getElementById("viz-grid");
for (const p of PAIRS) {
  const panel = document.createElement("div");
  panel.className = "panel";
  panel.innerHTML = `
    <h3>${p.label} — ${p.descripcion}</h3>
    <div class="slider-row">
      <span class="t-label">t = <b id="tval-${p.label}">0.5</b></span>
      <input type="range" min="0" max="10" value="5" id="slider-${p.label}">
    </div>
    <svg id="svg-${p.label}" viewBox="0 0 600 420" preserveAspectRatio="xMidYMid meet">
      <path id="pa-${p.label}"  fill="none" stroke="#1f77b4" stroke-width="1.4" stroke-dasharray="4 3"/>
      <path id="pb-${p.label}"  fill="none" stroke="#d62728" stroke-width="1.4" stroke-dasharray="4 3"/>
      <path id="pgt-${p.label}" fill="none" stroke="#ff7f0e" stroke-width="1.4" stroke-dasharray="4 3"/>
      <path id="pia-${p.label}" fill="none" stroke="#9c27b0" stroke-width="2.5" stroke-dasharray="6 4" opacity="0.8"/>
      <path id="pi-${p.label}"  fill="none" stroke="#2ca02c" stroke-width="2.2"/>
    </svg>
    <div class="legend">
      <span><span class="swatch" style="background:#1f77b4"></span>A</span>
      <span><span class="swatch" style="background:#d62728"></span>B</span>
      <span><span class="swatch" style="background:#ff7f0e"></span>GT</span>
      <span><span class="swatch" style="background:#2ca02c"></span>Int.(C++)</span>
      ${p.contour_ia && p.contour_ia.length ? `<span><span class="swatch" style="background:#9c27b0"></span>IA (t=0.5)</span>` : ""}
    </div>
  `;
  vizGrid.appendChild(panel);

  const allInterp = Object.values(p.interpolated).flat();
  const bb = bboxOf(p.contour_a, p.contour_b, p.contour_gt, allInterp, p.contour_ia);
  const W = 600, H = 420;
  document.getElementById(`pa-${p.label}`).setAttribute("d",  pathD(p.contour_a,  bb, W, H));
  document.getElementById(`pb-${p.label}`).setAttribute("d",  pathD(p.contour_b,  bb, W, H));
  document.getElementById(`pgt-${p.label}`).setAttribute("d", pathD(p.contour_gt, bb, W, H));
  
  if(p.contour_ia && p.contour_ia.length > 0) {
      document.getElementById(`pia-${p.label}`).setAttribute("d", pathD(p.contour_ia, bb, W, H));
  }

  const piEl   = document.getElementById(`pi-${p.label}`);
  const tLabel = document.getElementById(`tval-${p.label}`);
  const slider = document.getElementById(`slider-${p.label}`);

  function setT(idx) {
    const tStr = (idx / 10).toFixed(1);
    tLabel.textContent = tStr;
    const c = p.interpolated[tStr] || [];
    piEl.setAttribute("d", pathD(c, bb, W, H));
  }
  setT(5);
  slider.addEventListener("input", e => setT(parseInt(e.target.value, 10)));
}

const mean = arr => arr.length === 0 ? 0 : arr.reduce((a, b) => a + b, 0) / arr.length;
const summary = document.getElementById("summary");
const iaPairs = PAIRS.filter(p => p.metrics_ia);

let iaStats = "";
if(iaPairs.length > 0) {
    iaStats = `
    <b>Promedio IA (VFI) sobre ${iaPairs.length} pares evaluados:</b><br>
    Dice = ${mean(iaPairs.map(p => p.metrics_ia.dice)).toFixed(4)} ·
    IoU = ${mean(iaPairs.map(p => p.metrics_ia.iou)).toFixed(4)} ·
    Hausdorff = ${mean(iaPairs.map(p => p.metrics_ia.hausdorff)).toFixed(3)} ·
    areaErr = ${mean(iaPairs.map(p => p.metrics_ia.area_err)).toFixed(4)}
    <br><br>
    `;
}

summary.innerHTML = `
  <b>Promedio Lineal (C++) sobre ${PAIRS.length} pares:</b><br>
  Dice = ${mean(PAIRS.map(p => p.metrics.dice)).toFixed(4)} ·
  IoU = ${mean(PAIRS.map(p => p.metrics.iou)).toFixed(4)} ·
  Hausdorff = ${mean(PAIRS.map(p => p.metrics.hausdorff)).toFixed(3)} ·
  areaErr = ${mean(PAIRS.map(p => p.metrics.area_err)).toFixed(4)}
  <br><br>
  ${iaStats}
`;
</script>
</body>
</html>
"""

def main():
    if not BINARY.exists():
        print(f"ERROR: missing binary {BINARY}. Build first.", file=sys.stderr)
        sys.exit(1)
    if not BINARY_IA.exists():
        print(f"ERROR: missing IA binary {BINARY_IA}. Build first.", file=sys.stderr)
        sys.exit(1)

    resolved = [resolve_pair(spec) for spec in PAIRS_SPEC]
    results  = [process_pair(spec) for spec in resolved]

    payload = json.dumps(results, ensure_ascii=False, separators=(",", ":"))
    html    = HTML_TEMPLATE.replace("__PAIRS_JSON__", payload)

    fig3d         = build_3d_figure(results)
    section4_html = fig3d.to_html(
        full_html=False,
        include_plotlyjs=True,
        div_id="plotly-3d-main",
        config={"responsive": True, "displaylogo": False},
    )
    section4_html = (
        '<div style="width:100%;min-height:520px;display:block;">'
        + section4_html
        + '</div>'
    )
    html = html.replace("__SECTION4_HTML__", section4_html)

    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(html, encoding="utf-8")
    print(f"\nDashboard IA híbrido escrito en {OUT_HTML.relative_to(ROOT)}")

if __name__ == "__main__":
    main()