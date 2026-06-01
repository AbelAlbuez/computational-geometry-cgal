"""
Generate dashboard_lineal.html: interactive dashboard with metrics
(Dice, IoU, Hausdorff, areaErr) and a t-slider per pair.

For each pair (A, B, GT), runs the C++ binary `interpolation_lineal` 11 times
(t = 0.0, 0.1, ..., 1.0) and embeds all contour vertices in a single HTML
file. The slider redraws the green contour in real time and the t=0.5 result
is compared against GT to produce the metric row.

Usage:
    python scripts/generate_dashboard.py

Must be run from src/final-project-interpolation/ with the venv activated.
"""
import base64
import io
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from scipy.spatial.distance import directed_hausdorff
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs


ROOT     = Path(__file__).resolve().parent.parent
BINARY   = ROOT / "build" / "interpolation_lineal"
CONTOURS = ROOT / "data" / "contours"
OUT_HTML = ROOT / "src" / "interpolation-lineal" / "output" / "dashboard_lineal.html"

T_STEPS  = [round(i / 10.0, 1) for i in range(11)]      # 0.0 .. 1.0
T_INDEX_HALF = 5                                         # t=0.5 → metrics

PAIRS_SPEC = [
    {
        "case":     "BraTS-GLI-00008-100",
        "slice_a":  "slice_0070",
        "slice_b":  "slice_0072",
        "slice_gt": "slice_0071",
        "label":    "caso1_gt",
        "descripcion": "Caso 00008-100 — slices 70/72 con GT 71",
    },
    {
        "case":     "BraTS-GLI-00008-101",
        "slice_a":  "slice_0070",
        "slice_b":  "slice_0072",
        "slice_gt": "slice_0071",
        "label":    "caso2_gt",
        "descripcion": "Caso 00008-101 — slices 70/72 con GT 71",
    },
    {
        "case":     "BraTS-GLI-00009-101",
        "slice_a":  "slice_0070",
        "slice_b":  "slice_0072",
        "slice_gt": "slice_0071",
        "label":    "caso3_gt",
        "descripcion": "Caso 00009-101 — slices 70/72 con GT 71",
    },
    {
        "case":     "BraTS-GLI-00020-100",
        "slice_a":  "slice_0070",
        "slice_b":  "slice_0072",
        "slice_gt": "slice_0071",
        "label":    "caso4_gt",
        "descripcion": "Caso 00020-100 — slices 70/72 con GT 71",
    },
    {
        "case":     "BraTS-GLI-00528-101",
        "slice_a":  "slice_0070",
        "slice_b":  "slice_0072",
        "slice_gt": "slice_0071",
        "label":    "caso5_gt",
        "descripcion": "Caso 00528-101 — slices 70/72 con GT 71",
    },
]


# -----------------------------------------------------------------------------
# I/O helpers.
# -----------------------------------------------------------------------------

def read_obj_xy(path: Path):
    verts = []
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

    # Fallback: pick any 3 consecutive available slices.
    slices = sorted(base.glob("slice_*.obj"))
    for i in range(len(slices) - 2):
        s_a, s_gt, s_b = slices[i], slices[i + 1], slices[i + 2]
        try:
            na  = int(s_a.stem.split("_")[1])
            ngt = int(s_gt.stem.split("_")[1])
            nb  = int(s_b.stem.split("_")[1])
        except (IndexError, ValueError):
            continue
        if ngt == na + 1 and nb == na + 2:
            return {
                **spec,
                "slice_a":  s_a.stem,
                "slice_b":  s_b.stem,
                "slice_gt": s_gt.stem,
                "obj_a":    s_a,
                "obj_b":    s_b,
                "obj_gt":   s_gt,
            }
    raise FileNotFoundError(f"No consecutive triple found for {spec['case']}")


# -----------------------------------------------------------------------------
# Geometry / metrics.
# -----------------------------------------------------------------------------

def polygon_area(contour):
    if len(contour) < 3:
        return 0.0
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
    """Return a boolean mask of size n*n with True inside the polygon."""
    xs = np.linspace(xmin, xmax, n)
    ys = np.linspace(ymin, ymax, n)
    X, Y = np.meshgrid(xs, ys)
    pts = np.column_stack([X.ravel(), Y.ravel()])
    poly = MplPath(np.asarray(contour))
    mask = poly.contains_points(pts).reshape(n, n)
    return mask


def dice_iou(contour_a, contour_b):
    if len(contour_a) < 3 or len(contour_b) < 3:
        return 0.0, 0.0
    xmin, xmax, ymin, ymax = _bbox([contour_a, contour_b])
    ma = _rasterise(contour_a, xmin, xmax, ymin, ymax)
    mb = _rasterise(contour_b, xmin, xmax, ymin, ymax)
    inter = np.logical_and(ma, mb).sum()
    union = np.logical_or(ma, mb).sum()
    sum_  = ma.sum() + mb.sum()
    dice  = 0.0 if sum_  == 0 else 2.0 * inter / sum_
    iou   = 0.0 if union == 0 else inter / union
    return float(dice), float(iou)


def hausdorff(contour_a, contour_b):
    a = np.asarray(contour_a)
    b = np.asarray(contour_b)
    if a.size == 0 or b.size == 0:
        return float("nan")
    return float(max(
        directed_hausdorff(a, b)[0],
        directed_hausdorff(b, a)[0],
    ))


def area_err(contour, contour_gt):
    ag = polygon_area(contour_gt)
    if ag == 0:
        return float("nan")
    return abs(polygon_area(contour) - ag) / ag


# -----------------------------------------------------------------------------
# Binary invocation.
# -----------------------------------------------------------------------------

def run_binary(obj_a: Path, obj_b: Path, t: float, out_obj: Path):
    proc = subprocess.run(
        [str(BINARY), str(obj_a), str(obj_b), str(out_obj), f"{t:.3f}"],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    self_int = False
    for line in proc.stdout.splitlines():
        if line.startswith("Self-intersections detected:"):
            self_int = line.split(":")[1].strip().lower() == "yes"
    if proc.returncode != 0:
        raise RuntimeError(f"binary failed (t={t}): {proc.stderr}")
    return self_int


# -----------------------------------------------------------------------------
# Per-pair processing.
# -----------------------------------------------------------------------------

def process_pair(spec):
    print(f"=== {spec['label']}  ({spec['descripcion']}) ===")
    a  = read_obj_xy(spec["obj_a"])
    b  = read_obj_xy(spec["obj_b"])
    gt = read_obj_xy(spec["obj_gt"])
    print(f"  A:  {len(a)} v   B:  {len(b)} v   GT: {len(gt)} v")

    interpolated = {}
    self_int_05  = False
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for t in T_STEPS:
            out_obj = tmp_dir / f"t_{int(round(t * 100)):03d}.obj"
            si = run_binary(spec["obj_a"], spec["obj_b"], t, out_obj)
            interpolated[t] = read_obj_xy(out_obj)
            if abs(t - 0.5) < 1e-9:
                self_int_05 = si

    interp_05 = interpolated[0.5]
    dice, iou = dice_iou(interp_05, gt)
    hd        = hausdorff(interp_05, gt)
    aerr      = area_err(interp_05, gt)

    print(f"  Dice={dice:.4f}  IoU={iou:.4f}  "
          f"Hausdorff={hd:.3f}  areaErr={aerr:.4f}  "
          f"self-int={'yes' if self_int_05 else 'no'}")

    return {
        "label":       spec["label"],
        "descripcion": spec["descripcion"],
        "case":        spec["case"],
        "slice_a":     spec["slice_a"],
        "slice_b":     spec["slice_b"],
        "slice_gt":    spec["slice_gt"],
        "contour_a":   a,
        "contour_b":   b,
        "contour_gt":  gt,
        "interpolated": {f"{t:.1f}": interpolated[t] for t in T_STEPS},
        "metrics": {
            "dice":     dice,
            "iou":      iou,
            "hausdorff": hd,
            "area_err": aerr,
            "self_int": self_int_05,
        },
    }


# -----------------------------------------------------------------------------
# HTML rendering.
# -----------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<title>BraTS interpolation lineal — dashboard</title>
<style>
  body { font-family: system-ui, Arial, sans-serif; margin: 18px; color: #222; background: #fafafa; }
  h1   { font-size: 22px; margin-bottom: 4px; }
  h2   { font-size: 16px; margin-top: 28px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }
  .note { color: #555; font-size: 13px; max-width: 980px; }
  table { border-collapse: collapse; margin: 12px 0 24px; background: #fff; }
  th, td { border: 1px solid #ccc; padding: 6px 12px; text-align: right; font-size: 13px; }
  th { background: #f0f0f0; }
  td:first-child, th:first-child, td.txt, th.txt { text-align: left; }
  #metrics-table tbody tr:hover { background: #f3f8fe; }
  .good { background: #d4edda; }
  .mid  { background: #fff3cd; }
  .bad  { background: #f8d7da; }
  .panel { background: #fff; border: 1px solid #ddd; padding: 12px; border-radius: 6px; }
  .panel h3 { margin: 0 0 6px; font-size: 14px; }
  .summary { background: #fff; border: 1px solid #ddd; padding: 12px 16px; border-radius: 6px; max-width: 980px; }
  .summary-card { transition: border-color .12s ease; }
  .summary-card:hover { border-color: #378ADD !important; }
  .metric-card { background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 8px 10px; }
  .mc-label { font-size: 11px; color: #777; text-transform: uppercase; letter-spacing: .04em; }
  .mc-val   { font-size: 20px; font-weight: 700; color: #222; margin-top: 2px; }
</style>
__PLOTLY_BUNDLE__
</head>
<body>

<h1>BraTS glioma — interpolación lineal de contornos</h1>
<p class="note">
  Dashboard con métricas por par (Dice, IoU, Hausdorff, areaErr) comparando
  el contorno interpolado a <em>t</em>=0.5 contra el slice GT real intermedio.
  Haz clic en cualquier <em>card</em> de resumen o en una fila de la tabla para
  activar ese par en el visualizador; el <em>slider</em> recorre <em>t</em> ∈
  [0, 1] en pasos de 0.1 mostrando 11 fotogramas precalculados.
</p>

<h2>1. Resumen — 5 pares</h2>
<div id="cards-row" style="display:grid; grid-template-columns:repeat(5,1fr); gap:8px;"></div>
<p class="note">Card activa: borde azul. Click → activa el par en el visualizador (sección 3).</p>

<h2>2. Métricas por par</h2>
<table id="metrics-table">
  <thead>
    <tr>
      <th class="txt">Caso</th>
      <th class="txt">Slices A / GT / B</th>
      <th>Dice</th>
      <th>IoU</th>
      <th>Hausdorff</th>
      <th>areaErr</th>
      <th class="txt">Self-int</th>
    </tr>
  </thead>
  <tbody></tbody>
</table>
<p class="note">
  Verde: Dice ≥ 0.90 / IoU ≥ 0.80 / Hausdorff ≤ 3 / areaErr ≤ 0.05 ·
  Amarillo: zona intermedia · Rojo: peor desempeño. Fila activa resaltada en azul.
</p>

<h2>3. Visualizador de contornos</h2>
<div style="display:flex; gap:12px; align-items:flex-start; max-width:980px;">

  <div style="flex:1;">
    <img id="contour-img" style="width:100%; border-radius:8px; border:1px solid #ddd;" />
    <div style="display:flex; align-items:center; gap:10px; margin-top:8px;">
      <span style="font-size:13px; color:#555;">t =</span>
      <input id="t-slider" type="range" min="0" max="10" value="5"
             style="flex:1;" oninput="onSlider(this.value)">
      <span id="t-label" style="font-family:monospace; min-width:32px;">0.5</span>
    </div>
  </div>

  <div style="width:160px; display:flex; flex-direction:column; gap:8px;">
    <div class="metric-card"><div class="mc-label">Dice</div>
      <div id="mc-dice" class="mc-val"></div></div>
    <div class="metric-card"><div class="mc-label">IoU</div>
      <div id="mc-iou"  class="mc-val"></div></div>
    <div class="metric-card"><div class="mc-label">Hausdorff</div>
      <div id="mc-haus" class="mc-val"></div></div>
    <div class="metric-card"><div class="mc-label">Self-int</div>
      <div id="mc-si"   class="mc-val"></div></div>
  </div>

</div>

__STRESS_SECTION__

<h2>5. Visualización 3D — reconstrucción por mallado entre capas</h2>
<p class="note">
  Reconstrucción 3D (Plotly <code>mesh3d</code>) de los 11 contornos
  interpolados (<em>t</em>=0.0 a 1.0) apilados en el eje
  <em>z</em>=<em>t</em>×10. La superficie se obtiene triangulando entre
  capas consecutivas usando la correspondencia vértice a vértice ya
  garantizada por el <em>resampling</em> por longitud de arco; cada par
  de capas contiguas aporta 2<em>n</em> triángulos. La intensidad de
  color codifica <em>z</em> (escala RdYlGn). Rotar: clic izquierdo ·
  zoom: scroll · pan: clic derecho.
</p>
__FIG3D_HTML__

<h2>6. Resumen global</h2>
<div style="display:grid; grid-template-columns:repeat(4,1fr); gap:8px; max-width:980px;"
     id="summary-cards"></div>
<p class="note">
  Promedios sobre los 5 pares. GT = slice real intermedio (BraTS). El
  interpolado a t=0.5 se rasteriza en una grilla 512×512 dentro del bounding box
  común para calcular Dice/IoU; Hausdorff se computa directamente entre vértices.
</p>

<script>
var PAIRS = __PAIRS_JSON__;
var activeIdx = 0;

function cls(v, good, mid, lowerBetter) {
  if (typeof v !== "number" || isNaN(v)) return "";
  if (lowerBetter) {
    if (v <= good) return "good";
    if (v <= mid)  return "mid";
    return "bad";
  }
  if (v >= good) return "good";
  if (v >= mid)  return "mid";
  return "bad";
}

function fmt(v, d) { return (typeof v === "number" && !isNaN(v)) ? v.toFixed(d) : "–"; }

// -----------------------------------------------------------------------------
// Section 1: summary cards (one per pair).
// -----------------------------------------------------------------------------
var cardsRow = document.getElementById("cards-row");
PAIRS.forEach(function (p, i) {
  var m = p.metrics;
  var diceColor = m.dice >= 0.90 ? "#378ADD" : (m.dice >= 0.75 ? "#E8990C" : "#D64545");
  var siColor   = m.self_int ? "#E8990C" : "#2E9E5B";
  var siText    = m.self_int ? "self-int detectada" : "sin self-int";
  var card = document.createElement("div");
  card.className = "summary-card";
  card.style.cssText = "background:#fff;border:0.5px solid #ddd;border-radius:8px;padding:8px;cursor:pointer;";
  card.onclick = function () { activatePair(i); };
  card.innerHTML =
    '<img src="data:image/png;base64,' + p.frames["0.5"] + '" style="width:100%;border-radius:4px;display:block;">' +
    '<div style="font-size:11px;color:#777;margin-top:6px;">' + p.label + '</div>' +
    '<div style="font-size:24px;font-weight:700;line-height:1.1;color:' + diceColor + ';">' + m.dice.toFixed(3) + '</div>' +
    '<div style="font-size:10px;color:#999;text-transform:uppercase;letter-spacing:.04em;">Dice</div>' +
    '<div style="margin-top:6px;font-size:10px;font-weight:600;color:#fff;background:' + siColor +
      ';display:inline-block;padding:2px 6px;border-radius:10px;">' + siText + '</div>';
  cardsRow.appendChild(card);
});

// -----------------------------------------------------------------------------
// Section 2: clickable metrics table.
// -----------------------------------------------------------------------------
var tbody = document.querySelector("#metrics-table tbody");
PAIRS.forEach(function (p, i) {
  var m = p.metrics;
  var tr = document.createElement("tr");
  tr.style.cursor = "pointer";
  tr.onclick = function () { activatePair(i); };
  tr.innerHTML =
    '<td class="txt">' + p.label + '<br><span style="color:#777;font-size:11px">' + p.case + '</span></td>' +
    '<td class="txt">' + p.slice_a + ' / <b>' + p.slice_gt + '</b> / ' + p.slice_b + '</td>' +
    '<td class="' + cls(m.dice, 0.90, 0.75, false) + '">' + fmt(m.dice, 4) + '</td>' +
    '<td class="' + cls(m.iou, 0.80, 0.60, false) + '">' + fmt(m.iou, 4) + '</td>' +
    '<td class="' + cls(m.hausdorff, 3.0, 6.0, true) + '">' + fmt(m.hausdorff, 3) + '</td>' +
    '<td class="' + cls(m.area_err, 0.05, 0.15, true) + '">' + fmt(m.area_err, 4) + '</td>' +
    '<td class="txt">' + (m.self_int ? "yes" : "no") + '</td>';
  tbody.appendChild(tr);
});

// -----------------------------------------------------------------------------
// Section 3: slider + per-pair metric panel.
// -----------------------------------------------------------------------------
function onSlider(val) {
  var tStr = (parseInt(val, 10) / 10).toFixed(1);
  document.getElementById("t-label").textContent = tStr;
  var p = PAIRS[activeIdx];
  document.getElementById("contour-img").src =
    "data:image/png;base64," + p.frames[tStr];
}

function activatePair(i) {
  activeIdx = i;
  var p = PAIRS[i];
  var m = p.metrics;

  document.getElementById("t-slider").value = 5;
  onSlider(5);

  document.getElementById("mc-dice").textContent = fmt(m.dice, 3);
  document.getElementById("mc-iou").textContent  = fmt(m.iou, 3);
  document.getElementById("mc-haus").textContent = fmt(m.hausdorff, 2);
  document.getElementById("mc-si").textContent   = m.self_int ? "detectadas" : "no";

  document.querySelectorAll(".summary-card").forEach(function (c, j) {
    c.style.border = j === i ? "2px solid #378ADD" : "0.5px solid #ddd";
  });
  document.querySelectorAll("#metrics-table tbody tr").forEach(function (r, j) {
    r.style.background = j === i ? "#E6F1FB" : "";
  });
}

// -----------------------------------------------------------------------------
// Section 6: global summary cards (averages).
// -----------------------------------------------------------------------------
var mean = function (arr) { return arr.reduce(function (a, b) { return a + b; }, 0) / arr.length; };
var sc = document.getElementById("summary-cards");
[
  ["Dice promedio",      fmt(mean(PAIRS.map(function (p) { return p.metrics.dice;      })), 4)],
  ["IoU promedio",       fmt(mean(PAIRS.map(function (p) { return p.metrics.iou;       })), 4)],
  ["Hausdorff promedio", fmt(mean(PAIRS.map(function (p) { return p.metrics.hausdorff; })), 3)],
  ["areaErr promedio",   fmt(mean(PAIRS.map(function (p) { return p.metrics.area_err;  })), 4)],
].forEach(function (kv) {
  var d = document.createElement("div");
  d.className = "metric-card";
  d.innerHTML = '<div class="mc-label">' + kv[0] + '</div><div class="mc-val">' + kv[1] + '</div>';
  sc.appendChild(d);
});

// -----------------------------------------------------------------------------
// Init.
// -----------------------------------------------------------------------------
activatePair(0);
</script>

</body>
</html>
"""


STRESS_ROWS = [
    (10, 0, "no"),
    (5,  4, "no"),
    (3,  4, "no"),
    (1,  4, "no"),
]


def build_stress_section() -> str:
    stress_dir = ROOT / "src" / "interpolation-lineal" / "output" / "stress_test"
    img_block  = ""
    if stress_dir.exists():
        runs = sorted(p for p in stress_dir.iterdir() if p.is_dir())
        if runs:
            png = runs[-1] / "stress_comparison.png"
            if png.exists():
                b64_data = base64.b64encode(png.read_bytes()).decode()
                src = f"data:image/png;base64,{b64_data}"
                img_block = (
                    '<div class="panel" style="margin-top:14px;max-width:980px;">'
                    '<h3>Estrella (r=1) vs cuadrado — salida del binario</h3>'
                    f'<img src="{src}" alt="Comparativa stress test" '
                    'style="width:100%;height:auto;display:block;'
                    'border:1px solid #ddd;border-radius:4px;">'
                    '<p class="note" style="margin-top:6px;">'
                    'Izquierda: entradas A (azul) y B (rojo). Derecha: '
                    'interpolación producida por el binario (verde). Pese a la '
                    'geometría adversarial, el resultado es un polígono simple '
                    'sin cruces.</p></div>'
                )

    rows = "\n".join(
        f'    <tr><td>{r}</td><td>{c}</td>'
        f'<td class="txt good">{s}</td></tr>'
        for (r, c, s) in STRESS_ROWS
    )

    return f"""<h2>4. Experimento de stress — Robustez del pipeline</h2>
<p class="note">
  Los contornos ET reales de BraTS son empíricamente demasiado convexos
  para producir auto-intersecciones bajo interpolación lineal con
  <em>best_rotation</em>. Para forzar un escenario adversarial diseñé un
  experimento sintético: contorno <em>A</em> = estrella de 4 puntas
  (8 vértices, R=50) vs. contorno <em>B</em> = cuadrado de 8 vértices
  (esquinas + puntos medios, lado=60), ambos centrados en el mismo punto.
  Barrí el radio interior <em>r</em> de la estrella sobre {{10, 5, 3, 1}},
  y como control independiente verifiqué los cruces de la interpolación
  cruda vértice a vértice con una rutina O(n²) en Python.
</p>
<table id="stress-table">
  <thead>
    <tr>
      <th>r</th>
      <th>Cruces (Python, crudo)</th>
      <th class="txt">Self-int (binario)</th>
    </tr>
  </thead>
  <tbody>
{rows}
  </tbody>
</table>
<div class="summary" style="max-width:980px;">
  <b>Conclusión:</b> la pre-alineación (CCW + centroide + best_rotation)
  elimina los cruces antes de que el resolver los vea. El
  <code>SelfIntersectionResolver</code> actúa como salvaguarda de última
  línea para casos que el pipeline no anticipa.
</div>
{img_block}
"""


def _plot_poly(ax, contour, **kw):
    """Plot a closed polygon (repeating the first vertex to close the loop)."""
    if not contour:
        return
    xs = [p[0] for p in contour] + [contour[0][0]]
    ys = [p[1] for p in contour] + [contour[0][1]]
    ax.plot(xs, ys, **kw)


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def build_contour_images(results):
    """Render 11 matplotlib PNG frames (t = 0.0, 0.1, … 1.0) per pair and embed
    each as base64. A (blue dotted), B (red dotted) and GT (orange dashed) are
    fixed; only the green interpolated contour changes per frame. Axis limits
    are shared across frames so the polygon does not jump while sliding t.
    Uses ax.set_aspect("equal") + ax.invert_yaxis() (image coordinates).

    Returns a list of dicts:
        {label, descripcion, case, slice_a/b/gt, metrics, frames={t: b64}}.
    """
    out = []
    for r in results:
        a, b, gt    = r["contour_a"], r["contour_b"], r["contour_gt"]
        interp_lyr  = r["interpolated"]

        allpts = a + b + gt + [p for c in interp_lyr.values() for p in c]
        xs = [p[0] for p in allpts]
        ys = [p[1] for p in allpts]
        pad  = 0.05 * max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
        xlim = (min(xs) - pad, max(xs) + pad)
        ylim = (min(ys) - pad, max(ys) + pad)

        frames = {}
        for i in range(11):
            t_str  = f"{i / 10:.1f}"
            interp = interp_lyr.get(t_str, [])

            fig, ax = plt.subplots(figsize=(4.6, 4.6))
            _plot_poly(ax, a,  color="#1f77b4", ls=":",  lw=1.5, label="A")
            _plot_poly(ax, b,  color="#d62728", ls=":",  lw=1.5, label="B")
            _plot_poly(ax, gt, color="#ff7f0e", ls="--", lw=1.5, label="GT")
            _plot_poly(ax, interp, color="#2ca02c", ls="-", lw=2.6,
                       label=f"interp t={t_str}")
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_title(f"{r['label']} — t = {t_str}", fontsize=11)
            ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
            ax.tick_params(labelsize=8)
            frames[t_str] = _fig_to_b64(fig)

        out.append({
            "label":       r["label"],
            "descripcion": r["descripcion"],
            "case":        r["case"],
            "slice_a":     r["slice_a"],
            "slice_b":     r["slice_b"],
            "slice_gt":    r["slice_gt"],
            "metrics":     r["metrics"],
            "frames":      frames,
        })
    return out


def build_3d_figure(results):
    """Build a single plotly.graph_objects.Figure with one Mesh3d trace per
    pair. Each trace stacks the 11 interpolated contours at z=t*10 and
    triangulates between consecutive layers (correspondence guaranteed by
    the arc-length resampling already done by the C++ binary)."""
    fig = go.Figure()
    for r in results:
        label  = r["label"]
        layers = r["interpolated"]
        ns = [len(pts) for pts in layers.values()]
        if len(set(ns)) != 1 or ns[0] == 0:
            print(f"[SKIP] {label}: capas con n distinto {set(ns)}",
                  file=sys.stderr)
            continue
        n = ns[0]

        pts0 = layers["0.0"]
        cx = sum(p[0] for p in pts0) / n
        cy = sum(p[1] for p in pts0) / n

        sorted_layers = sorted(layers.items(), key=lambda kv: float(kv[0]))
        xs, ys, zs = [], [], []
        for t_str, pts in sorted_layers:
            z = float(t_str) * 10.0
            for p in pts:
                xs.append(p[0] - cx)
                ys.append(p[1] - cy)
                zs.append(z)

        num_layers = len(sorted_layers)
        ii, jj, kk = [], [], []
        for layer in range(num_layers - 1):
            base = layer * n
            for v in range(n):
                j = (v + 1) % n
                ii.append(base + v);     jj.append(base + j);         kk.append(base + n + v)
                ii.append(base + j);     jj.append(base + n + j);     kk.append(base + n + v)

        fig.add_trace(go.Mesh3d(
            x=xs, y=ys, z=zs,
            i=ii, j=jj, k=kk,
            intensity=zs,
            colorscale="RdYlGn",
            reversescale=True,
            opacity=0.7,
            name=label,
            showscale=True,
            colorbar=dict(title="z = t×10"),
            flatshading=False,
        ))

    fig.update_layout(
        title="Reconstrucción 3D — contornos interpolados",
        paper_bgcolor="#1a1a2e",
        plot_bgcolor="#1a1a2e",
        font=dict(color="#ddd"),
        scene=dict(
            xaxis=dict(title="x (mm)", backgroundcolor="#1a1a2e",
                       gridcolor="#444466", zerolinecolor="#666688"),
            yaxis=dict(title="y (mm)", backgroundcolor="#1a1a2e",
                       gridcolor="#444466", zerolinecolor="#666688"),
            zaxis=dict(title="z = t×10", backgroundcolor="#1a1a2e",
                       gridcolor="#444466", zerolinecolor="#666688"),
            bgcolor="#1a1a2e",
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, t=40, b=0),
        height=520,
    )
    return fig


def main():
    if not BINARY.exists():
        print(f"ERROR: missing binary {BINARY}. Build first.", file=sys.stderr)
        sys.exit(1)

    resolved = [resolve_pair(spec) for spec in PAIRS_SPEC]
    results  = [process_pair(spec) for spec in resolved]

    # Section 3 viewer: 11 base64 PNG frames per pair (no Plotly needed —
    # the slider just swaps <img> sources). This payload also drives the
    # section-1 cards, the section-2 table and the section-6 averages.
    print("Rendering contour frames…")
    contour_images = build_contour_images(results)
    payload = json.dumps(contour_images, ensure_ascii=False, separators=(",", ":"))
    html    = HTML_TEMPLATE.replace("__PAIRS_JSON__", payload)

    # Only the 3D figure (section 5) needs Plotly. Embed the full bundle once,
    # inline in <head>, exactly like the working report.html export — fully
    # self-contained (no CDN, no local paths) and verified to render.
    html = html.replace(
        "__PLOTLY_BUNDLE__",
        f'<script type="text/javascript">{get_plotlyjs()}</script>',
    )

    fig3d     = build_3d_figure(results)
    fig3d.update_layout(height=820)
    fig3d_html = fig3d.to_html(
        full_html=False,
        include_plotlyjs=False,
        div_id="plotly-3d-main",
        config={"responsive": True, "displaylogo": False},
    )
    html = html.replace("__FIG3D_HTML__", fig3d_html)

    html = html.replace("__STRESS_SECTION__", build_stress_section())
    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(html, encoding="utf-8")
    print(f"\nDashboard escrito en {OUT_HTML.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
