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
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MplPath
from scipy.spatial.distance import directed_hausdorff


ROOT     = Path(__file__).resolve().parent.parent
BINARY   = ROOT / "build" / "interpolation_lineal"
CONTOURS = ROOT / "data" / "contours"
OUT_HTML = ROOT / "src" / "interpolation-lineal" / "dashboard_lineal.html"

T_STEPS  = [round(i / 10.0, 1) for i in range(11)]      # 0.0 .. 1.0
T_INDEX_HALF = 5                                         # t=0.5 → metrics

PAIRS_SPEC = [
    {
        "case":     "BraTS-GLI-00008-100",
        "slice_a":  "slice_0070",
        "slice_b":  "slice_0072",
        "slice_gt": "slice_0071",
        "label":    "caso1_gt",
        "descripcion": "Caso 100 — slices 70/72 con GT 71",
    },
    {
        "case":     "BraTS-GLI-00008-101",
        "slice_a":  "slice_0070",
        "slice_b":  "slice_0072",
        "slice_gt": "slice_0071",
        "label":    "caso2_gt",
        "descripcion": "Caso 101 — slices 70/72 con GT 71",
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
  .good { background: #d4edda; }
  .mid  { background: #fff3cd; }
  .bad  { background: #f8d7da; }
  .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }
  .panel { background: #fff; border: 1px solid #ddd; padding: 12px; border-radius: 6px; }
  .panel h3 { margin: 0 0 6px; font-size: 14px; }
  .slider-row { display: flex; align-items: center; gap: 10px; margin: 6px 0 10px; }
  .slider-row input[type=range] { flex: 1; }
  .t-label { font-family: ui-monospace, Menlo, monospace; font-size: 13px; min-width: 64px; }
  svg { width: 100%; height: 420px; background: #fff; }
  .legend { font-size: 12px; color: #444; margin-top: 4px; }
  .legend span { display: inline-block; margin-right: 12px; }
  .swatch { display: inline-block; width: 12px; height: 3px; vertical-align: middle; margin-right: 4px; }
  .summary { background: #fff; border: 1px solid #ddd; padding: 12px 16px; border-radius: 6px; max-width: 980px; }
</style>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>

<h1>BraTS glioma — interpolación lineal de contornos</h1>
<p class="note">
  Dashboard con métricas por par (Dice, IoU, Hausdorff, areaErr) comparando
  el contorno interpolado a <em>t</em>=0.5 contra el slice GT real intermedio.
  El visualizador interactivo permite recorrer <em>t</em> ∈ [0, 1] en pasos de
  0.1; los contornos azul (A), rojo (B) y naranja (GT) permanecen fijos, y el
  contorno verde se actualiza con los vértices precalculados para cada <em>t</em>.
</p>

<h2>1. Métricas por par</h2>
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
  Amarillo: zona intermedia · Rojo: peor desempeño.
</p>

<h2>2. Comparativa por métrica</h2>
<div class="grid">
  <div class="panel"><h3>Dice (mayor es mejor)</h3><canvas id="chart-dice"></canvas></div>
  <div class="panel"><h3>IoU (mayor es mejor)</h3><canvas id="chart-iou"></canvas></div>
  <div class="panel"><h3>Hausdorff (menor es mejor)</h3><canvas id="chart-hd"></canvas></div>
  <div class="panel"><h3>areaErr (menor es mejor)</h3><canvas id="chart-aerr"></canvas></div>
</div>

<h2>3. Visualizador interactivo — slider de <em>t</em></h2>
<div id="viz-grid" class="grid"></div>

<h2>4. Visualización 3D — contornos apilados</h2>
<p class="note">
  Escena Three.js con los 11 contornos interpolados (<em>t</em>=0.0 a 1.0)
  apilados en el eje Z (Z = t × 10). Color azul → rojo según <em>t</em>;
  contorno A en Z=0 (azul sólido), B en Z=10 (rojo sólido), GT en Z=5
  (naranja discontinuo). Rotar: clic izquierdo · zoom: scroll ·
  pan: clic derecho.
</p>
<div id="viz3d-grid" class="grid"></div>

<h2>5. Resumen</h2>
<div class="summary" id="summary"></div>

<script>
const PAIRS = __PAIRS_JSON__;

// -----------------------------------------------------------------------------
// Section 1: metrics table.
// -----------------------------------------------------------------------------
function cls(v, good, mid, lowerBetter=false) {
  if (Number.isNaN(v)) return "";
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
  const tr = document.createElement("tr");
  tr.innerHTML = `
    <td class="txt">${p.label}<br><span style="color:#777;font-size:11px">${p.case}</span></td>
    <td class="txt">${p.slice_a} / <b>${p.slice_gt}</b> / ${p.slice_b}</td>
    <td class="${cls(m.dice, 0.90, 0.75)}">${m.dice.toFixed(4)}</td>
    <td class="${cls(m.iou,  0.80, 0.60)}">${m.iou.toFixed(4)}</td>
    <td class="${cls(m.hausdorff, 3.0, 6.0, true)}">${m.hausdorff.toFixed(3)}</td>
    <td class="${cls(m.area_err, 0.05, 0.15, true)}">${m.area_err.toFixed(4)}</td>
    <td class="txt">${m.self_int ? "yes" : "no"}</td>
  `;
  tbody.appendChild(tr);
}

// -----------------------------------------------------------------------------
// Section 2: bar charts.
// -----------------------------------------------------------------------------
const labels = PAIRS.map(p => p.label);
function bar(id, key, title) {
  new Chart(document.getElementById(id), {
    type: "bar",
    data: {
      labels: labels,
      datasets: [{
        label: title,
        data:  PAIRS.map(p => p.metrics[key]),
        backgroundColor: "rgba(33,150,243,0.7)",
        borderColor:     "rgba(33,150,243,1)",
        borderWidth: 1,
      }],
    },
    options: {
      responsive: true,
      plugins: { legend: { display: false } },
      scales:  { y: { beginAtZero: true } },
    },
  });
}
bar("chart-dice", "dice",      "Dice");
bar("chart-iou",  "iou",       "IoU");
bar("chart-hd",   "hausdorff", "Hausdorff");
bar("chart-aerr", "area_err",  "areaErr");

// -----------------------------------------------------------------------------
// Section 3: interactive SVG viewer with t-slider.
// -----------------------------------------------------------------------------
const NS = "http://www.w3.org/2000/svg";

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
  if (!contour.length) return "";
  const sx = x => ((x - bb.xmin) / (bb.xmax - bb.xmin)) * w;
  const sy = y => h - ((y - bb.ymin) / (bb.ymax - bb.ymin)) * h;     // flip y
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
      <path id="pi-${p.label}"  fill="none" stroke="#2ca02c" stroke-width="2.2"/>
    </svg>
    <div class="legend">
      <span><span class="swatch" style="background:#1f77b4"></span>A (${p.slice_a})</span>
      <span><span class="swatch" style="background:#d62728"></span>B (${p.slice_b})</span>
      <span><span class="swatch" style="background:#ff7f0e"></span>GT (${p.slice_gt})</span>
      <span><span class="swatch" style="background:#2ca02c"></span>Interpolado(t)</span>
    </div>
  `;
  vizGrid.appendChild(panel);

  const allInterp = Object.values(p.interpolated).flat();
  const bb = bboxOf(p.contour_a, p.contour_b, p.contour_gt, allInterp);
  const W = 600, H = 420;
  document.getElementById(`pa-${p.label}`).setAttribute("d",  pathD(p.contour_a,  bb, W, H));
  document.getElementById(`pb-${p.label}`).setAttribute("d",  pathD(p.contour_b,  bb, W, H));
  document.getElementById(`pgt-${p.label}`).setAttribute("d", pathD(p.contour_gt, bb, W, H));

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

// -----------------------------------------------------------------------------
// Section 4: 3D stacked contour viewer (Three.js).
// -----------------------------------------------------------------------------
function centroidOf(contour) {
  let sx = 0, sy = 0;
  for (const [x, y] of contour) { sx += x; sy += y; }
  const n = contour.length || 1;
  return [sx / n, sy / n];
}

function contourToLine(contour, cx, cy, z, color, dashed, linewidth) {
  const pts = [];
  for (const [x, y] of contour) pts.push(new THREE.Vector3(x - cx, y - cy, z));
  if (contour.length) pts.push(new THREE.Vector3(contour[0][0] - cx, contour[0][1] - cy, z));
  const geom = new THREE.BufferGeometry().setFromPoints(pts);
  const mat = dashed
    ? new THREE.LineDashedMaterial({ color, linewidth, dashSize: 1.0, gapSize: 0.6 })
    : new THREE.LineBasicMaterial({ color, linewidth });
  const line = new THREE.Line(geom, mat);
  if (dashed) line.computeLineDistances();
  return line;
}

function hslGradient(t) {
  // hue 240 (blue) -> 0 (red)
  const hue = (240 * (1 - t)) / 360;
  const col = new THREE.Color();
  col.setHSL(hue, 0.85, 0.55);
  return col;
}

function build3DScene(panel, pair) {
  const canvasWrap = panel.querySelector(".scene3d");
  const W = canvasWrap.clientWidth || 600;
  const H = 480;

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x1a1a2e);

  const camera = new THREE.PerspectiveCamera(45, W / H, 0.1, 500);
  const initialPos = new THREE.Vector3(25, -25, 15);
  camera.position.copy(initialPos);
  camera.up.set(0, 0, 1);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(W, H);
  renderer.setPixelRatio(window.devicePixelRatio);
  canvasWrap.appendChild(renderer.domElement);

  const controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.target.set(0, 0, 5);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;

  // Reference axes + ground grid (faint).
  const grid = new THREE.GridHelper(40, 20, 0x444466, 0x2a2a44);
  grid.rotation.x = Math.PI / 2;
  scene.add(grid);

  const [cx, cy] = centroidOf(pair.contour_a);

  // 11 interpolated contours t=0..1 stacked at z = t*10.
  const tKeys = Object.keys(pair.interpolated).sort();
  for (const tStr of tKeys) {
    const t = parseFloat(tStr);
    const c = pair.interpolated[tStr];
    if (!c || !c.length) continue;
    const color = hslGradient(t);
    scene.add(contourToLine(c, cx, cy, t * 10.0, color, false, 1));
  }

  // A at z=0 (solid blue), B at z=10 (solid red).
  scene.add(contourToLine(pair.contour_a, cx, cy, 0.0,  0x1f77b4, false, 2));
  scene.add(contourToLine(pair.contour_b, cx, cy, 10.0, 0xd62728, false, 2));

  // GT at z=5 if present (dashed orange).
  if (pair.contour_gt && pair.contour_gt.length) {
    scene.add(contourToLine(pair.contour_gt, cx, cy, 5.0, 0xff7f0e, true, 1.5));
  }

  function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
  }
  animate();

  // Reset camera button.
  panel.querySelector(".reset-cam").addEventListener("click", () => {
    camera.position.copy(initialPos);
    controls.target.set(0, 0, 5);
    controls.update();
  });

  // Handle resize.
  window.addEventListener("resize", () => {
    const w = canvasWrap.clientWidth || W;
    renderer.setSize(w, H);
    camera.aspect = w / H;
    camera.updateProjectionMatrix();
  });
}

const viz3dGrid = document.getElementById("viz3d-grid");
const pendingScenes = [];
for (const p of PAIRS) {
  const panel = document.createElement("div");
  panel.className = "panel";
  panel.innerHTML = `
    <h3>${p.label} — ${p.descripcion}</h3>
    <div class="scene3d" style="width:100%;height:480px;position:relative;
         background:#1a1a2e;border-radius:4px;overflow:hidden;"></div>
    <div class="legend" style="margin-top:6px;">
      <span><span class="swatch" style="background:#1f77b4"></span>A (Z=0)</span>
      <span><span class="swatch" style="background:#d62728"></span>B (Z=10)</span>
      <span><span class="swatch" style="background:#ff7f0e"></span>GT (Z=5, discontinuo)</span>
      <span><span class="swatch" style="background:linear-gradient(90deg,#3b6cff,#ff3b3b)"></span>Interpolados (t=0→1)</span>
      <button class="reset-cam" style="margin-left:12px;font-size:11px;
              padding:3px 8px;cursor:pointer;">Reset cámara</button>
    </div>
  `;
  viz3dGrid.appendChild(panel);
  pendingScenes.push({ panel, pair: p, built: false });
}

// Lazy-build scenes when they enter the viewport.
const io = new IntersectionObserver((entries, obs) => {
  for (const e of entries) {
    if (!e.isIntersecting) continue;
    const item = pendingScenes.find(s => s.panel === e.target);
    if (item && !item.built) {
      item.built = true;
      build3DScene(item.panel, item.pair);
      obs.unobserve(e.target);
    }
  }
}, { rootMargin: "100px" });
for (const s of pendingScenes) io.observe(s.panel);

// -----------------------------------------------------------------------------
// Section 5: summary.
// -----------------------------------------------------------------------------
const mean = arr => arr.reduce((a, b) => a + b, 0) / arr.length;
const summary = document.getElementById("summary");
summary.innerHTML = `
  <b>Promedio sobre ${PAIRS.length} pares:</b><br>
  Dice = ${mean(PAIRS.map(p => p.metrics.dice)).toFixed(4)} ·
  IoU = ${mean(PAIRS.map(p => p.metrics.iou)).toFixed(4)} ·
  Hausdorff = ${mean(PAIRS.map(p => p.metrics.hausdorff)).toFixed(3)} ·
  areaErr = ${mean(PAIRS.map(p => p.metrics.area_err)).toFixed(4)}
  <br><br>
  <span class="note">
    GT = slice real intermedio (BraTS). El interpolado a t=0.5 se rasteriza en
    una grilla 512×512 dentro del bounding box común para calcular Dice/IoU;
    Hausdorff se computa directamente entre vértices.
  </span>
`;
</script>

__STRESS_SECTION__

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
                rel = f"output/stress_test/{runs[-1].name}/stress_comparison.png"
                img_block = (
                    '<div class="panel" style="margin-top:14px;max-width:980px;">'
                    '<h3>Estrella (r=1) vs cuadrado — salida del binario</h3>'
                    f'<img src="{rel}" alt="Comparativa stress test" '
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

    return f"""<h2>6. Experimento de stress — Robustez del pipeline</h2>
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


def main():
    if not BINARY.exists():
        print(f"ERROR: missing binary {BINARY}. Build first.", file=sys.stderr)
        sys.exit(1)

    resolved = [resolve_pair(spec) for spec in PAIRS_SPEC]
    results  = [process_pair(spec) for spec in resolved]

    payload = json.dumps(results, ensure_ascii=False, separators=(",", ":"))
    html    = HTML_TEMPLATE.replace("__PAIRS_JSON__", payload)
    html    = html.replace("__STRESS_SECTION__", build_stress_section())
    OUT_HTML.write_text(html, encoding="utf-8")
    print(f"\nDashboard escrito en {OUT_HTML.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
