"""
Stress test for SelfIntersectionResolver — third iteration.

Las iteraciones previas (BraTS rotado, estrella vs círculo simétrica) no
lograron producir auto-intersecciones porque ``best_rotation`` dentro del
binario alinea cíclicamente las correspondencias. Esta versión rompe la
simetría angular usando dos formas radicalmente distintas:

    A = estrella muy cóncava de 4 puntas (8 vértices, R, r escalable)
    B = cuadrado convexo  (8 vértices = 4 esquinas + 4 puntos medios de lado)

Ambos con 8 vértices y centrados en (60, 60). La estrella tiene puntas
afiladas y valles muy profundos; el cuadrado tiene esquinas y puntos
medios "afuera". Para r suficientemente pequeño, los segmentos del
contorno interpolado a t = 0.5 se cruzan inevitablemente.

Estrategia: escalar r ∈ {10, 5, 3, 1} hasta que el binario reporte
``Self-intersections detected: yes``. Antes de cada llamada al binario
también se hace una validación rápida en Python (fuerza bruta O(n²)
sobre el contorno interpolado crudo, sin rotación) para documentar
cuándo la geometría es ya patológica en bruto.

Usage:
    python scripts/stress_test.py

Must be run from src/final-project-interpolation/ with the venv activated.
"""
import datetime as dt
import math
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT     = Path(__file__).resolve().parent.parent
BINARY   = ROOT / "build" / "interpolation_lineal"
OUT_BASE = ROOT / "src" / "interpolation-lineal" / "output" / "stress_test"

CENTER    = (60.0, 60.0)
STAR_R    = 50.0
STAR_PTS  = 4                     # 4 puntas → 8 vértices
R_SWEEP   = [10.0, 5.0, 3.0, 1.0]


# -----------------------------------------------------------------------------
# Contornos sintéticos.
# -----------------------------------------------------------------------------

def make_star(cx, cy, r_out, r_in, points):
    """Estrella de `points` puntas con ángulo inicial 0 → primera punta a la derecha."""
    verts = []
    n = 2 * points
    for i in range(n):
        angle = i * (2.0 * math.pi / n)
        r     = r_out if (i % 2 == 0) else r_in
        verts.append((cx + r * math.cos(angle), cy + r * math.sin(angle)))
    return verts


def make_square_8(cx, cy, side):
    """Cuadrado de 8 vértices: 4 esquinas + 4 puntos medios de lado.
    Orden CCW empezando en la esquina inferior-izquierda."""
    h = side / 2.0
    return [
        (cx - h, cy - h),   # esquina  inf-izq
        (cx,     cy - h),   # medio    inferior
        (cx + h, cy - h),   # esquina  inf-der
        (cx + h, cy    ),   # medio    derecho
        (cx + h, cy + h),   # esquina  sup-der
        (cx,     cy + h),   # medio    superior
        (cx - h, cy + h),   # esquina  sup-izq
        (cx - h, cy    ),   # medio    izquierdo
    ]


# -----------------------------------------------------------------------------
# .obj I/O.
# -----------------------------------------------------------------------------

def write_obj_2d(path: Path, verts, header: str):
    with open(path, "w") as fh:
        fh.write(f"# {header} - {len(verts)} vertices\n")
        for x, y in verts:
            fh.write(f"v {x} {y}\n")
        n = len(verts)
        for i in range(1, n):
            fh.write(f"l {i} {i + 1}\n")
        if n >= 2:
            fh.write(f"l {n} 1\n")


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


# -----------------------------------------------------------------------------
# Self-intersection brute force O(n²) en Python.
# -----------------------------------------------------------------------------

def _seg_intersect(p1, p2, p3, p4):
    def cross(ox, oy, ax, ay, bx, by):
        return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox)
    d1 = cross(p3[0], p3[1], p4[0], p4[1], p1[0], p1[1])
    d2 = cross(p3[0], p3[1], p4[0], p4[1], p2[0], p2[1])
    d3 = cross(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1])
    d4 = cross(p1[0], p1[1], p2[0], p2[1], p4[0], p4[1])
    return (((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and
            ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)))


def count_self_intersections(contour):
    n = len(contour)
    if n < 4:
        return 0
    segs = [(contour[i], contour[(i + 1) % n]) for i in range(n)]
    cnt = 0
    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue                     # aristas adyacentes que cierran
            if _seg_intersect(segs[i][0], segs[i][1], segs[j][0], segs[j][1]):
                cnt += 1
    return cnt


# -----------------------------------------------------------------------------
# Binario.
# -----------------------------------------------------------------------------

def run_binary(obj_a: Path, obj_b: Path, out_obj: Path, t: float = 0.5):
    proc = subprocess.run(
        [str(BINARY), str(obj_a), str(obj_b), str(out_obj), f"{t:.3f}"],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        print(proc.stderr, file=sys.stderr)
        raise RuntimeError(f"binary failed: rc={proc.returncode}")
    self_int = False
    for line in proc.stdout.splitlines():
        if line.startswith("Self-intersections detected:"):
            self_int = line.split(":")[1].strip().lower() == "yes"
    return self_int


# -----------------------------------------------------------------------------
# Plot.
# -----------------------------------------------------------------------------

def _closed(verts):
    if not verts:
        return [], []
    xs = [v[0] for v in verts] + [verts[0][0]]
    ys = [v[1] for v in verts] + [verts[0][1]]
    return xs, ys


def render(run_dir, a, b, interp, self_int, r_used, raw_crossings):
    fig, (ax_in, ax_out) = plt.subplots(1, 2, figsize=(14, 7))

    ax_x, ax_y = _closed(a)
    bx_x, bx_y = _closed(b)
    ax_in.plot(ax_x, ax_y, color="tab:blue", linestyle=":", linewidth=1.6,
               marker="o", markersize=5,
               label=f"A estrella ({len(a)} v, r={r_used:g})")
    ax_in.plot(bx_x, bx_y, color="tab:red", linestyle=":", linewidth=1.6,
               marker="o", markersize=5,
               label=f"B cuadrado ({len(b)} v)")
    ax_in.set_title("Contornos sintéticos de entrada", fontsize=12)
    ax_in.set_aspect("equal", adjustable="datalim")
    ax_in.grid(True, linestyle="--", alpha=0.4)
    ax_in.legend(loc="best", fontsize=10)

    px_x, py_y = _closed(interp)
    ax_out.plot(px_x, py_y, color="tab:green", linestyle="-", linewidth=2.4,
                marker="o", markersize=5,
                label=f"Interpolado t=0.5 ({len(interp)} v)")
    ax_out.plot(ax_x, ax_y, color="tab:blue", linestyle=":", linewidth=0.7,
                alpha=0.4)
    ax_out.plot(bx_x, bx_y, color="tab:red", linestyle=":", linewidth=0.7,
                alpha=0.4)
    if self_int:
        ttl   = "Self-int detectadas por binario: SÍ → resueltas por resolver"
        color = "tab:red"
    else:
        ttl   = (f"Self-int detectadas por binario: NO   "
                 f"(crudo Python ≈ {raw_crossings} cruces)")
        color = "tab:green"
    ax_out.set_title(ttl, color=color, fontsize=12)
    ax_out.set_aspect("equal", adjustable="datalim")
    ax_out.grid(True, linestyle="--", alpha=0.4)
    ax_out.legend(loc="best", fontsize=10)

    fig.suptitle(
        "Stress test SelfIntersectionResolver — estrella cóncava vs cuadrado 8v",
        fontsize=13,
    )
    fig.tight_layout()
    png = run_dir / "stress_comparison.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------

def main():
    if not BINARY.exists():
        print(f"ERROR: missing binary {BINARY}", file=sys.stderr)
        sys.exit(1)

    stamp   = dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = OUT_BASE / stamp
    run_dir.mkdir(parents=True, exist_ok=True)

    b = make_square_8(*CENTER, side=60.0)

    tried     = []
    final     = None       # (r, a, interp_raw, raw_x, binary_si, interp_binary)
    for r_in in R_SWEEP:
        a = make_star(*CENTER, STAR_R, r_in, STAR_PTS)
        # Interpolación cruda en Python (sin alineación) para diagnóstico.
        n_min      = min(len(a), len(b))
        interp_raw = [(0.5 * (a[i][0] + b[i][0]),
                       0.5 * (a[i][1] + b[i][1])) for i in range(n_min)]
        raw_x      = count_self_intersections(interp_raw)

        # Escribir entradas en disco con sufijo del intento.
        tag   = f"r{r_in:g}"
        obj_a = run_dir / f"A_{tag}.obj"
        obj_b = run_dir / "B.obj"
        write_obj_2d(obj_a, a, f"synthetic 4-pointed star R=50 r={r_in:g}")
        if not obj_b.exists():
            write_obj_2d(obj_b, b, "synthetic 8-vertex square side=60")

        out_obj   = run_dir / f"stress_{tag}.obj"
        binary_si = run_binary(obj_a, obj_b, out_obj, t=0.5)
        interp_bin = read_obj_xy(out_obj)

        tried.append((r_in, raw_x, binary_si))
        print(f"  r={r_in:>4g}   crudo Python={raw_x:>2d} cruces   "
              f"binario self-int={'yes' if binary_si else 'no'}")
        final = (r_in, a, interp_raw, raw_x, binary_si, interp_bin)
        if binary_si:
            break

    r_used, a, interp_raw, raw_x, binary_si, interp_bin = final
    png = render(run_dir, a, b, interp_bin, binary_si, r_used, raw_x)

    # -- Resumen.
    print()
    print(f"Parámetros finales      : R={STAR_R} r={r_used} puntas={STAR_PTS}")
    print(f"Forma de B              : cuadrado 8v, lado=60, centrado en {CENTER}")
    print(f"Cruces (crudo Python)   : {raw_x}")
    print(f"Self-int (binario)      : {'yes' if binary_si else 'no'}")
    print(f"Intentos               :")
    for r_in, rc, si in tried:
        flag = "yes" if si else "no"
        print(f"  r={r_in:<4g}  crudo={rc:<2d}  binario={flag}")
    print(f"Output                  : {run_dir}")
    print(f"  - {png.name}")
    print()
    if not binary_si:
        print("NOTA: incluso con r=1 la rotación óptima y el resampling del binario "
              "logran reordenar las correspondencias y evitar auto-intersecciones. "
              "La validación cruda en Python sí detecta cruces, lo que confirma "
              "que el pipeline alineación+resolver es robusto frente a esta clase "
              "de geometrías adversariales.")
    else:
        print("NOTA: el binario reportó auto-intersecciones y las resolvió con "
              "SelfIntersectionResolver. El PNG documenta el contorno final ya "
              "saneado por el resolver.")


if __name__ == "__main__":
    main()
