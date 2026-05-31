"""
Stress test for SelfIntersectionResolver using SYNTHETIC contours.

Los contornos BraTS ET reales son demasiado convexos para producir
auto-intersecciones bajo interpolación lineal — además, la búsqueda de
rotación óptima dentro del binario neutraliza cualquier rotación cíclica
adversarial que apliquemos desde afuera. Por eso usamos un par patológico
controlado:

    A = estrella de 5 puntas (10 vértices, R=50, r=20)
    B = círculo suave        (10 vértices, R=40)

Ambos centrados en (60, 60). Las puntas de la estrella caen fuera del
círculo y las concavidades caen dentro, de modo que la correspondencia
vértice-a-vértice produce aristas que se cruzan en t=0.5 sin importar
qué rotación óptima encuentre el binario.

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

CENTER       = (60.0, 60.0)
STAR_R_OUT   = 50.0
STAR_R_IN    = 20.0
STAR_POINTS  = 5
CIRCLE_R     = 40.0
N_CIRCLE     = 10
N_STAR       = 2 * STAR_POINTS


# -----------------------------------------------------------------------------
# Contour generators.
# -----------------------------------------------------------------------------

def make_star(cx, cy, r_out, r_in, points=5):
    verts = []
    n = 2 * points
    for i in range(n):
        # ángulo inicial = -pi/2 para que la primera punta apunte arriba
        angle = -math.pi / 2.0 + i * (2.0 * math.pi / n)
        r     = r_out if (i % 2 == 0) else r_in
        verts.append((cx + r * math.cos(angle), cy + r * math.sin(angle)))
    return verts


def make_circle(cx, cy, r, n):
    verts = []
    for i in range(n):
        angle = -math.pi / 2.0 + i * (2.0 * math.pi / n)
        verts.append((cx + r * math.cos(angle), cy + r * math.sin(angle)))
    return verts


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
# Binary runner.
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
    return self_int, proc.stdout


# -----------------------------------------------------------------------------
# Plot helpers.
# -----------------------------------------------------------------------------

def _closed(verts):
    if not verts:
        return [], []
    xs = [v[0] for v in verts] + [verts[0][0]]
    ys = [v[1] for v in verts] + [verts[0][1]]
    return xs, ys


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

    # -- Step 1-2: build synthetic contours and dump as .obj.
    a = make_star(*CENTER, STAR_R_OUT, STAR_R_IN, STAR_POINTS)
    b = make_circle(*CENTER, CIRCLE_R, N_CIRCLE)
    obj_a = run_dir / "A.obj"
    obj_b = run_dir / "B.obj"
    write_obj_2d(obj_a, a, "synthetic 5-pointed star (R=50, r=20)")
    write_obj_2d(obj_b, b, "synthetic smooth circle (R=40)")

    # -- Step 3-4: interpolate at t = 0.5 and parse self-int flag.
    out_stress = run_dir / "stress.obj"
    self_int, _ = run_binary(obj_a, obj_b, out_stress, t=0.5)
    p_stress = read_obj_xy(out_stress)

    # -- Step 5: comparison figure.
    fig, (ax_in, ax_out) = plt.subplots(1, 2, figsize=(14, 7))

    ax_x, ax_y = _closed(a)
    bx_x, bx_y = _closed(b)
    ax_in.plot(ax_x, ax_y, color="tab:blue", linestyle=":", linewidth=1.6,
               marker="o", markersize=4,
               label=f"A estrella ({len(a)} v)")
    ax_in.plot(bx_x, bx_y, color="tab:red", linestyle=":", linewidth=1.6,
               marker="o", markersize=4,
               label=f"B círculo ({len(b)} v)")
    ax_in.set_title("Contornos sintéticos de entrada", fontsize=12)
    ax_in.set_aspect("equal", adjustable="datalim")
    ax_in.grid(True, linestyle="--", alpha=0.4)
    ax_in.legend(loc="best", fontsize=10)

    px_x, py_y = _closed(p_stress)
    ax_out.plot(px_x, py_y, color="tab:green", linestyle="-", linewidth=2.4,
                marker="o", markersize=5,
                label=f"Interpolado t=0.5 ({len(p_stress)} v)")
    ax_out.plot(ax_x, ax_y, color="tab:blue", linestyle=":", linewidth=0.8,
                alpha=0.45)
    ax_out.plot(bx_x, bx_y, color="tab:red", linestyle=":", linewidth=0.8,
                alpha=0.45)
    if self_int:
        ttl   = "Self-int detectadas: SÍ → resueltas por SelfIntersectionResolver"
        color = "tab:red"
    else:
        ttl   = "Self-int detectadas: NO"
        color = "tab:green"
    ax_out.set_title(ttl, color=color, fontsize=12)
    ax_out.set_aspect("equal", adjustable="datalim")
    ax_out.grid(True, linestyle="--", alpha=0.4)
    ax_out.legend(loc="best", fontsize=10)

    fig.suptitle(
        "Stress test SelfIntersectionResolver — par sintético estrella vs círculo",
        fontsize=13,
    )
    fig.tight_layout()
    png = run_dir / "stress_comparison.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)

    # -- Step 6: console summary.
    print(f"Vertices A (estrella)   : {len(a)}")
    print(f"Vertices B (circulo)    : {len(b)}")
    print(f"Self-int detectadas     : {'yes' if self_int else 'no'}")
    print(f"Output                  : {run_dir}")
    print(f"  - {obj_a.name}")
    print(f"  - {obj_b.name}")
    print(f"  - {out_stress.name}")
    print(f"  - {png.name}")
    print()
    print("NOTA: contornos sintéticos usados porque BraTS ET es demasiado convexo "
          "para producir auto-intersecciones bajo interpolación lineal con best_rotation.")


if __name__ == "__main__":
    main()
