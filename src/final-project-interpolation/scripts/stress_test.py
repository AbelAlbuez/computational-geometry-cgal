"""
Stress test for SelfIntersectionResolver.

Forces a self-intersection at t = 0.5 by rotating B by n//2 vertices
(so the optimal-rotation search in the binary cannot recover the
original correspondence). Runs the binary twice (normal vs stress),
parses stdout for the self-intersection flag and renders a side-by-side
comparison PNG.

Usage:
    python scripts/stress_test.py

Must be run from src/final-project-interpolation/ with the venv activated.
"""
import datetime as dt
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT     = Path(__file__).resolve().parent.parent
BINARY   = ROOT / "build" / "interpolation_lineal"
OBJ_A    = ROOT / "data" / "contours" / "BraTS-GLI-00008-100" / "slice_0070.obj"
OBJ_B    = ROOT / "data" / "contours" / "BraTS-GLI-00008-100" / "slice_0072.obj"
OUT_BASE = ROOT / "src" / "interpolation-lineal" / "output" / "stress_test"


# -----------------------------------------------------------------------------
# .obj I/O.
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


def write_obj_2d(path: Path, verts):
    with open(path, "w") as fh:
        fh.write(f"# stress test - rotated contour, {len(verts)} vertices\n")
        for x, y in verts:
            fh.write(f"v {x} {y}\n")
        n = len(verts)
        for i in range(1, n):
            fh.write(f"l {i} {i + 1}\n")
        if n >= 2:
            fh.write(f"l {n} 1\n")


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
# Plot.
# -----------------------------------------------------------------------------

def _closed(verts):
    if not verts:
        return [], []
    xs = [v[0] for v in verts] + [verts[0][0]]
    ys = [v[1] for v in verts] + [verts[0][1]]
    return xs, ys


def plot_pair(ax, a, b, p, label_b, title, color_title):
    ax_x, ax_y = _closed(a)
    bx_x, bx_y = _closed(b)
    px_x, py_y = _closed(p)
    ax.plot(ax_x, ax_y, color="tab:blue", linestyle=":", linewidth=1.3,
            label=f"A ({len(a)} v)")
    ax.plot(bx_x, bx_y, color="tab:red",  linestyle=":", linewidth=1.3,
            label=f"{label_b} ({len(b)} v)")
    ax.plot(px_x, py_y, color="tab:green", linestyle="-", linewidth=2.2,
            label=f"Interpolado t=0.5 ({len(p)} v)")
    ax.set_title(title, color=color_title, fontsize=12)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="best", fontsize=9)


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------

def main():
    if not BINARY.exists():
        print(f"ERROR: missing binary {BINARY}", file=sys.stderr)
        sys.exit(1)
    for p in (OBJ_A, OBJ_B):
        if not p.exists():
            print(f"ERROR: missing input {p}", file=sys.stderr)
            sys.exit(1)

    stamp   = dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = OUT_BASE / stamp
    run_dir.mkdir(parents=True, exist_ok=True)

    # -- Step 1: read B, compute rotation = n // 2 and write rotated copy.
    b_verts  = read_obj_xy(OBJ_B)
    n        = len(b_verts)
    rotation = n // 2
    b_rot    = b_verts[rotation:] + b_verts[:rotation]
    b_rot_path = run_dir / "B_rotado.obj"
    write_obj_2d(b_rot_path, b_rot)

    # -- Step 2: two binary runs.
    out_normal = run_dir / "normal.obj"
    out_stress = run_dir / "stress.obj"
    si_normal, _ = run_binary(OBJ_A, OBJ_B,     out_normal)
    si_stress, _ = run_binary(OBJ_A, b_rot_path, out_stress)

    # -- Step 3: comparison PNG.
    a       = read_obj_xy(OBJ_A)
    b       = read_obj_xy(OBJ_B)
    p_norm  = read_obj_xy(out_normal)
    p_stress = read_obj_xy(out_stress)

    fig, (ax_n, ax_s) = plt.subplots(1, 2, figsize=(14, 7))
    ttl_n = ("Normal — sin self-int" if not si_normal
             else "Normal — SELF-INT detectada")
    ttl_s = ("Stress (B rotado n/2) — sin self-int" if not si_stress
             else "Stress (B rotado n/2) — SELF-INT detectada")
    color_n = "tab:green" if not si_normal else "tab:red"
    color_s = "tab:green" if not si_stress else "tab:red"
    plot_pair(ax_n, a, b,     p_norm,   "B",         ttl_n, color_n)
    plot_pair(ax_s, a, b_rot, p_stress, "B (rotado)", ttl_s, color_s)
    fig.suptitle(
        f"Stress test SelfIntersectionResolver — "
        f"BraTS-GLI-00008-100  (rotation = n//2 = {rotation})",
        fontsize=13,
    )
    fig.tight_layout()
    png = run_dir / "stress_comparison.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)

    # -- Step 4: console summary.
    print(f"n vertices B            : {n}")
    print(f"rotation aplicada       : {rotation}")
    print(f"Self-int corrida normal : {'yes' if si_normal else 'no'}")
    print(f"Self-int corrida stress : {'yes' if si_stress else 'no'}")
    print(f"Outputs                 : {run_dir}")
    print(f"  - {b_rot_path.name}")
    print(f"  - {out_normal.name}")
    print(f"  - {out_stress.name}")
    print(f"  - {png.name}")


if __name__ == "__main__":
    main()
