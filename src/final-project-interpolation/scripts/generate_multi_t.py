"""
Generate 11 interpolated contours (t = 0.0 .. 1.0, step 0.1) with REAL z
coordinates so they can be stacked into a 3D volume in VTK / ParaView.

Usage:
    python scripts/generate_multi_t.py <obj_A> <obj_B> <z_a> <z_b>

Example:
    python scripts/generate_multi_t.py \\
        data/contours/BraTS-GLI-00008-100/slice_0070.obj \\
        data/contours/BraTS-GLI-00008-100/slice_0072.obj \\
        70 72

Must be run from src/final-project-interpolation/ with the venv activated.
"""
import csv
import datetime as dt
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d proj)


ROOT     = Path(__file__).resolve().parent.parent
BINARY   = ROOT / "build" / "interpolation_lineal"
OUT_BASE = ROOT / "src" / "interpolation-lineal" / "output" / "multi_t"

T_STEPS  = [round(i / 10.0, 1) for i in range(11)]


def read_obj_xy(path: Path):
    verts, edges = [], []
    with open(path, "r") as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == "v":
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        verts.append((float(parts[1]), float(parts[2])))
                    except ValueError:
                        pass
            elif line[0] == "l":
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        edges.append((int(parts[1]), int(parts[2])))
                    except ValueError:
                        pass
    return verts, edges


def rewrite_obj_with_z(path_in: Path, path_out: Path, z: float, t: float):
    verts, edges = read_obj_xy(path_in)
    with open(path_out, "w") as fh:
        fh.write(f"# contorno t={t:.2f} z={z:.4f}\n")
        for x, y in verts:
            fh.write(f"v {x} {y} {z:.4f}\n")
        if edges:
            for i, j in edges:
                fh.write(f"l {i} {j}\n")
        else:
            n = len(verts)
            for i in range(1, n):
                fh.write(f"l {i} {i + 1}\n")
            if n >= 2:
                fh.write(f"l {n} 1\n")
    return verts


def run_binary(obj_a: Path, obj_b: Path, out_obj: Path, t: float):
    proc = subprocess.run(
        [str(BINARY), str(obj_a), str(obj_b), str(out_obj), f"{t:.3f}"],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        return None, proc.stderr
    self_int = False
    for line in proc.stdout.splitlines():
        if line.startswith("Self-intersections detected:"):
            self_int = line.split(":")[1].strip().lower() == "yes"
    return self_int, proc.stdout


def main():
    if len(sys.argv) != 5:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    obj_a = Path(sys.argv[1]).resolve()
    obj_b = Path(sys.argv[2]).resolve()
    z_a   = float(sys.argv[3])
    z_b   = float(sys.argv[4])

    if not BINARY.exists():
        print(f"ERROR: missing binary {BINARY}", file=sys.stderr)
        sys.exit(1)
    for p in (obj_a, obj_b):
        if not p.exists():
            print(f"ERROR: missing input {p}", file=sys.stderr)
            sys.exit(1)

    stamp   = dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = OUT_BASE / stamp
    run_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    contours_3d = []   # for the preview
    for t in T_STEPS:
        tmp_obj = run_dir / f"_tmp_t{int(round(t * 100)):03d}.obj"
        si, msg = run_binary(obj_a, obj_b, tmp_obj, t)
        if si is None:
            print(f"WARN: binary failed for t={t:.2f}: {msg.strip()}",
                  file=sys.stderr)
            tmp_obj.unlink(missing_ok=True)
            continue

        z_real    = z_a + t * (z_b - z_a)
        final_obj = run_dir / f"contour_t{t:.2f}_z{z_real:.1f}.obj"
        verts     = rewrite_obj_with_z(tmp_obj, final_obj, z_real, t)
        tmp_obj.unlink(missing_ok=True)

        rows.append({
            "t":          f"{t:.2f}",
            "z":          f"{z_real:.4f}",
            "filename":   final_obj.name,
            "n_vertices": len(verts),
            "self_int":   "yes" if si else "no",
        })
        contours_3d.append((t, z_real, verts))
        print(f"  t={t:.2f}  z={z_real:7.3f}  n={len(verts):3d}  "
              f"self-int={'yes' if si else 'no'}  → {final_obj.name}")

    # -- index.csv
    csv_path = run_dir / "index.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["t", "z", "filename", "n_vertices", "self_int"])
        w.writeheader()
        w.writerows(rows)

    # -- 3D preview
    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111, projection="3d")
    cmap = cm.get_cmap("RdYlGn_r")
    for t, z_real, verts in contours_3d:
        if not verts:
            continue
        xs = [v[0] for v in verts] + [verts[0][0]]
        ys = [v[1] for v in verts] + [verts[0][1]]
        zs = [z_real] * len(xs)
        ax.plot(xs, ys, zs, color=cmap(t), linewidth=1.4,
                label=f"t={t:.1f} z={z_real:.1f}")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("z (mm)")
    ax.set_title(f"Multi-t stack ({len(contours_3d)} contornos)\n"
                 f"{obj_a.name}  z={z_a}  →  {obj_b.name}  z={z_b}")
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.05, 1.0))
    png_path = run_dir / "preview_3d.png"
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)

    print(f"\nGeneradas {len(rows)} capas en {run_dir}")
    print(f"index   : {csv_path.name}")
    print(f"preview : {png_path.name}")


if __name__ == "__main__":
    main()
