"""
Visualiza contorno A, contorno B y contorno interpolado en un PNG.

Uso:
    python visualize_results.py <obj_A> <obj_B> <obj_interpolado> <output_png>
"""
import os
import sys

import matplotlib.pyplot as plt


def read_obj(path):
    """Devuelve la lista de vértices (x, y) de un .obj con líneas 'v x y [z]'."""
    verts = []
    with open(path, "r") as fh:
        for line in fh:
            if not line or line[0] != "v":
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                x = float(parts[1])
                y = float(parts[2])
            except ValueError:
                continue
            verts.append((x, y))
    return verts


def close_loop(verts):
    if not verts:
        return [], []
    xs = [v[0] for v in verts] + [verts[0][0]]
    ys = [v[1] for v in verts] + [verts[0][1]]
    return xs, ys


def main():
    if len(sys.argv) < 5:
        print(
            "Uso: visualize_results.py <obj_A> <obj_B> "
            "<obj_interpolado> <output_png>"
        )
        sys.exit(1)

    path_a, path_b, path_i, out_png = sys.argv[1:5]

    a = read_obj(path_a)
    b = read_obj(path_b)
    p = read_obj(path_i)

    ax_x, ax_y = close_loop(a)
    bx_y, by_y = close_loop(b)
    px_x, py_y = close_loop(p)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(
        ax_x, ax_y,
        color="tab:blue", linestyle=":", linewidth=1.2,
        label=f"Slice A ({len(a)} v)",
    )
    ax.plot(
        bx_y, by_y,
        color="tab:red", linestyle=":", linewidth=1.2,
        label=f"Slice B ({len(b)} v)",
    )
    ax.plot(
        px_x, py_y,
        color="tab:green", linestyle="-", linewidth=2.2,
        label=f"Interpolado t=0.5 ({len(p)} v)",
    )
    ax.scatter(
        [v[0] for v in p], [v[1] for v in p],
        color="tab:green", s=10, zorder=3,
    )

    case = os.path.basename(os.path.dirname(path_a))
    name_a = os.path.splitext(os.path.basename(path_a))[0]
    name_b = os.path.splitext(os.path.basename(path_b))[0]
    ax.set_title(f"{case}\n{name_a}  →  {name_b}")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="best", fontsize=9)

    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"PNG guardado en {out_png}")


if __name__ == "__main__":
    main()
