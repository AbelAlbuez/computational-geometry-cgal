"""
Corre 4 pares representativos de slices BraTS, genera .obj interpolados y PNGs.

Uso:
    python scripts/run_examples.py

Debe correrse desde src/final-project-interpolation/ con el venv activado.
"""
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
BINARY = ROOT / "build" / "contour_interpolator"
CONTOURS = ROOT / "data" / "contours"
OUTPUT = ROOT / "output"
VIZ_SCRIPT = ROOT / "scripts" / "visualize_results.py"


def load_viz():
    spec = importlib.util.spec_from_file_location("visualize_results",
                                                  VIZ_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _parse_stats(stdout):
    stats = {"va": None, "vb": None, "vp": None, "self_int": None}
    for line in stdout.splitlines():
        ls = line.strip()
        if ls.startswith("Contour A:"):
            stats["va"] = int(ls.split(":")[1].split()[0])
        elif ls.startswith("Contour B:"):
            stats["vb"] = int(ls.split(":")[1].split()[0])
        elif ls.startswith("Interpolated contour:"):
            stats["vp"] = int(ls.split(":")[1].split()[0])
        elif ls.startswith("Self-intersections detected:"):
            stats["self_int"] = ls.split(":")[1].strip().lower() == "yes"
    return stats


def run_pair(pair, viz):
    label = pair["label"]
    obj_a = pair["obj_a"]
    obj_b = pair["obj_b"]
    out_obj = OUTPUT / f"{label}_interpolated.obj"
    out_png = OUTPUT / f"{label}.png"

    print(f"=== {label}  ({pair['descripcion']}) ===")
    print(f"  A: {obj_a.relative_to(ROOT)}")
    print(f"  B: {obj_b.relative_to(ROOT)}")

    proc = subprocess.run(
        [str(BINARY), str(obj_a), str(obj_b)],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    print(proc.stdout, end="")
    if proc.returncode != 0:
        print(proc.stderr, end="")
        return {**pair, "status": "error"}

    stats = _parse_stats(proc.stdout)
    shutil.copy2(OUTPUT / "contour_interpolated.obj", out_obj)

    sys.argv = [
        "visualize_results.py",
        str(obj_a), str(obj_b), str(out_obj), str(out_png),
    ]
    viz.main()
    print()
    return {**pair, "status": "ok", "png": out_png, "obj": out_obj, **stats}


def _resolve_pair(spec):
    """Resuelve los paths reales (obj_a, obj_b) y los nombres de slice."""
    label = spec["label"]
    sa_desc = spec["slice_a"]
    sb_desc = spec["slice_b"]

    if spec["label"] == "caso4_cruzado":
        # Cruzado: primer slice de dos casos distintos.
        case_a = "BraTS-GLI-00008-100"
        case_b = "BraTS-GLI-00008-101"
        obj_a = sorted((CONTOURS / case_a).glob("slice_*.obj"))[0]
        obj_b = sorted((CONTOURS / case_b).glob("slice_*.obj"))[0]
        return {
            "label": label,
            "case": f"{case_a} vs {case_b}",
            "descripcion": spec["descripcion"],
            "slice_a": obj_a.stem,
            "slice_b": obj_b.stem,
            "obj_a": obj_a,
            "obj_b": obj_b,
        }

    case = spec["case"]
    slices = sorted((CONTOURS / case).glob("slice_*.obj"))
    if sa_desc == "primer slice disponible":
        obj_a = slices[0]
        sa = obj_a.stem
    else:
        obj_a = CONTOURS / case / f"{sa_desc}.obj"
        sa = sa_desc
    if sb_desc == "último slice disponible":
        obj_b = slices[-1]
        sb = obj_b.stem
    else:
        obj_b = CONTOURS / case / f"{sb_desc}.obj"
        sb = sb_desc
    return {
        "label": label,
        "case": case,
        "descripcion": spec["descripcion"],
        "slice_a": sa,
        "slice_b": sb,
        "obj_a": obj_a,
        "obj_b": obj_b,
    }


PAIRS = [
    {
        "case": "BraTS-GLI-00008-100",
        "slice_a": "slice_0070",
        "slice_b": "slice_0071",
        "label": "caso1_consecutivos",
        "descripcion": "Slices consecutivos",
    },
    {
        "case": "BraTS-GLI-00008-100",
        "slice_a": "slice_0066",
        "slice_b": "slice_0072",
        "label": "caso2_separados",
        "descripcion": "Slices separados 6",
    },
    {
        "case": "BraTS-GLI-00008-101",
        "slice_a": "primer slice disponible",
        "slice_b": "último slice disponible",
        "label": "caso3_extremos",
        "descripcion": "Par extremo mismo caso",
    },
    {
        "case": "BraTS-GLI-00008-100 vs BraTS-GLI-00008-101",
        "slice_a": "primer slice de BraTS-GLI-00008-100",
        "slice_b": "primer slice de BraTS-GLI-00008-101",
        "label": "caso4_cruzado",
        "descripcion": "Par cruzado entre casos",
    },
]


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    viz = load_viz()
    resolved = [_resolve_pair(p) for p in PAIRS]
    results = [run_pair(p, viz) for p in resolved]

    print("=" * 92)
    print("Resumen:")
    header = (f"  {'label':<22} {'caso':<40} "
              f"{'A':>4} {'B':>4} {'self-int':>9}")
    print(header)
    print("  " + "-" * (len(header) - 2))
    for r in results:
        if r["status"] != "ok":
            print(f"  {r['label']:<22} {r['case']:<40} [error]")
            continue
        si = "yes" if r.get("self_int") else "no"
        print(f"  {r['label']:<22} {r['case']:<40} "
              f"{str(r.get('va','?')):>4} {str(r.get('vb','?')):>4} "
              f"{si:>9}")

    ok = [r for r in results if r["status"] == "ok"]
    print(f"\n{len(ok)} PNGs y {len(ok)} .obj generados en "
          f"{OUTPUT.relative_to(ROOT.parent.parent)}/")


if __name__ == "__main__":
    main()
