"""
Corre múltiples pares de slices de BraTS, genera .obj interpolados y PNGs.

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


def first_two_slices(case):
    objs = sorted((CONTOURS / case).glob("*.obj"))
    if len(objs) < 2:
        return None
    return objs[0].stem, objs[1].stem


def find_extreme_pairs():
    """Primer y último slice de cada caso con al menos 3 slices."""
    pairs = []
    for case_dir in sorted(CONTOURS.iterdir()):
        if not case_dir.is_dir():
            continue
        slices = sorted(case_dir.glob("slice_*.obj"))
        if len(slices) >= 3:
            pairs.append({
                "case": case_dir.name,
                "slice_a": slices[0].stem,
                "slice_b": slices[-1].stem,
                "label": f"{case_dir.name}_extremos",
            })
    return pairs


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
    case = pair["case"]
    sa = pair["slice_a"]
    sb = pair["slice_b"]
    label = pair["label"]
    obj_a = Path(pair.get("obj_a", CONTOURS / case / f"{sa}.obj"))
    obj_b = Path(pair.get("obj_b", CONTOURS / case / f"{sb}.obj"))
    out_obj = OUTPUT / f"{label}_interpolated.obj"
    out_png = OUTPUT / f"{label}.png"

    print(f"=== {label}  ({case}: {sa} -> {sb}) ===")
    if not obj_a.exists() or not obj_b.exists():
        missing = [str(p) for p in (obj_a, obj_b) if not p.exists()]
        print(f"  [SKIP] archivos faltantes: {missing}")
        return {"label": label, "case": case, "sa": sa, "sb": sb,
                "status": "skipped", "reason": "missing"}

    # 1) Ejecutar el binario C++.
    proc = subprocess.run(
        [str(BINARY), str(obj_a), str(obj_b)],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    print(proc.stdout, end="")
    if proc.returncode != 0:
        print(proc.stderr, end="")
        return {"label": label, "case": case, "sa": sa, "sb": sb,
                "status": "error", "reason": "binary"}

    stats = _parse_stats(proc.stdout)

    # 2) Copiar el .obj resultante con el nombre del caso.
    src_obj = OUTPUT / "contour_interpolated.obj"
    shutil.copy2(src_obj, out_obj)

    # 3) Generar el PNG llamando a visualize_results como módulo.
    sys.argv = [
        "visualize_results.py",
        str(obj_a), str(obj_b), str(out_obj), str(out_png),
    ]
    viz.main()
    print()
    return {"label": label, "case": case, "sa": sa, "sb": sb,
            "status": "ok", "png": out_png, "obj": out_obj, **stats}


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)

    PAIRS = [
        {"case": "BraTS-GLI-00008-100",
         "slice_a": "slice_0070", "slice_b": "slice_0071",
         "label": "caso1_consec"},
        # Slices separados ~6 dentro del mismo caso (los 60/65/75/80 pedidos
        # originalmente no fueron extraídos por el pipeline de Python para
        # este caso; usamos los disponibles más distantes).
        {"case": "BraTS-GLI-00008-100",
         "slice_a": "slice_0066", "slice_b": "slice_0072",
         "label": "caso2_sep6"},
        # Slices muy separados en el caso 102 (rango 0048..0059, 11 de gap).
        {"case": "BraTS-GLI-00008-102",
         "slice_a": "slice_0048", "slice_b": "slice_0059",
         "label": "caso3_sep11"},
    ]

    # Casos 4 y 5: primeros dos slices con ET disponibles.
    for case, label in [("BraTS-GLI-00008-101", "caso4_nuevo"),
                        ("BraTS-GLI-00008-102", "caso5_nuevo")]:
        pair = first_two_slices(case)
        if pair is None:
            print(f"[WARN] {case} no tiene al menos dos .obj; se omite.")
            continue
        sa, sb = pair
        PAIRS.append({"case": case, "slice_a": sa, "slice_b": sb,
                      "label": label})

    # Pares extremos (primer vs último slice) de cada caso con >=3 slices.
    extreme = find_extreme_pairs()
    print(f"[INFO] Pares extremos detectados: {len(extreme)}")
    PAIRS.extend(extreme)

    viz = load_viz()
    results = [run_pair(p, viz) for p in PAIRS]

    # ------------------------------------------------------------------
    # Tabla resumen.
    # ------------------------------------------------------------------
    print("=" * 92)
    print("Resumen:")
    header = f"  {'label':<38} {'caso':<22} {'slices':<24} {'A':>4} {'B':>4} {'self-int':>9}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    ok_results = []
    for r in results:
        if r["status"] != "ok":
            print(f"  {r['label']:<38} {r['case']:<22} "
                  f"{r['sa']+'  '+r['sb']:<24} "
                  f"[{r['status']}]")
            continue
        ok_results.append(r)
        si = "yes" if r.get("self_int") else "no"
        print(f"  {r['label']:<38} {r['case']:<22} "
              f"{r['sa']+'  '+r['sb']:<24} "
              f"{str(r.get('va','?')):>4} {str(r.get('vb','?')):>4} {si:>9}")

    # ------------------------------------------------------------------
    # Detección de auto-intersecciones; si ninguno, probar pares cruzados.
    # ------------------------------------------------------------------
    intersected = [r for r in ok_results if r.get("self_int")]
    if intersected:
        print("\nPares con auto-intersecciones:")
        for r in intersected:
            print(f"  * {r['label']}  ({r['case']}: {r['sa']} -> {r['sb']})")
    else:
        print("\nAVISO: ningún par generó auto-intersecciones. "
              "Considerar pares entre casos distintos.")
        print("Probando pares cruzados (primer slice de X vs último de Y)...\n")

        cases = sorted(
            d for d in CONTOURS.iterdir()
            if d.is_dir() and len(sorted(d.glob("slice_*.obj"))) >= 2
        )
        cross_pairs = []
        for i, ca in enumerate(cases[:5]):
            cb = cases[(i + 1) % len(cases[:5])]
            if ca == cb:
                continue
            sa_path = sorted(ca.glob("slice_*.obj"))[0]
            sb_path = sorted(cb.glob("slice_*.obj"))[-1]
            cross_pairs.append({
                "case": f"{ca.name}__x__{cb.name}",
                "slice_a": sa_path.stem,
                "slice_b": sb_path.stem,
                "label": f"cross_{i:02d}_{ca.name}_vs_{cb.name}",
                "obj_a": sa_path,
                "obj_b": sb_path,
            })

        cross_results = [run_pair(p, viz) for p in cross_pairs]
        cross_ok = [r for r in cross_results if r["status"] == "ok"]
        cross_int = [r for r in cross_ok if r.get("self_int")]
        print("-" * 92)
        print("Resumen pares cruzados:")
        for r in cross_ok:
            si = "yes" if r.get("self_int") else "no"
            print(f"  {r['label']:<60} A={r.get('va','?'):>4} "
                  f"B={r.get('vb','?'):>4}  self-int={si}")
        if cross_int:
            print("\nPares cruzados con auto-intersecciones:")
            for r in cross_int:
                print(f"  * {r['label']}")
        else:
            print("\n(Ningún par cruzado generó auto-intersecciones tampoco.)")

        ok_results.extend(cross_ok)

    ok_pngs = [r["png"] for r in ok_results]
    if ok_pngs:
        print(f"\n{len(ok_pngs)} PNGs generados en {OUTPUT}/")


if __name__ == "__main__":
    main()
