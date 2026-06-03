#!/usr/bin/env python3
# =============================================================================
# Genera figuras estáticas (PNG) para el informe LaTeX:
#   - loo_cdf.png : CDF empírica de Dice por método
#   - loo_box.png : distribución (box) de Dice por método
#   - interp_methods_<caso>.png : estilo "Abel" en 2D — para un corte oculto k,
#     se dibujan los vecinos A=k-1 y B=k+1, el corte real k y la predicción de
#     cada método (lineal / spline / sdf) en un panel 1x3.
# =============================================================================
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import evaluate as ev  # noqa: E402  (reutiliza list_slices/read_obj/predict_*)

M_COLOR = {"linear": "#1f77b4", "spline": "#ff7f0e", "sdf": "#9467bd"}
M_ES = {"linear": "lineal", "spline": "spline", "sdf": "sdf"}


def closed(c):
    return np.vstack([c, c[0]]) if len(c) else c


def load_raw(res):
    rows = []
    for fp in glob.glob(str(res / "loo_*" / "results_raw.csv")):
        for r in csv.DictReader(open(fp)):
            rows.append(r)
    return rows


def fig_cdf(res, out):
    raw = load_raw(res)
    plt.figure(figsize=(6, 4))
    for m in ("linear", "spline", "sdf"):
        d = np.sort([float(r["dice"]) for r in raw if r["method"] == m])
        if len(d):
            plt.plot(d, np.arange(1, len(d) + 1) / len(d), color=M_COLOR[m],
                     lw=2, label=M_ES[m])
    plt.axvline(0.8, ls="--", color="gray")
    plt.xlim(0, 1); plt.ylim(0, 1)
    plt.xlabel("Dice"); plt.ylabel("fracción de cortes ≤ Dice")
    plt.title("Distribución de Dice (CDF empírica)")
    plt.legend(loc="upper left"); plt.grid(alpha=.3)
    plt.tight_layout(); plt.savefig(out, dpi=140); plt.close()
    print("wrote", out)


def fig_box(res, out):
    raw = load_raw(res)
    data = [[float(r["dice"]) for r in raw if r["method"] == m]
            for m in ("linear", "spline", "sdf")]
    plt.figure(figsize=(6, 4))
    bp = plt.boxplot(data, labels=[M_ES[m] for m in ("linear", "spline", "sdf")],
                     patch_artist=True, showmeans=True)
    for patch, m in zip(bp["boxes"], ("linear", "spline", "sdf")):
        patch.set_facecolor(M_COLOR[m]); patch.set_alpha(.55)
    plt.ylabel("Dice"); plt.title("Distribución de Dice por método")
    plt.grid(alpha=.3, axis="y"); plt.tight_layout()
    plt.savefig(out, dpi=140); plt.close()
    print("wrote", out)


def fig_interp(binp, case_dir, out, k=None):
    slices = ev.list_slices(case_dir)
    if len(slices) < 5:
        print("skip (few slices):", case_dir); return False
    if k is None:
        k = len(slices) // 2
    z0, f0 = slices[k - 1]; zt, ft = slices[k]; z2, f2 = slices[k + 1]
    A, B, T = ev.read_obj(f0), ev.read_obj(f2), ev.read_obj(ft)
    t = (zt - z0) / (z2 - z0) if z2 > z0 else 0.5
    preds = {}
    with tempfile.TemporaryDirectory() as tmp:
        for m in ("linear", "sdf"):
            p = os.path.join(tmp, f"{m}.obj")
            ev.predict_pairwise(binp, m, f0, f2, p, t)
            preds[m] = ev.read_obj(p)
        try:
            sp = ev.predict_spline(binp, slices, k, os.path.join(tmp, "spl"))
            preds["spline"] = ev.read_obj(sp)
        except Exception:
            preds["spline"] = np.empty((0, 2))

    fig, axs = plt.subplots(1, 3, figsize=(11, 4.0))
    case = os.path.basename(case_dir)
    for ax, m in zip(axs, ("linear", "spline", "sdf")):
        ax.plot(*closed(A).T, color="#888", ls="--", lw=1, label=f"vecino A (z={int(z0)})")
        ax.plot(*closed(B).T, color="#bbb", ls="-.", lw=1, label=f"vecino B (z={int(z2)})")
        ax.plot(*closed(T).T, color="black", ls=":", lw=1.4, label=f"real (z={int(zt)})")
        if len(preds[m]):
            ax.plot(*closed(preds[m]).T, color=M_COLOR[m], lw=2.2, label=f"{M_ES[m]} (t={t:.2f})")
        ax.set_title(M_ES[m]); ax.set_aspect("equal"); ax.invert_yaxis()
        ax.tick_params(labelsize=8)
    axs[0].legend(loc="upper left", fontsize=8, framealpha=.9)
    fig.suptitle(f"Interpolación del corte oculto — {case}", y=1.02)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close()
    print("wrote", out)
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--contours", required=True, help="data/contours root")
    ap.add_argument("--out", required=True, help="output dir for PNGs")
    args = ap.parse_args()
    res, out = Path(args.results), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    fig_cdf(res, str(out / "loo_cdf.png"))
    fig_box(res, str(out / "loo_box.png"))
    for case in ("BraTS-GLI-00528-101", "BraTS-GLI-02066-105"):
        cdir = Path(args.contours) / case
        if cdir.is_dir():
            fig_interp(args.bin, str(cdir), str(out / f"interp_methods_{case}.png"))
        else:
            print("missing contours:", cdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
