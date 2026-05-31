# RUN_LOCAL.md — Building & Running the Final Project on This Machine

> Audited, fixed, and **verified end-to-end on 2026-05-31**. The full pipeline
> (contour extraction → 3 interpolations → Poisson reconstruction → metrics →
> leave-one-slice-out) runs to completion. The **only** thing not present is the
> real BraTS 2024 GLI dataset — a synthetic BraTS-format generator stands in so
> everything downstream is exercised; swap in real data via `BRATS_DIR`.

## 1. Detected environment

| Item | Detected |
|------|----------|
| Host OS | Windows 11 Pro (10.0.26200) |
| **Build/run environment** | **WSL2 — Ubuntu 24.04.3 LTS** (kernel 6.6.87.2). The Windows host has no C++ toolchain; all building/running happens inside WSL. |
| CGAL | **5.6** (`libcgal-dev 5.6-1build3`, header-only) — *not* 6.x |
| Eigen | system: absent; **vendored** at `~/opt/eigen3/usr/include/eigen3` (no-sudo) |
| Compiler | g++ **13.3.0** |
| CMake / Make | **3.28.3** / GNU Make 4.3 |
| Python (WSL) | **3.12.3** |
| Python (Windows host) | 3.11.9 — not used for the build |
| BraTS dataset | not present (synthetic substitute provided) |

## 2. Dependency status (after setup)

| Dependency | Required for | Status |
|-----------|--------------|--------|
| CGAL ≥ 5.6 | C++ core | ✅ 5.6 |
| Eigen3 | Poisson + jet normals | ✅ vendored (`-DEIGEN3_INCLUDE_DIR`) |
| g++ (C++17) / CMake | build | ✅ 13.3.0 / 3.28.3 |
| numpy / matplotlib | evaluate.py | ✅ 2.4.6 / 3.10.9 |
| nibabel / scikit-image / tqdm | extract_contours.py | ✅ installed (5.4.2 / 0.26.0 / 4.67.3) |
| BraTS 2024 GLI data | extraction | ❌ not present → synthetic stand-in |

## 3. Install commands for dependencies

### This machine (WSL / Ubuntu) — what was done
```bash
# Python deps (no sudo; Ubuntu 24.04 Python is PEP-668 externally-managed)
python3 -m pip install --user --break-system-packages \
        nibabel scikit-image tqdm numpy matplotlib
# …or use an isolated venv instead:
sudo apt install -y python3-venv                 # one-time, needs sudo
python3 -m venv venv && . venv/bin/activate && pip install -r requirements.txt

# CGAL is already present. If ever missing:  sudo apt install -y libcgal-dev
# Eigen vendored without sudo:               bash scripts/setup_eigen.sh
```

### Windows-native alternative (vcpkg) — only if abandoning WSL
```powershell
git clone https://github.com/microsoft/vcpkg
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg\vcpkg install cgal eigen3
cmake -S src\final-project-interpolation -B build `
  -DCMAKE_TOOLCHAIN_FILE=C:/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
# also needs MSVC Build Tools + a Windows Python with the pip deps above
```

### macOS (Homebrew) equivalent
```bash
brew install cgal eigen cmake        # CGAL pulls boost, gmp, mpfr
python3 -m venv venv && . venv/bin/activate && pip install -r requirements.txt
```

## 4. Verified command sequence

All commands run **inside WSL**. Repo (from WSL):
`/mnt/c/Users/LENOVO LOQ/Desktop/geometria_computacional/taller1/computational-geometry-cgal/src/final-project-interpolation`

### One-shot wrapper (recommended) — ✅ VERIFIED end-to-end
```bash
cd <repo>/src/final-project-interpolation
# Synthetic data (no BraTS needed):
bash scripts/run_pipeline.sh
# Real data:
BRATS_DIR=/path/to/BraTS  bash scripts/run_pipeline.sh
```
This builds (if needed), extracts contours, reconstructs (linear + spline),
prints the metrics, and runs the leave-one-slice-out study.

### Or step by step

**(a) Configure + build out-of-source** — ✅ VERIFIED (in-source builds are forbidden)
```bash
cmake -S <repo>/src/final-project-interpolation -B ~/builds/contour \
      -DCMAKE_BUILD_TYPE=Release \
      -DEIGEN3_INCLUDE_DIR=$HOME/opt/eigen3/usr/include/eigen3
cmake --build ~/builds/contour -j
ctest --test-dir ~/builds/contour --output-on-failure        # 5/5 pass
```

**(b) Contour extraction on a small subset** — ✅ VERIFIED
```bash
# Synthetic BraTS-format dataset (only substitute for the real download):
python3 scripts/make_synthetic_case.py --out ~/synth_brats --cases 2 --spacing 1 1 1
python3 scripts/extract_contours.py --dataset ~/synth_brats --cases 2
#   -> data/contours/<case>/slice_XXXX.obj  +  data/contours/index.csv
#      (index.csv now carries dx_mm,dy_mm,dz_mm read from the NIfTI header)
# With REAL data, point --dataset (or $BRATS_DIR) at the folder holding either
# the extracted 'BraTS2024-BraTS-GLI-TrainingData/' tree or the dataset zip.
```

**(c) Run each phase on the contours** — ✅ VERIFIED
```bash
BIN=~/builds/contour/contour_interpolator
C=data/contours/BraTS-GLI-00000-000
$BIN linear $C/slice_0010.obj $C/slice_0011.obj 0.5 /tmp/lin.obj          # M1
$BIN sdf    $C/slice_0010.obj $C/slice_0011.obj 0.5 /tmp/sdf.obj          # M3
$BIN series --kind spline --window $C --upsample 4 /tmp/series_out        # M2
$BIN reconstruct --contours $C --method linear --out /tmp/mesh.off --dz 1 --upsample 3
$BIN metrics --mesh /tmp/mesh.off --contours $C --dz 1 --grid 40
python3 scripts/evaluate.py --contours $C --bin $BIN --out /tmp/loo       # LOO table + plots
```

### Helper scripts (`scripts/`)
| Script | Purpose |
|--------|---------|
| `run_pipeline.sh` | full end-to-end driver (real or synthetic data); also renders the meshes |
| `setup_eigen.sh` | vendor Eigen with no sudo |
| `build_all.sh` | configure + build + ctest |
| `make_synthetic_case.py` | synthetic BraTS-format dataset generator |
| `render_mesh.py` | headless 3D PNG snapshots of a reconstructed `.off` surface |
| `audit_env.sh` | environment / dependency audit |
| `run_all_demo.sh`, `run_reconstruct_demo.sh`, `run_metrics_demo.sh` | per-stage demos |

### 3D output formats & images
Reconstruction writes **`.off`** and **`.ply`** (open directly in MeshLab, Blender,
ParaView, Windows 3D Viewer). `scripts/render_mesh.py` produces **PNG previews**
of the Poisson surface (two shaded views), and `run_pipeline.sh` calls it
automatically (e.g. `data/contours/mesh_<case>_<method>.png`). If you need
another mesh format (`.stl`, `.glb`, `.3ds`), convert from `.ply`/`.off` with
`meshio convert mesh.ply mesh.stl` or via Blender/MeshLab — no pipeline change.

## 5. Known issues / fixes applied

| # | Issue | Resolution |
|---|-------|-----------|
| 1 | Hard-coded macOS dataset path `extract_contours.py:43-44` | ✅ **Fixed:** reads `--dataset` / `$BRATS_DIR`; added `--cases N`; clean error if unset. |
| 2 | Slice spacing not in `index.csv`; volumes assumed 1 mm | ✅ **Fixed:** spacing taken from `nib…header.get_zooms()`; contour coords scaled to mm; `dx_mm,dy_mm,dz_mm` written to `index.csv`. |
| 3 | Python deps `nibabel`/`scikit-image`/`tqdm` missing | ✅ **Installed** (`pip --user --break-system-packages`). |
| 4 | Eigen not system-installed | ✅ Vendored no-sudo (`scripts/setup_eigen.sh`); pass `-DEIGEN3_INCLUDE_DIR`. |
| 5 | Build env is WSL, not the Windows host | ✅ Documented; all build/run via `wsl`. |
| 6 | CGAL **5.6** not 6.x | ✅ Handled in code (`AABB_traits`, hand-rolled cotangent-Laplacian curvature). |
| 7 | PEP-668 externally-managed Python | ✅ `--break-system-packages` or venv. |
| 8 | In-source build forbidden | ✅ Always `-S/-B` with a separate build dir. |
| 9 | No BraTS data on disk | ⏩ **Remaining input only.** `make_synthetic_case.py` substitutes for testing; set `BRATS_DIR` for the real run. |

## 6. Smoke test results (verified this session)

| Step | Result |
|------|--------|
| CMake configure + build (WSL) | ✅ success |
| `ctest` | ✅ **5/5 pass** (linear, series, sdf, reconstruct, metrics) |
| `extract_contours.py` on synthetic dataset | ✅ 102 contours, `index.csv` with spacing |
| reconstruct (linear/spline) | ✅ closed, volume-bounding; stack cross-check −0.5% |
| metrics | ✅ symmetry Dice 0.996, curvature stats, stack cross-check |
| leave-one-slice-out (`evaluate.py`) | ✅ Dice linear 0.985 / spline 0.984 / sdf 0.952 |
| `extract_contours.py` with no dataset | ✅ clean message, exit 1 (no traceback) |

**Bottom line:** everything builds and runs end-to-end today. The single
remaining input is the real BraTS 2024 GLI dataset — provide it via
`BRATS_DIR=/path/to/BraTS bash scripts/run_pipeline.sh`.
