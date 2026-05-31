#!/bin/bash
# End-to-end driver: (build) -> contour extraction -> reconstruction -> metrics
# -> leave-one-slice-out study. Uses real BraTS if $BRATS_DIR is set, otherwise
# generates a synthetic BraTS-format dataset so the whole chain still runs.
#
#   BRATS_DIR=/path/to/BraTS  bash scripts/run_pipeline.sh        # real data
#   bash scripts/run_pipeline.sh                                  # synthetic
set -e

SRC="/mnt/c/Users/LENOVO LOQ/Desktop/geometria_computacional/taller1/computational-geometry-cgal/src/final-project-interpolation"
BUILD="$HOME/builds/contour_all"
BIN="$BUILD/contour_interpolator"
EIGEN_INC="$HOME/opt/eigen3/usr/include/eigen3"
CASES="${CASES:-2}"

echo "==================== build ===================="
if [ ! -x "$BIN" ]; then
  cmake -S "$SRC" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
        -DEIGEN3_INCLUDE_DIR="$EIGEN_INC" >/dev/null
  cmake --build "$BUILD" -j >/dev/null
fi
echo "binary: $BIN"

echo; echo "==================== dataset ===================="
if [ -n "$BRATS_DIR" ]; then
  DATASET="$BRATS_DIR"; echo "using real BraTS: $DATASET"
else
  DATASET="$HOME/synth_brats"
  echo "no BRATS_DIR -> generating synthetic dataset at $DATASET"
  python3 "$SRC/scripts/make_synthetic_case.py" --out "$DATASET" --cases "$CASES" \
          --spacing 1 1 1
fi

echo; echo "==================== extract contours ===================="
python3 "$SRC/scripts/extract_contours.py" --dataset "$DATASET" --cases "$CASES"

CONT="$SRC/data/contours"
IDX="$CONT/index.csv"
CASE="$(ls "$CONT" | grep '^BraTS-GLI-' | head -1)"
CDIR="$CONT/$CASE"
DZ="$(python3 -c "import csv,sys;r=list(csv.DictReader(open(sys.argv[1])));print(r[0].get('dz_mm','1.0') if r else '1.0')" "$IDX")"
echo "first case: $CASE   slices: $(ls "$CDIR"/*.obj 2>/dev/null | wc -l)   dz=$DZ mm (from index.csv)"

echo; echo "==================== reconstruct + metrics ===================="
for METHOD in linear spline; do
  MESH="$CONT/mesh_${CASE}_${METHOD}.off"
  "$BIN" reconstruct --contours "$CDIR" --method "$METHOD" \
         --out "$MESH" --dz "$DZ" --upsample 3
  "$BIN" metrics --mesh "$MESH" --contours "$CDIR" --dz "$DZ" --grid 36
  # 3D preview image of the Poisson surface (headless PNG snapshots)
  python3 "$SRC/scripts/render_mesh.py" "$MESH" \
          "$CONT/mesh_${CASE}_${METHOD}.png" --title "$CASE $METHOD" || true
  echo "---"
done

echo; echo "==================== leave-one-slice-out ===================="
python3 "$SRC/scripts/evaluate.py" --contours "$CDIR" --bin "$BIN" \
        --out "$CONT/loo_${CASE}"
echo "results: $CONT/loo_${CASE}/results.csv  +  comparison.png"

echo; echo "PIPELINE OK"
