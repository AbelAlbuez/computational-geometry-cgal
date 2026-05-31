#!/bin/bash
# Show the numeric output of all three test binaries, then exercise the CLI for
# every interpolation method (linear / series-spline / series-poly / sdf).
set -e
BUILD="$HOME/builds/contour_all"
BIN="$BUILD/contour_interpolator"
TMP="$(mktemp -d)"

echo "############ unit tests (verbose) ############"
echo "----- test_linear -----";  "$BUILD/test_linear"
echo "----- test_series -----";  "$BUILD/test_series"
echo "----- test_sdf -----";     "$BUILD/test_sdf"

gen() { # cx cy r m path
  python3 - "$@" <<'PY'
import math,sys
cx,cy,r,m,path=float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]),sys.argv[5]
with open(path,'w') as f:
    f.write(f"# circle\n# {m} vertices\n")
    for i in range(m):
        a=2*math.pi*i/m
        f.write(f"v {cx+r*math.cos(a):.6f} {cy+r*math.sin(a):.6f}\n")
    for i in range(1,m): f.write(f"l {i} {i+1}\n")
    f.write(f"l {m} 1\n")
PY
}

echo; echo "############ CLI: linear (M1) ############"
gen 0 0 10 64 "$TMP/A.obj"
gen 6 2 18 80 "$TMP/B.obj"
"$BIN" linear "$TMP/A.obj" "$TMP/B.obj" 0.5 "$TMP/lin.obj"

echo; echo "############ CLI: sdf (M3) ############"
"$BIN" sdf "$TMP/A.obj" "$TMP/B.obj" 0.5 "$TMP/sdf.obj"
"$BIN" sdf "$TMP/A.obj" "$TMP/B.obj" 0.5 "$TMP/sdf_na.obj" --no-align

echo; echo "############ CLI: series (M2) ############"
# Build a 5-slice stack of a growing, drifting circle: slice_000..slice_004
mkdir -p "$TMP/stack"
for z in 0 1 2 3 4; do
  cxv=$(python3 -c "import math;print(8*math.sin(0.5*$z))")
  rv=$(python3 -c "print(15+2*$z)")
  gen "$cxv" 0 "$rv" 50 "$TMP/stack/slice_000$z.obj"
done
echo "--- spline, upsample 4 ---"
"$BIN" series --kind spline --window "$TMP/stack" --upsample 4 "$TMP/out_spline"
echo "--- polynomial, upsample 4 ---"
"$BIN" series --kind poly --window "$TMP/stack" --upsample 4 "$TMP/out_poly"
echo "produced files:"
ls "$TMP/out_spline" | head -3; echo "  ... ($(ls "$TMP/out_spline" | wc -l) total)"

rm -rf "$TMP"
echo; echo "ALL DEMOS OK"
