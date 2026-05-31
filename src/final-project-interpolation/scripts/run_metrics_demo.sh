#!/bin/bash
# Phase 5 demo: verbose metrics test, a `metrics` CLI run on a reconstructed
# case (with the stack cross-check), and the leave-one-slice-out study.
set -e
BUILD="$HOME/builds/contour_all"
BIN="$BUILD/contour_interpolator"
HERE="$(cd "$(dirname "$0")" && pwd)"
TMP="$(mktemp -d)"
mkdir -p "$TMP/case"

echo "############ test_metrics (verbose) ############"
"$BUILD/test_metrics"

# Generate a lobed-ellipsoid contour stack.
python3 - "$TMP/case" <<'PY'
import math,sys,os
d=sys.argv[1]; a=22.0; c=26.0
for z in range(2, int(2*c)-1):
    u=(z-c)/c
    base=a*math.sqrt(max(0.0,1-u*u))
    if base<3: continue
    m=max(24,int(2*math.pi*base/2))
    with open(os.path.join(d,f"slice_{z:04d}.obj"),'w') as f:
        f.write(f"# slice {z}\n# n\n")
        pts=[(base*(1+0.12*math.cos(3*2*math.pi*i/m))*math.cos(2*math.pi*i/m),
              base*(1+0.12*math.cos(3*2*math.pi*i/m))*math.sin(2*math.pi*i/m)) for i in range(m)]
        for x,y in pts: f.write(f"v {x:.6f} {y:.6f}\n")
        for i in range(1,len(pts)): f.write(f"l {i} {i+1}\n")
        f.write(f"l {len(pts)} 1\n")
print("generated", len(os.listdir(d)), "slices")
PY

echo; echo "############ reconstruct + metrics CLI ############"
"$BIN" reconstruct --contours "$TMP/case" --method linear --out "$TMP/mesh.off" --dz 1 --upsample 3
echo "---"
"$BIN" metrics --mesh "$TMP/mesh.off" --contours "$TMP/case" --dz 1 --grid 36

echo; echo "############ leave-one-slice-out study (evaluate.py) ############"
if python3 -c "import numpy, matplotlib" 2>/dev/null; then
  python3 "$HERE/evaluate.py" --contours "$TMP/case" --bin "$BIN" --out "$TMP/loo"
  echo "--- results.csv ---"; cat "$TMP/loo/results.csv"
  ls -la "$TMP/loo/comparison.png" | awk '{print "plot:", $5, "bytes", $9}'
else
  echo "numpy/matplotlib not available in WSL python3 -> skipping evaluate.py"
  echo "(install with: python3 -m pip install --user -r requirements.txt)"
fi

rm -rf "$TMP"
echo; echo "METRICS DEMO DONE"
