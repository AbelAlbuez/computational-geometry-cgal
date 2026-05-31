#!/bin/bash
# Verbose reconstruct test, then a full CLI run on a generated contour stack
# (a "tumor-like" blobby ellipsoid), for the linear method + Boissonnat band.
set -e
BUILD="$HOME/builds/contour_all"
BIN="$BUILD/contour_interpolator"
TMP="$(mktemp -d)"
mkdir -p "$TMP/case"

echo "############ test_reconstruct (verbose) ############"
"$BUILD/test_reconstruct"

# Generate a stack: a lobed ellipsoid sampled as axial contours slice_0000..NNNN
python3 - "$TMP/case" <<'PY'
import math,sys,os
d=sys.argv[1]; a=22.0; c=26.0
for z in range(2, int(2*c)-1):
    u=(z-c)/c
    base=a*math.sqrt(max(0.0,1-u*u))
    if base<3: continue
    m=max(24,int(2*math.pi*base/2))
    with open(os.path.join(d,f"slice_{z:04d}.obj"),'w') as f:
        f.write(f"# slice {z}\n")
        pts=[]
        for i in range(m):
            ang=2*math.pi*i/m
            # lobed radius (non-trivial shape) so it is not a trivial sphere
            r=base*(1.0+0.12*math.cos(3*ang))
            pts.append((r*math.cos(ang), r*math.sin(ang)))
        f.write(f"# {len(pts)} vertices\n")
        for x,y in pts: f.write(f"v {x:.6f} {y:.6f}\n")
        for i in range(1,len(pts)): f.write(f"l {i} {i+1}\n")
        f.write(f"l {len(pts)} 1\n")
print("generated", len(os.listdir(d)), "slices")
PY

echo; echo "############ CLI reconstruct (linear, +boissonnat) ############"
"$BIN" reconstruct --contours "$TMP/case" --method linear \
       --out "$TMP/mesh_linear.off" --dz 1 --upsample 3 \
       --boissonnat "$TMP/band.off"

echo; echo "############ CLI reconstruct (spline) ############"
"$BIN" reconstruct --contours "$TMP/case" --method spline \
       --out "$TMP/mesh_spline.off" --dz 1 --upsample 3

echo; echo "############ output files ############"
ls -la "$TMP"/*.off "$TMP"/*.ply 2>/dev/null | awk '{print $5, $9}'

rm -rf "$TMP"
echo; echo "RECONSTRUCT DEMO OK"
