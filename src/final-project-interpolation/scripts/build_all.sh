#!/bin/bash
# Configure + build every target (interpolator + 3 tests) in a WSL-local build
# dir, then run the full ctest suite. Interpolation core needs no Eigen.
set -e
SRC="/mnt/c/Users/LENOVO LOQ/Desktop/geometria_computacional/taller1/computational-geometry-cgal/src/final-project-interpolation"
BUILD="$HOME/builds/contour_all"
EIGEN_INC="$HOME/opt/eigen3/usr/include/eigen3"   # vendored, no sudo (setup_eigen.sh)

echo "=== configure ==="
cmake -S "$SRC" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DEIGEN3_INCLUDE_DIR="$EIGEN_INC" 2>&1 | tail -6

echo "=== build ==="
cmake --build "$BUILD" -j 2>&1 | tail -20

echo "=== ctest ==="
( cd "$BUILD" && ctest --output-on-failure )
