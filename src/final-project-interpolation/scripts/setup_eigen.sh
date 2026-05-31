#!/bin/bash
# Vendors Eigen (header-only) into WSL home WITHOUT sudo, by downloading the
# Ubuntu .deb as a normal user and unpacking it locally. CGAL's build is then
# pointed at it via Eigen3_DIR. This does not touch the system CGAL install.
set -e

DEST="$HOME/opt/eigen3"
WORK="$HOME/opt/eigen_pkg"

echo "=== preparing $WORK ==="
rm -rf "$WORK" "$DEST"
mkdir -p "$WORK" "$DEST"
cd "$WORK"

echo "=== downloading libeigen3-dev .deb as normal user ==="
apt-get download libeigen3-dev

DEB=$(ls *.deb | head -1)
echo "got: $DEB"

echo "=== unpacking .deb into $DEST ==="
dpkg -x "$DEB" "$DEST"

echo "=== locating headers and cmake config ==="
HDR=$(find "$DEST" -type d -name Eigen | head -1)
CFG=$(find "$DEST" -name Eigen3Config.cmake | head -1)
echo "EIGEN_HEADERS_DIR=$(dirname "$HDR")"
echo "EIGEN3_CONFIG=$CFG"
echo "EIGEN3_CONFIG_DIR=$(dirname "$CFG")"

echo "=== sanity: Dense present? ==="
test -f "$(dirname "$HDR")/Eigen/Dense" && echo "EIGEN_OK" || echo "EIGEN_MISSING"
