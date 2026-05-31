#!/bin/bash
# Environment audit for building/running the final project on this machine.
echo "==================== OS ===================="
uname -a
. /etc/os-release 2>/dev/null && echo "distro: $PRETTY_NAME"
echo "WSL_DISTRO_NAME=${WSL_DISTRO_NAME:-<none>}"

echo; echo "==================== CGAL ===================="
if [ -f /usr/include/CGAL/version.h ]; then
  grep "define CGAL_VERSION " /usr/include/CGAL/version.h
else
  echo "no /usr/include/CGAL/version.h"
fi
dpkg -l libcgal-dev 2>/dev/null | grep '^ii' || echo "libcgal-dev: not via dpkg"
echo -n "pkg-config cgal: "; pkg-config --modversion cgal 2>/dev/null || echo "n/a (CGAL is header-only)"

echo; echo "==================== Eigen ===================="
test -f /usr/include/eigen3/Eigen/Dense && echo "system eigen3: present" || echo "system eigen3: ABSENT"
test -f "$HOME/opt/eigen3/usr/include/eigen3/Eigen/Dense" \
  && echo "vendored eigen3: $HOME/opt/eigen3/usr/include/eigen3 (present)" \
  || echo "vendored eigen3: ABSENT"

echo; echo "==================== Compiler / CMake ===================="
g++ --version | head -1
cmake --version | head -1
echo -n "make: "; make --version 2>/dev/null | head -1

echo; echo "==================== Python ===================="
python3 --version
echo "--- requirements.txt deps ---"
for mod in nibabel numpy matplotlib skimage tqdm; do
  v=$(python3 -c "import $mod,sys; print(getattr($mod,'__version__','?'))" 2>/dev/null)
  if [ -n "$v" ]; then printf "  %-12s OK    %s\n" "$mod" "$v"
  else printf "  %-12s MISSING\n" "$mod"; fi
done

echo; echo "==================== BraTS data ===================="
RAW="/mnt/c/Users/LENOVO LOQ/Desktop/geometria_computacional/taller1/computational-geometry-cgal/src/final-project-interpolation/data/raw"
if [ -d "$RAW" ]; then echo "data/raw: present ($(ls "$RAW" | wc -l) entries)"; else echo "data/raw: ABSENT (no dataset extracted)"; fi
echo "DONE"
