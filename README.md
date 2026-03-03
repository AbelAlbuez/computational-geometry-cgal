# computational-geometry-cgal

**Author:** Abel Albuez Sanchez

---

## Description

This academic project focuses on **geometric intersection** of line segments and polygons. It compares a **manual implementation** (no external geometry library) with implementations using **CGAL** (Computational Geometry Algorithms Library), and extends the base code provided by the course instructor.

Main goals:

- **Segment intersection:** compute intersections between line segments from input data.
- **Comparison:** contrast manual (brute-force) algorithms with CGAL-based solutions.
- **Extensibility:** serve as a base for further work (e.g. convex polygon intersection in the workshop).

---

## Project Structure

```
computational-geometry-cgal/
├── src/
│   ├── no_cgal/
│   │   └── brute_force.cxx
│   ├── polygon_intersection/
│   │   ├── compare.cxx
│   │   ├── brute_force.cxx
│   │   └── brute_force_2.cxx
│   └── cgal/
│       ├── compare.cxx
│       ├── brute_force.cxx
│       └── brute_force_2.cxx
├── lib/
│   └── pujCGAL/
│       ├── SegmentsIntersection.h
│       ├── SegmentsIntersection.hxx
│       ├── IO.h
│       └── IO.hxx
├── data/
│   ├── input_00.obj
│   ├── input_01.obj
│   ├── input_02.obj
│   ├── input_03.obj
│   └── *_intersections.obj
├── CMakeLists.txt
├── baseConf.cmake
├── LICENSE
└── README.md
```

- **`src/no_cgal/`** — implementations without CGAL (manual segment intersection).
- **`src/polygon_intersection/`** — polygon-related intersection code (compare and brute-force variants).
- **`src/cgal/`** — CGAL-based segment intersection (compare and brute-force variants).
- **`lib/pujCGAL/`** — shared headers for segment intersection and I/O (`.obj` read/write).
- **`data/`** — sample inputs (`input_XX.obj`) and optional intersection outputs (`*_intersections.obj`).

---

## Requirements

- **C++17** or later (project is configured for C++23 in `baseConf.cmake`).
- **CMake** ≥ 3.6.
- **CGAL** (with Core component).
- **Boost** (required by CGAL).

---

## Build Instructions

From the project root:

```bash
cd computational-geometry-cgal
mkdir build
cd build
cmake ..
make
```

Generated executables are placed in the **`build`** directory (same as the binary tree). Typical names include targets such as `no_cgal_brute_force`, `cgal_brute_force`, `cgal_brute_force_2`, and `cgal_compare`, depending on the CMake configuration.

---

## Running the Executables

Run any binary from the `build` directory:

```bash
cd build
./<executable_name>
```

To use the sample data as input, pass the path to an `.obj` file. For example, if running from `build`:

```bash
./<executable_name> ../data/input_00.obj
```

Input `.obj` files use a simple format: lines starting with `v` define 2D vertices, and lines starting with `l` define segments by vertex indices. Output intersection results may be written as `*_intersections.obj` in the same format.

---

## Troubleshooting

### CGAL not found

- Install CGAL and its dependencies (Boost, GMP, MPFR). On macOS: `brew install cgal`. On Ubuntu/Debian: `sudo apt-get install libcgal-dev`.
- Ensure CMake can find CGAL (e.g. set `CGAL_DIR` if you use a non-standard installation).

### Boost linking issues

- CGAL depends on Boost. Install Boost (e.g. `brew install boost` or `sudo apt-get install libboost-all-dev`) and ensure it is visible to CMake.
- If linking fails, check that the same compiler and standard library are used for Boost and CGAL.

### CMake errors

- Do not run CMake from the project root as the build directory; use a separate `build` folder (e.g. `mkdir build && cd build && cmake ..`). The project explicitly disallows in-source builds.
- If the C++ standard or compiler is wrong, adjust `baseConf.cmake` or pass `-DCMAKE_CXX_STANDARD=17` (or higher) when configuring.

---

## Workshop Note

The workshop will extend this codebase to implement **intersection of convex polygons** using CGAL. The current structure (segment intersection, I/O, and CGAL usage) is intended as a base for that task.

---

## License

See the [LICENSE](LICENSE) file in the project root.
