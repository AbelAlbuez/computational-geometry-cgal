# computational-geometry-cgal

**Autor:** Abel Albuez Sanchez

---

## Descripción

Proyecto académico de geometría computacional que implementa **intersección de segmentos y polígonos convexos** en C++ con CGAL. Compara una implementación manual (Sutherland-Hodgman) con la solución nativa de CGAL.

---

## Estructura del proyecto

```
computational-geometry-cgal/
├── src/
│   ├── no_cgal/
│   │   └── brute_force.cxx              — intersección sin CGAL
│   ├── polygon_intersection/
│   │   ├── compare.cxx                  — compara BruteForce vs BentleyOttmann
│   │   ├── brute_force.cxx              — fuerza bruta standalone
│   │   ├── brute_force_2.cxx            — fuerza bruta usando pujCGAL/IO
│   │   ├── main.cxx                     — [TALLER] Sutherland-Hodgman manual
│   │   └── main_cgal.cxx                — [TALLER] intersección con CGAL nativo
│   └── cgal/
│       ├── compare.cxx
│       ├── brute_force.cxx
│       └── brute_force_2.cxx
├── lib/
│   └── pujCGAL/
│       ├── SegmentsIntersection.h / .hxx
│       └── IO.h / .hxx
├── data/
│   ├── input_00.obj .. input_03.obj     — segmentos de prueba del profesor
│   ├── poly_P.obj                       — [TALLER] polígono P de prueba
│   └── poly_Q.obj                       — [TALLER] polígono Q de prueba
├── CMakeLists.txt
└── baseConf.cmake
```

---

## Requisitos

- **C++23** (configurado en `baseConf.cmake`)
- **CMake** ≥ 3.6
- **CGAL** con componente Core
- **Boost** (requerido por CGAL)
- **GMP y MPFR** (requeridos por el kernel exacto de CGAL)

### Instalación de dependencias

**Ubuntu / Debian:**
```bash
sudo apt-get install cmake libcgal-dev libgmp-dev libmpfr-dev libboost-all-dev
```

**macOS:**
```bash
brew install cmake cgal boost gmp mpfr
```

---

## Compilar

```bash
# Desde la raíz del proyecto
mkdir build
cd build
cmake ..
make
```

Los ejecutables quedan en `build/`.

Para compilar solo un target específico:
```bash
make polygon_intersection
make polygon_intersection_cgal
```

---

## Ejecutables del taller

### `polygon_intersection` — Sutherland-Hodgman manual

Calcula la intersección de dos polígonos convexos usando el algoritmo de recorte por semiplanos.

```bash
./polygon_intersection <poligono_P.obj> <poligono_Q.obj> <salida.obj>
```

**Ejemplo:**
```bash
./polygon_intersection ../data/poly_P.obj ../data/poly_Q.obj ../data/result.obj
```

---

### `polygon_intersection_cgal` — CGAL nativo

Misma operación usando `CGAL::Polygon_2` y `CGAL::intersection` con kernel de aritmética exacta.

```bash
./polygon_intersection_cgal <poligono_P.obj> <poligono_Q.obj> <salida.obj>
```

**Ejemplo:**
```bash
./polygon_intersection_cgal ../data/poly_P.obj ../data/poly_Q.obj ../data/result_cgal.obj
```

---

## Formato de los archivos `.obj`

Los archivos de entrada y salida usan un subconjunto del formato OBJ:

```
# comentario
v 0 0       ← vértice 2D (x y)
v 2 0
v 2 2
v 0 2
l 1 2       ← arista: índice base-1
l 2 3
l 3 4
l 4 1
```

Los vértices deben estar en orden **CCW (antihorario)** para que el predicado de semiplano funcione correctamente.

---

## Casos manejados

| Caso | Resultado esperado |
|---|---|
| Intersección parcial | Polígono convexo recortado |
| Sin intersección | 0 vértices (archivo vacío) |
| P contenido en Q | Devuelve P completo |
| Q contenido en P | Devuelve Q completo |

---

## Otros ejecutables del proyecto base

```bash
# Fuerza bruta sin CGAL
./no_cgal_brute_force <input.obj> <output.obj>

# Fuerza bruta con CGAL
./cgal_brute_force <input.obj> <output.obj>

# Comparar BruteForce vs BentleyOttmann
./cgal_compare <input.obj> <bf_output.obj> <bo_output.obj>
```

---

## Solución de problemas

**CGAL no encontrado:**
```bash
# Ubuntu
sudo apt-get install libcgal-dev
# o pasar la ruta manualmente
cmake .. -DCGAL_DIR=/ruta/a/cgal
```

**Error de GMP / MPFR al compilar `polygon_intersection_cgal`:**
```bash
sudo apt-get install libgmp-dev libmpfr-dev
```

**Error de in-source build:**
No ejecutar `cmake` desde la raíz del proyecto. Siempre usar un directorio `build/` separado.

---

## Licencia

Ver el archivo [LICENSE](LICENSE) en la raíz del proyecto.