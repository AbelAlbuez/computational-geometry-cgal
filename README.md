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

## Taller 3 — Simplificación de Heightmap por Vecindad de Orden k

### Descripción

Lee una imagen PNG como heightmap, construye una triangulación de Delaunay 2D sobre los píxeles y **simplifica** la malla eliminando vértices redundantes (zonas planas), conservando los vértices importantes (bordes, crestas, valles).

**Algoritmo:** para cada vértice `p` se calcula un error de planitud = promedio de `|h(p) - h(q)|` sobre sus vecinos `q`. Se usa un min-heap (estilo Dijkstra) para procesar primero los vértices con menor error. Si `error(p) < ε`, el vértice es eliminado con `T.remove(p)` y CGAL re-triangula el hueco automáticamente (Delaunay local). Los errores de los vecinos de `p` en la vecindad de orden k se recalculan y reinsertan en el heap (lazy update).

**Nota:** los valores de altura quedan normalizados a `[0, 1]` por la clase `pujCGAL::Heightmap`, por lo que `epsilon` debe estar en ese rango. Un valor típico es `0.02`–`0.10`.

### Conceptos del curso aplicados

| Concepto | Aplicación |
|---|---|
| **DCEL / half-edge** (Clase 6) | `Face_circulator` de CGAL recorre la estrella de p en orden circular (twin→next) |
| **Triangulación de Delaunay** (Clase 5) | La malla de partida y la re-triangulación automática tras cada remoción |
| **Min-heap (Dijkstra)** | Orden de procesamiento: menor error de planitud sale primero |
| **Dualidad Voronoi-Delaunay** (Clase 5) | Eliminar un vértice equivale a fusionar celdas de Voronoi vecinas |
| **Face-vertex mesh** (Clase 6) | La triangulación CGAL es exactamente esta estructura |

### Compilar (Taller 3)

```bash
cd src/taller-3-codigo-base
mkdir -p build && cd build
cmake ..
make
```

Requiere: `cmake`, `CGAL`, `libpng`.

```bash
# macOS
brew install cmake cgal libpng

# Ubuntu / Debian
sudo apt-get install cmake libcgal-dev libpng-dev
```

### Ejecutar

```bash
./taller3 input.png output.obj [epsilon] [orden_k]
```

| Parámetro | Descripción | Valor por defecto |
|---|---|---|
| `input.png` | Imagen PNG de entrada (heightmap) | — |
| `output.obj` | Malla simplificada en formato OBJ 3D | — |
| `epsilon` | Umbral de planitud en `[0, 1]` | `10.0` |
| `orden_k` | Radio de vecindad para actualización | `2` |

**Ejemplo:**
```bash
./taller3 ../../../data/input_00.png output_simplificado.obj 0.05 2
```

### Resultados — `data/input_00.png` (500×333 px, RGB)

| epsilon | orden\_k | Vértices antes | Vértices después | Reducción |
|---|---|---|---|---|
| 0.02 | 2 | 166 500 | 14 727 | 91.2 % |
| 0.05 | 2 | 166 500 | 13 831 | 91.7 % |
| 10.0 | 2 | 166 500 | 4 | ~100 % (degenera) |

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