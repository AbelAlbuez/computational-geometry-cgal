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

Lee una imagen PNG como heightmap, construye una triangulación de Delaunay 2D
sobre los píxeles y simplifica la malla eliminando vértices redundantes en zonas
planas, conservando los que definen el relieve del terreno (bordes, crestas, valles).

### Algoritmo

Para cada vértice `p` se calcula un **error de planitud**:

```
error(p) = promedio de |p.z - qi.z| para cada vecino qi
```

Se usa un **min-heap (Dijkstra)** para procesar primero los vértices de menor error.
Si `error(p) < ε`, el vértice se elimina con `T.remove(p)` y CGAL re-triangula el
hueco automáticamente (Delaunay local). Los errores de los vecinos en la vecindad
de orden k se recalculan y reinsertan en el heap (**lazy update**).

### Conceptos del curso aplicados

| Concepto | Aplicación |
|---|---|
| **DCEL / half-edge** (Clase 6) | `Face_circulator` recorre la estrella de p con twin→next |
| **Triangulación de Delaunay** (Clase 5) | Malla base y re-triangulación automática tras cada remoción |
| **Dualidad Voronoi-Delaunay** (Clase 5) | Eliminar un vértice fusiona celdas de Voronoi vecinas |
| **Face-vertex mesh** (Clase 6) | Estructura interna de la triangulación CGAL |
| **Min-heap (Dijkstra)** | Orden de procesamiento: menor error de planitud sale primero |

### Estructura del taller

```
src/taller-3-codigo-base/
├── src/
│   └── heightmap.cxx          ← algoritmo principal
├── lib/
│   └── pujCGAL/               ← librería del profesor
├── data/                      ← imágenes PNG de entrada
├── output/                    ← mallas .obj + visualizaciones por imagen
│   └── test-[nombre]/
│       ├── original.obj
│       ├── simplificado.obj
│       ├── paso_01.png .. paso_07.png
│       ├── resultado.gif
│       └── paso_a_paso.gif
├── diagramas/                 ← diagramas explicativos del algoritmo
│   ├── 01_criterio_planitud.png
│   ├── 02_face_circulator.png
│   ├── 03_min_heap.png
│   ├── 04_remove_retrigulacion.png
│   ├── 05_lazy_deletion.png
│   ├── 06_vecindad_orden_k.png
│   ├── 07_resultados.png
│   └── algoritmo_completo.gif
├── generar_diagramas.py       ← script que genera los diagramas
├── visualizer.py              ← genera GIFs paso a paso por imagen
├── requirements.txt           ← dependencias Python
└── CMakeLists.txt
```

### Compilar

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
# Desde build/
./taller3 <input.png> <epsilon> <orden_k>
```

| Parámetro | Descripción | Rango útil |
|---|---|---|
| `input.png` | Imagen PNG de entrada (heightmap) | — |
| `epsilon` | Umbral de planitud normalizado a [0,1] | 0.01 – 0.15 |
| `orden_k` | Radio de vecindad para actualización | 1 – 3 |

**Ejemplo:**
```bash
./taller3 "../../../data/SRTM_US_scaled_512.png" 0.05 2
```

### Resultados (ε = 0.05, orden_k = 2)

| Imagen | Vértices antes | Vértices después | Reducción |
|---|---:|---:|---:|
| SRTM_US_scaled_256.png | 32,768 | 4,840 | **85.2%** |
| SRTM_US_scaled_512.png | 131,072 | 14,190 | **89.2%** |
| NormalMap-second.png | 196,608 | 22,352 | **88.6%** |
| NormalMap.png | 262,144 | 1,647 | **99.4%** |
| SRTM_US_scaled_1024.png | 524,288 | 37,188 | **92.9%** |
| SRTM_US_scaled_2048.png | 2,097,152 | 81,758 | **96.1%** |

### Generar diagramas explicativos

```bash
cd src/taller-3-codigo-base

# Crear entorno virtual
python3 -m venv .venv
source .venv/bin/activate

# Instalar dependencias
pip install -r requirements.txt

# Generar los 7 PNGs + GIF animado
python3 generar_diagramas.py

deactivate
```

Los archivos quedan en `src/taller-3-codigo-base/diagramas/`.

### Visualizador paso a paso

```bash
# Se ejecuta automáticamente después de cada taller3
# También puede ejecutarse manualmente:
python3 visualizer.py --output output/test-[nombre] --name [nombre]
```

Genera un GIF de 7 pasos que muestra:
1. Malla original densa
2. Heatmap de error de planitud (verde=plano, rojo=borde)
3. Vértice más plano resaltado (Dijkstra)
4. Estrella con Face_circulator
5. Expansión a orden k
6. Estado intermedio de eliminación
7. Resultado final simplificado

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