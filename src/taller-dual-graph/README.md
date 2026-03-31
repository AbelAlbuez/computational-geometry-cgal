# Taller 2 — Grafo Dual de un Polígono Simple

**Autores:**
- Abel Albuez Sánchez — aa-albuezs@javeriana.edu.co
- Santiago Gil Gallego — santiago_gil@javeriana.edu.co

**Docente:** Leonardo Flórez-Valencia — florez-l@javeriana.edu.co

Maestría en Ingeniería de Sistemas y Computación — Pontificia Universidad Javeriana

---

## Descripción

Implementación del algoritmo de grafo dual de un polígono simple. Dado un polígono de entrada, el programa lo triangula con CGAL y construye el grafo dual: cada triángulo se convierte en un nodo ubicado en su baricentro, y dos nodos se conectan si los triángulos correspondientes comparten una diagonal interior. Los triángulos que tocan el borde del polígono se conectan a un nodo especial P∞ que representa la cara exterior.

---

## Estructura del proyecto

```
taller-dual-graph/
├── CMakeLists.txt              — build: dual_graph + visualizer (CGAL + VTK)
├── requirements.txt            — dependencias Python (solo polygon_drawer.py)
├── output/                     — carpetas auto-numeradas con los resultados
│   ├── result-001/
│   │   ├── 00_polygon.png
│   │   ├── 01_triangulation.png
│   │   ├── 02_barycenters.png
│   │   ├── 03_internal_edges.png
│   │   ├── 04_dual_complete.png
│   │   └── dual_graph.gif      ← solo si se usó --gif
│   └── result-002/ …
├── data/
│   ├── 01_triangulo.obj
│   ├── 02_cuadrado_ccw.obj
│   ├── 03_cuadrado_cw.obj
│   ├── 04_pentagono_regular.obj
│   ├── 05_L.obj
│   ├── 06_T.obj
│   ├── 07_diente.obj
│   ├── 08_flecha.obj
│   ├── 09_escalera.obj
│   ├── poly_00.obj
│   └── poly_01.obj             — polígono de 21 vértices (dibujado con polygon_drawer)
├── lib/
│   └── pujCGAL/
│       ├── Polygon.h / .hxx
│       ├── Triangulation.h / .hxx
│       ├── IO.h / .hxx
│       ├── DualGraph.h / .hxx       ← implementación propia
│       └── IO_DualGraph.h / .hxx   ← implementación propia
└── src/
    ├── triangulate.cxx         — referencia del profesor
    ├── dual_graph.cxx          — ejecutable principal: triangula + construye grafo dual
    ├── visualizer.cxx          — visualizador C++/VTK: PNGs + ventana interactiva + GIF
    └── polygon_drawer.py       — herramienta Python/VTK para dibujar polígonos
```

---

## Requisitos

### C++

- CMake ≥ 3.14
- C++17
- CGAL
- VTK ≥ 9.0
- **ffmpeg** en el PATH (solo para `--gif`)

```bash
# Ubuntu / Debian
sudo apt-get install cmake libcgal-dev libgmp-dev libmpfr-dev libvtk9-dev ffmpeg

# macOS
brew install cmake cgal gmp mpfr vtk ffmpeg
```

### Python (opcional — solo `polygon_drawer.py`)

**Paso 1 — dependencias del sistema:**

```bash
# Ubuntu / Debian
sudo apt-get install libpango1.0-dev libcairo2-dev pkg-config python3-dev python3-venv

# macOS
brew install pango cairo pkg-config
```

**Paso 2 — crear y activar el entorno virtual:**

```bash
cd taller-dual-graph
python3 -m venv .venv
source .venv/bin/activate    # Linux / macOS
# .venv\Scripts\Activate.ps1  # Windows PowerShell
```

**Paso 3 — instalar dependencias:**

```bash
pip install -r requirements.txt
```

**Desactivar cuando termines:**

```bash
deactivate
```

---

## Compilar

```bash
cd taller-dual-graph
mkdir -p build && cd build
cmake ..
make dual_graph visualizer
```

Los ejecutables quedan en `build/`.

---

## Uso

### 1. Dibujar un polígono (opcional)

Si quieres crear un polígono nuevo en lugar de usar los de `data/`:

```bash
source .venv/bin/activate    # activar el venv
cd src
python polygon_drawer.py
```

Controles:
- Clic izquierdo — agregar vértice
- Clic derecho — cerrar polígono
- Clic del medio — deshacer último vértice
- `s` — guardar como `polygon.obj`
- `c` — limpiar y empezar de nuevo
- `q` / ESC — salir

### 2. Calcular el grafo dual

```bash
./build/dual_graph  <input.obj>  <triangulation.obj>  <dual.obj>
```

**Ejemplo:**

```bash
./build/dual_graph  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj
```

El programa imprime en consola el número de triángulos, aristas y la matriz de adyacencia.

Salidas:
- `<triangulation.obj>` — triangulación CGAL del polígono
- `<dual.obj>` — grafo dual: baricentros + P∞ + aristas (formato `v` y `l`)

### 3. Visualizar y generar imágenes

El visualizador tiene tres modos:

#### Modo por defecto — genera 5 PNGs

```bash
./build/visualizer  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj
```

Crea automáticamente una carpeta nueva `output/result-NNN/` con 5 PNGs que
muestran el proceso de construcción paso a paso:

```
output/result-001/
├── 00_polygon.png          — polígono de entrada
├── 01_triangulation.png    — + triangulación
├── 02_barycenters.png      — + baricentros (nodos del dual)
├── 03_internal_edges.png   — + aristas internas
└── 04_dual_complete.png    — grafo dual completo
```

Cada ejecución crea una carpeta nueva (`result-002`, `result-003`, …) sin
pisar los resultados anteriores.

#### Modo `--gif` — PNGs + animación

```bash
./build/visualizer  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj  --gif
```

Genera los mismos 5 PNGs más `dual_graph.gif` en la misma carpeta.
Requiere `ffmpeg` instalado en el PATH.

#### Modo `--interactive` — ventana VTK interactiva

```bash
./build/visualizer  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj  --interactive
```

Abre una ventana con las tres capas. Controles:
- `1` — mostrar/ocultar polígono
- `2` — mostrar/ocultar triangulación
- `3` — mostrar/ocultar grafo dual
- `s` — captura de pantalla → `screenshot_dual.png` (directorio actual)
- `q` / ESC — salir

---

## Formato de los archivos `.obj`

Los archivos de entrada usan el mismo subconjunto OBJ del resto del proyecto:

```
# comentario
v 0.0 0.0 0.0    ← vértice 2D (x y z, z siempre 0)
v 4.0 0.0 0.0
v 4.0 4.0 0.0
v 0.0 4.0 0.0
f 1 2 3 4        ← cara (índices base 1)
```

Los vértices deben estar en orden **CCW (antihorario)**. Si están en CW,
`dual_graph` los invierte automáticamente con `guarantee_CCW()`.

El archivo del grafo dual usa `v` para los baricentros (P∞ como último
vértice) y `l` para las aristas:

```
# Dual graph: 4 nodes, 2 internal edges, 3 external edges
v 3.33 2.00 0
v 2.00 0.67 0
v 0.67 2.00 0
v 10.0 -6.00 0    ← P∞ (último vértice)
l 1 2
l 2 3
l 1 4
```

---

## Flujo completo de una sola vez

```bash
# desde la raíz del taller
mkdir -p build && cd build && cmake .. && make dual_graph visualizer && cd ..

# calcular el grafo dual
./build/dual_graph  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj

# generar PNGs (modo por defecto)
./build/visualizer  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj

# generar PNGs + GIF
./build/visualizer  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj  --gif
```

---

## Solución de problemas

**CGAL no encontrado:**
```bash
sudo apt-get install libcgal-dev
# o especificar la ruta:
cmake .. -DCGAL_DIR=/ruta/a/cgal
```

**VTK no encontrado:**
```bash
# Ubuntu
sudo apt-get install libvtk9-dev
# macOS
brew install vtk
```

**`--gif` falla o no genera el GIF:**
Instala `ffmpeg` y asegúrate de que esté en el PATH:
```bash
which ffmpeg    # debe mostrar una ruta
# si no: brew install ffmpeg  /  sudo apt-get install ffmpeg
```
Si ffmpeg falla, los PNGs quedan en `output/result-NNN/` y puedes unirlos manualmente.

**Error de in-source build:**
No ejecutar `cmake` desde la raíz del taller. Usar siempre `build/` separado.

**`polygon_drawer` no arranca (ModuleNotFoundError):**
El venv no está activo. Ejecutar `source .venv/bin/activate` antes de correr cualquier script Python.
