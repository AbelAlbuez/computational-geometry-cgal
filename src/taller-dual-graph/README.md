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
├── CMakeLists.txt          — build: dual_graph + visualizer (CGAL + VTK)
├── requirements.txt        — dependencias Python (solo polygon_drawer)
├── data/
│   ├── 01_triangulo.obj
│   ├── 02_cuadrado_ccw.obj
│   ├── 03_cuadrado_cw.obj
│   ├── …
│   ├── poly_00.obj
│   └── poly_01.obj         — polígono de 21 vértices (dibujado con polygon_drawer)
├── lib/
│   └── pujCGAL/
│       ├── Polygon.h / .hxx
│       ├── Triangulation.h / .hxx
│       ├── IO.h / .hxx
│       ├── DualGraph.h / .hxx
│       └── IO_DualGraph.h / .hxx
└── src/
    ├── triangulate.cxx     — referencia del profesor
    ├── dual_graph.cxx      — triangulación + grafo dual → .obj
    ├── visualizer.cxx      — VTK: ventana 3 capas + GIF (--gif) vía ffmpeg
    └── polygon_drawer.py   — herramienta interactiva para dibujar polígonos (Python + VTK)
```

---

## Requisitos

### C++

- CMake ≥ 3.14
- C++17
- CGAL
- VTK ≥ 9.0 (visualizador)
- **ffmpeg** en el PATH (solo para `./visualizer … --gif`)

```bash
# Ubuntu / Debian
sudo apt-get install cmake libcgal-dev libgmp-dev libmpfr-dev libvtk9-dev ffmpeg

# macOS
brew install cmake cgal gmp mpfr vtk ffmpeg
```

### Python (opcional — solo `polygon_drawer.py`)

```bash
cd taller-dual-graph
python3 -m venv .venv
source .venv/bin/activate   # Linux / macOS
pip install -r requirements.txt
```

---

## Compilar

```bash
cd taller-dual-graph
mkdir -p build && cd build
cmake ..
make dual_graph visualizer
```

Ejecutables: `build/dual_graph`, `build/visualizer`.

---

## Uso

### 1. Dibujar un polígono (opcional)

```bash
source .venv/bin/activate   # si usas venv
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

**Salida:**
- `<triangulation.obj>` — triangulación CGAL del polígono
- `<dual.obj>` — grafo dual: baricentros + P∞ + aristas (`v` y `l`)

### 3. Visualizar (VTK, ejecutable C++)

Desde la raíz del taller (o ajusta rutas):

```bash
./build/visualizer  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj
```

Controles:
- `1` — mostrar/ocultar polígono
- `2` — mostrar/ocultar triangulación
- `3` — mostrar/ocultar grafo dual
- `s` — captura → `screenshot_dual.png` (directorio de trabajo actual)
- `q` / ESC — salir

### 4. GIF de la secuencia (off-screen + ffmpeg)

Genera unos cuantos fotogramas PNG en memoria de disco, llama a **ffmpeg** y escribe `dual_graph.gif` (y borra los PNG temporales si ffmpeg tuvo éxito).

```bash
./build/visualizer  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj  --gif
```

Si `ffmpeg` no está instalado, quedan los archivos `dual_viz_frame_*.png` para unirlos a mano.

---

## Formato de los archivos `.obj`

Los archivos de entrada usan el mismo subconjunto OBJ del resto del proyecto:

```
# comentario
v 0.0 0.0 0.0    ← vértice 2D (x y z, z suele ser 0)
v 4.0 0.0 0.0
v 4.0 4.0 0.0
v 0.0 4.0 0.0
f 1 2 3 4        ← cara (índices base 1)
```

Los vértices deben estar en orden **CCW (antihorario)** cuando aplica. Si están en CW, `dual_graph` puede corregir con `guarantee_CCW()`.

El archivo del grafo dual usa `v` para los baricentros (más P∞ como último vértice) y `l` para las aristas:

```
# Dual graph: …
v …
l 1 2
…
```

---

## Flujo completo de una sola vez

```bash
cd taller-dual-graph
mkdir -p build && cd build && cmake .. && make dual_graph visualizer && cd ..
./build/dual_graph  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj
./build/visualizer  data/poly_01.obj  build/out_triang.obj  build/out_dual.obj
```

Para dibujar polígonos nuevos con Python, activa el venv e instala `requirements.txt` antes de `python polygon_drawer.py`.

---

## Solución de problemas

**CGAL no encontrado:**
```bash
sudo apt-get install libcgal-dev
# o: cmake .. -DCGAL_DIR=/ruta/a/cgal
```

**VTK no encontrado (CMake):**
```bash
brew install vtk
# Ubuntu: sudo apt-get install libvtk9-dev
```

**`--gif` no genera `dual_graph.gif`:**
Instala `ffmpeg` y vuelve a ejecutar. Revisa mensajes en consola.

**Error de in-source build:**
Usar siempre un directorio `build/` separado.

**`ModuleNotFoundError: No module named 'vtk'` (polygon_drawer):**
Activa el entorno virtual e instala `pip install -r requirements.txt`.
