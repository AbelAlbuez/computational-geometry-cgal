# Interpolación Geométrica de Contornos en Imágenes Médicas

Proyecto final — Geometría Computacional, Maestría en Ingeniería de Sistemas
y Computación, Pontificia Universidad Javeriana.

Objetivo: interpolar geométricamente contornos 2D del tumor activo (label ET)
entre slices axiales consecutivos del dataset **BraTS 2024 GLI** usando CGAL.

## Estructura

```
final-project-interpolation/
├── CMakeLists.txt
├── README.md
├── requirements.txt
├── data/
│   ├── raw/          # 100 casos BraTS (solo *-seg.nii.gz)
│   └── contours/     # contornos 2D por slice en formato .obj + index.csv
├── scripts/
│   └── extract_contours.py
├── src/
│   ├── main.cxx
│   ├── ContourInterpolator.h
│   ├── ContourInterpolator.hxx
│   ├── ContourResampler.h
│   ├── ContourResampler.hxx
│   ├── LinearInterpolator.h
│   ├── LinearInterpolator.hxx
│   ├── SelfIntersectionResolver.h
│   └── SelfIntersectionResolver.hxx
└── output/
```

## Pipeline

1. **Resampling por longitud de arco** (`ContourResampler`) — iguala
   `|A|` y `|B|` redistribuyendo `n` vértices uniformemente sobre los
   segmentos del contorno cerrado.
2. **Interpolación lineal** (`LinearInterpolator`) —
   `P(t)[i] = (1-t)·A[i] + t·B[i]` para cada vértice `i`.
3. **Construcción de segmentos CGAL** — el contorno interpolado se convierte
   en un `std::vector<Segment_2>` conectando vértices consecutivos en
   orden circular.
4. **Detección de auto-intersecciones con Bentley-Ottmann**
   (`SelfIntersectionResolver`) — usa
   `pujCGAL::SegmentsIntersection::BentleyOttmann` para reportar todos los
   puntos de cruce reales (excluyendo extremos compartidos entre aristas
   consecutivas).
5. **Resolución de cruces** — el contorno se subdivide en los puntos de
   cruce y se extraen los sub-lazos mediante recorrido con pila; se
   retiene el lazo con mayor área absoluta.
6. **Escritura del contorno limpio** en `output/contour_interpolated.obj`
   (formato 3D con `z = 0.0` para compatibilidad con visores como
   ParaView, Online 3D Viewer e ImageToSTL).

## Módulos implementados

| Clase | Archivo | Responsabilidad |
|---|---|---|
| `ContourResampler` | `ContourResampler.h/.hxx` | Resampling por longitud de arco |
| `LinearInterpolator` | `LinearInterpolator.h/.hxx` | Interpolación lineal vértice a vértice |
| `SelfIntersectionResolver` | `SelfIntersectionResolver.h/.hxx` | Bentley-Ottmann + resolución de cruces |

## Formato `.obj` de los contornos de entrada

```
# Contorno tumor ET - caso <id> - slice axial <z>
# <N> vertices
v x0 y0
v x1 y1
...
v x{N-1} y{N-1}
l 1 2
l 2 3
...
l N 1
```

Los contornos de entrada son 2D (sin componente `z`); las aristas se
indexan en base 1 y cierran el polígono con `l N 1`.

## Formato `.obj` de salida (interpolado)

```
# Contorno interpolado - 32 vertices
v 63.0000 50.5000 0.0000
v 62.3092 50.2268 0.0000
...
l 1 2
l 2 3
...
l 32 1
```

Nota: los vértices de salida incluyen `z=0.0` para compatibilidad con
visores 3D como ParaView, Online 3D Viewer e ImageToSTL.

## Kernel CGAL

`CGAL::Exact_predicates_inexact_constructions_kernel` (suficiente para
predicados robustos sin sacrificar velocidad en construcciones).

## Requisitos

- CGAL (`find_package(CGAL REQUIRED COMPONENTS Core)`)
- **C++20** (`pujCGAL/SegmentsIntersection.hxx` usa `std::iter_value_t`).
  El `CMakeLists.txt` ya incluye
  `target_compile_features(contour_interpolator PRIVATE cxx_std_20)`.
- Python 3 con `numpy`, `nibabel`, `scikit-image` (ver `requirements.txt`).

## Ejecución

### 1. Crear y activar el entorno virtual

Este proyecto usa un entorno virtual Python. Siempre activarlo antes de
correr cualquier script.

```bash
# Crear el entorno (solo la primera vez)
python3 -m venv venv

# Activar (cada vez que abras una terminal nueva)
source venv/bin/activate   # macOS / Linux
venv\Scripts\activate      # Windows
```

### 2. Instalar dependencias (solo la primera vez)

```bash
pip install -r requirements.txt
```

### 3. Extraer contornos desde BraTS

```bash
python scripts/extract_contours.py
```

Al finalizar verás:

```
OK Casos procesados:   X/100
OK Contornos exportados: Y archivos .obj
OK Index guardado en:    data/contours/index.csv
```

### 4. Compilar el módulo C++/CGAL

```bash
mkdir build && cd build
cmake .. && make
```

### 5. Ejecutar el interpolador con dos slices

```bash
./contour_interpolator ../data/contours/BraTS-GLI-00008-100/slice_0070.obj \
                       ../data/contours/BraTS-GLI-00008-100/slice_0071.obj
```

Output esperado:

```
Contour A: 32 vertices
Contour B: 30 vertices
Resampled to: 32 vertices
Interpolated contour: 32 vertices
Self-intersections detected: no
Result written to output/contour_interpolated.obj
```

> **Nota:** las carpetas `venv/`, `data/` y `output/` están en `.gitignore`
> y no se suben al repositorio.
