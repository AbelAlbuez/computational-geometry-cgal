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
│   └── ContourInterpolator.hxx
└── output/
```

## Pipeline

### 1. Extraer contornos desde BraTS (Python)

```bash
cd src/final-project-interpolation
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python scripts/extract_contours.py
```

Este paso:
- Copia los primeros 100 casos de `BraTS2024-BraTS-GLI-TrainingData.zip` (o
  de la carpeta ya descomprimida) a `data/raw/`, conservando solo el
  archivo de segmentación `*-seg.nii.gz`.
- Para cada slice axial con label ET (=3) y al menos 10 píxeles, extrae el
  contorno exterior más grande con `skimage.measure.find_contours` y lo
  escribe como polilínea cerrada en `data/contours/<caso>/slice_XXXX.obj`.
- Genera `data/contours/index.csv` con `case_id, slice_z, n_vertices, obj_path`.

### 2. Compilar el módulo C++ / CGAL

```bash
mkdir build && cd build
cmake ..
make
```

### 3. Probar el stub con dos slices consecutivos

```bash
./contour_interpolator \
  ../data/contours/BraTS-GLI-00000-000/slice_0070.obj \
  ../data/contours/BraTS-GLI-00000-000/slice_0071.obj
```

Por ahora solo confirma que ambos contornos se cargaron y reporta
`[TODO] Interpolacion pendiente`. La implementación geométrica vive en
`ContourInterpolator.{h,hxx}`.

## Formato `.obj` de los contornos

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

Los contornos son 2D (sin componente `z`); las aristas se indexan en base 1
y cierran el polígono con `l N 1`.

## Kernel CGAL

`CGAL::Exact_predicates_inexact_constructions_kernel` (suficiente para
predicados robustos sin sacrificar velocidad en construcciones).

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
!! Casos sin ET:         Z
```

### 4. Compilar el módulo C++/CGAL

```bash
mkdir build && cd build
cmake .. && make
```

### 5. Ejecutar el interpolador con dos slices

```bash
./contour_interpolator ../data/contours/BraTS-GLI-00000-000/slice_0070.obj \
                       ../data/contours/BraTS-GLI-00000-000/slice_0071.obj
```

> **Nota:** las carpetas `venv/`, `data/` y `output/` están en `.gitignore`
> y no se suben al repositorio.
