# Interpolación Geométrica e Híbrida de Contornos en Imágenes Médicas

Proyecto final — Geometría Computacional, Maestría en Ingeniería de Sistemas y Computación, Pontificia Universidad Javeriana.

Este proyecto implementa un pipeline híbrido para la interpolación de contornos tumorales 2D extraídos de resonancias magnéticas axiales (**BraTS 2024 GLI**). El sistema combina técnicas de geometría computacional avanzada en C++ (CGAL) con modelos de aprendizaje profundo para la interpolación de fotogramas de video (Video Frame Interpolation - VFI).

---

## 🛠️ Estado Actual del Proyecto y División del Trabajo

El desarrollo se encuentra dividido en módulos funcionales distribuidos eficientemente entre el equipo:

1. **🟡 Jesús  — Correspondencia de Vértices e Integración de IA (Orquestador Principal):**
   * **Correspondencia por Ángulo Polar:** Implementación en `ContourInterpolator.hxx` que calcula el centroide de cada contorno y reordena los vértices circularmente basándose en su ángulo polar (usando `std::atan2`). Esto sincroniza el "punto cero" del recorrido de los polígonos, eliminando por completo el efecto de *twisting* (colapso geométrico/retorcimiento) antes de la interpolación lineal.
   * **Soporte Diferencial de Índices:** Manejo proporcional de índices para permitir la interpolación directa entre slices con diferente resolución de vértices (p. ej., un slice de 32 vértices frente a uno de 30).
   * **Integración del Baseline de IA:** Creación del script ejecutor en Python que rasteriza la geometría discreta a máscaras de píxeles, gestiona las dimensiones con acolchado dinámico (multiplos de 32) e inyecta los tensores en una red de flujo óptico neuronal (`BaselineIFNet`).
   * **Visualización Inteligente (Auto-Zoom):** Incorporación de un algoritmo de detección de Bounding Box con margen adaptativo de 15 píxeles que recorta el tumor interpolado por IA (de tamaño milimétrico en el lienzo nativo de 240x240) y lo amplifica 8x mediante interpolación por vecino más cercano, facilitando su evaluación visual.

2. **🟢 Abel — Interpolación y Validación Geométrica:**
   * Recepción de los contornos sincronizados por ángulo polar para realizar la interpolación lineal paramétrica en el espacio continuo (t en [0, 1]).
   * Implementación de la validación topológica mediante el algoritmo de Bentley-Ottmann (CGAL) para la detección y limpieza de auto-intersecciones en geometrías complejas no convexas.

3. **🔵 Santiago — Reconstrucción Superficial y Visualización 3D:**
   * Pipeline de apilamiento secuencial de los archivos `.obj` intermedios e interpolados respetando su altura espacial indexada.
   * Construcción de mallas poligonales complejas de la superficie tumoral utilizando el toolkit de visualización VTK.

---

## 📁 Estructura del Proyecto

```text
final-project-interpolation/
├── CMakeLists.txt
├── README.md
├── requirements.txt
├── .gitignore          # Configurado para excluir entornos virtuales y archivos .pth pesados
├── data/
│   ├── raw/            # Datos crudos de segmentación (*-seg.nii.gz)
│   └── contours/       # Contornos extraídos en formato .obj indexados por caso
├── modelos/
│   └── modelo_baseline_100_definitivo.pth  # Pesos locales de la red neuronal (omitidos en Git)
├── scripts/
│   ├── extract_contours.py       # Etapa 1: Extractor de contornos de la MRI a .obj 3D
│   └── inferencia_ia_contours.py # Etapa 3: Rasterizador y ejecutor del modelo neuronal
└── src/
    ├── main.cxx                  # Orquestador C++ (Carga datos, calcula geometría e invoca IA)
    ├── ContourInterpolator.h     # Interfaz pública de la clase y definición del Kernel CGAL
    └── ContourInterpolator.hxx   # Implementación de I/O, Centroides y Correspondencia Polar
```

---

## 🚀 Requisitos del Sistema

### Dependencias Globales de C++
Es necesario instalar la librería de desarrollo de CGAL y las herramientas de compilación base directamente en el sistema operativo:
```bash
sudo apt update
sudo apt install build-essential cmake libcgal-dev
```

### Dependencias del Entorno Virtual (Python)
El entorno requiere versiones específicas que garanticen la compatibilidad de la interfaz binaria (ABI) entre los tensores avanzados de PyTorch y las estructuras de arrays numéricos:
* **Python 3.10 o superior**
* **PyTorch >= 2.3.0** (Garantiza soporte nativo estable para las ramas más recientes de NumPy)

---

## 💻 Guía de Ejecución Paso a Paso

### 1. Configurar el Entorno Virtual e Instalar Dependencias
Desde la raíz del proyecto (`final-project-interpolation/`), crea el entorno virtual e instala los paquetes necesarios. Si cuentas con una GPU NVIDIA (arquitectura LOQ), instala la distribución CUDA para acelerar el procesamiento de la IA:

```bash
# Crear el entorno virtual
python3 -m venv venv

# Activar el entorno virtual
source venv/bin/activate

# Instalar PyTorch optimizado para CUDA (Recomendado para GPU NVIDIA)
pip install torch torchvision torchaudio --index-url [https://download.pytorch.org/whl/cu118](https://download.pytorch.org/whl/cu118)

# Instalar el resto de dependencias del ecosistema médico y visión
pip install -r requirements.txt
```

### 2. Extracción de Contornos Compatibles con Visores 3D
Ejecuta el script de extracción. Este proceso genera archivos `.obj` estructurados con tres coordenadas `(X, Y, Z)` por cada vértice, asignando a la componente `Z` la altura real del slice axial del cerebro. Esto previene fallos de parseo en herramientas de análisis volumétrico como **ParaView** o la jerarquía ejecutable de **VTK**:

```bash
python scripts/extract_contours.py
```

### 3. Compilación del Módulo de Geometría en C++
Utiliza CMake para generar el MakeFile nativo y compilar el binario orquestador:

```bash
# Crear y entrar al directorio de construcción
mkdir -p build && cd build

# Configurar el proyecto con CMake
cmake ..

# Compilar el binario
make
```

### 4. Ejecución del Pipeline Híbrido de Prueba
Una vez compilado el ejecutable `contour_interpolator`, córrelo pasándole como argumentos dos archivos `.obj` que sean slices consecutivos y contengan datos válidos del tumor (por ejemplo, los cortes 70 y 71 del caso de prueba `BraTS-GLI-00008-100`):

```bash
./contour_interpolator \
  ../data/contours/BraTS-GLI-00008-100/slice_0070.obj \
  ../data/contours/BraTS-GLI-00008-100/slice_0071.obj
```

---

## 📊 Formato de Resultados Generados

Tras la ejecución exitosa, los resultados se almacenarán de forma automatizada directamente en el directorio del paciente evaluado (`data/contours/BraTS-GLI-00008-100/`):

1. **`geom_interp_mid.obj`:** Contorno de malla poligonal generado al 50% de la transición temporal (t=0.5). Contiene los vértices ordenados por ángulo polar y procesados geométricamente de forma continua. Puede abrirse directamente en **ParaView** para su renderizado espacial.
2. **`baseline_mid_slice_0070_slice_0071.png`:** Máscara de píxeles en resolución absoluta de 240x240 que representa el espacio físico real donde la red neuronal infiere que se posicionará la masa tumoral en el plano intermedio.
3. **`baseline_mid_slice_0070_slice_0071_zoomed.png`:** Imagen de inspección visual recortada dinámicamente sobre la zona de interés del tumor y amplificada 8x. Muestra los detalles de interpolación fina y fronteras generados por el modelo de aprendizaje profundo, ideal para realizar las comparativas cualitativas frente al contorno geométrico de CGAL.