# Taller 3 - Simplificacion de Heightmap por Vecindad de Orden k

## Descripcion
Este programa lee un heightmap PNG, construye una triangulacion de Delaunay 2D y simplifica la malla eliminando vertices redundantes en zonas planas.
El criterio es el error de planitud local, procesado con min-heap (estilo Dijkstra) y vecindad de orden k con Face circulator (DCEL).

## Conceptos del curso aplicados
- DCEL / half-edge: recorrido de estrella con Face circulator.
- Dijkstra: min-heap para priorizar menor error de planitud.
- Delaunay: malla base y retriangulacion automatica tras T.remove.
- Face-vertex mesh: representacion topologica de vertices, caras y adyacencias.

## Compilacion
```bash
cd src/taller-3-codigo-base
mkdir -p build && cd build
cmake ..
make
```

## Dependencias Python (GIF)
```bash
pip3 install matplotlib Pillow numpy
```

## Ejecucion
Firma actual:
```bash
./taller3 input.png [epsilon] [orden_k]
```
Ejemplo:
```bash
./taller3 ../../../data/input_00.png 0.05 2
```

## Outputs generados
- output/original.obj: malla original antes de simplificar.
- output/simplificado.obj: malla final simplificada.
- output/resultado.gif: animacion comparando original vs simplificada.

## Visualizacion GIF
Al terminar, el ejecutable invoca automaticamente:
```bash
python3 ../visualizer.py
```
El script renderiza dos frames PNG y crea output/resultado.gif en loop infinito.

## Resultado de referencia (input_00.png)
| epsilon | orden_k | vertices antes | vertices despues | reduccion |
|---|---:|---:|---:|---:|
| 0.05 | 2 | 166500 | 13831 | 91.7% |
