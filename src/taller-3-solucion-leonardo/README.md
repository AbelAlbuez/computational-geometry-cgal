# Taller 3 - Solucion Leonardo

En este proyecto uso batch_runner para procesar todos los PNG del directorio data, construir la malla original, decimarla y dejar resultados numericos y visuales por imagen.

## Entorno Python

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

## Compilacion

cmake -S . -B build
cmake --build build -j

## Ejecucion

./build/batch_runner [data_dir] [order] [gamma]

Si no paso argumentos, usa data del proyecto, order=2 y gamma=1.0.

## Salidas

Por cada imagen genero:
- output/batch-<nombre>/original.obj
- output/batch-<nombre>/simplificado.obj
- output/batch-<nombre>/00_heightmap.png
- output/batch-<nombre>/01_original.png
- output/batch-<nombre>/02_simplificado.png
- output/batch-<nombre>/03_comparativo.png

Al final tambien dejo:
- output/batch-resumen/resumen.csv
- output/batch-resumen/bitacora.md
