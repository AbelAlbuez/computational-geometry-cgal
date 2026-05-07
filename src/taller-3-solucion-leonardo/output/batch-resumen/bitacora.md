# Bitacora -- Batch runner solucion Leonardo (Taller 3)

Fecha: 2026-05-06
Parametros: order=2, gamma=1

## Pipeline
En este batch replique el pipeline del profesor en dos etapas. Primero construyo el Delaunay completo desde el heightmap y guardo original.obj. Luego leo ese OBJ como Surface_mesh y aplico el criterio estadistico |h - mu| < gamma*sigma sobre vecindad de orden k para filtrar vertices. Con los sobrevivientes reconstruyo Delaunay y guardo simplificado.obj.

## Imagenes procesadas
- SRTM_US_scaled_512.png (512x256 = 131072 pixeles)

## Resultados por imagen

### SRTM_US_scaled_512
| vertices antes -> despues | reduccion % | tiempo etapa 1 (ms) | tiempo etapa 2 (ms) |
|---:|---:|---:|---:|
| 131072 -> 49581 | 62.17 | 284.67 | 499.52 |

En esta imagen consegui una reduccion de 62.17%. La etapa 2 tomo mas tiempo que la etapa 1. El resultado me parece razonable para el criterio de filtrado usado.

PNGs: 4/4 generados (00_heightmap.png, 01_original.png, 02_simplificado.png, 03_comparativo.png)

## Observaciones globales

- La mayor reduccion fue en SRTM_US_scaled_512 con 62.17%.
- La menor reduccion fue en SRTM_US_scaled_512 con 62.17%.
- La imagen mas grande fue SRTM_US_scaled_512 (131072 pixeles).
- La imagen mas pequena fue SRTM_US_scaled_512 (131072 pixeles).
- Tiempo promedio etapa 1: 284.67 ms; etapa 2: 499.52 ms.
- Imagenes con error durante el batch: 0.
