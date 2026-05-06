## =========================================================================
## @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
## Taller 3 — Simplificación de Heightmap por Vecindad de Orden k
## =========================================================================

cmake_minimum_required(VERSION 3.14)
project(Taller3Heightmap)

set(CMAKE_CXX_STANDARD 17)

# -- pujCGAL headers
include_directories(
  ${CMAKE_SOURCE_DIR}/lib
  )

# -- Packages
find_package(CGAL REQUIRED)
find_package(PNG REQUIRED)

# -- Executables
add_executable(taller3 src/heightmap.cxx)
target_link_libraries(taller3 CGAL::CGAL PNG::PNG)

add_executable(comparador src/comparador.cxx)
target_link_libraries(comparador CGAL::CGAL PNG::PNG)

## eof - CMakeLists.txt
