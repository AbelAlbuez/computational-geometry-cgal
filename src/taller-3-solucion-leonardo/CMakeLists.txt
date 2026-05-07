## =========================================================================
## @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
## =========================================================================

# -- Basic configuration
cmake_minimum_required(VERSION 3.10)
project(DelaunayDecimation)
include("./baseConf.cmake")

include_directories(
  ${CMAKE_BINARY_DIR}/lib
  ${CMAKE_SOURCE_DIR}/lib
  )

# -- Packages
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(PNG REQUIRED)

# -- Executables
set(
  _execs
  build_delaunay_from_heightmap
  decimate_mesh
  )
foreach(_e ${_execs})
  add_executable(pujCGAL_${_e} src/${_e}.cxx)
  target_link_libraries(pujCGAL_${_e} PNG::PNG CGAL::CGAL)
endforeach()

add_executable(batch_runner src/batch_runner.cxx)
target_link_libraries(batch_runner CGAL::CGAL PNG::PNG)

## eof - CMakeLists.txt
