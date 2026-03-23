## =========================================================================
## @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
## =========================================================================

# -- Basic configuration
cmake_minimum_required(VERSION 3.10)
project(SegmentsIntersections)
include("./baseConf.cmake")

include_directories(
  ${CMAKE_BINARY_DIR}/lib
  ${CMAKE_SOURCE_DIR}/lib
  )

# -- Packages
find_package(CGAL REQUIRED COMPONENTS Core)

# -- Executables
add_executable(pujCGAL_triangulate src/triangulate.cxx)
target_link_libraries(pujCGAL_triangulate CGAL::CGAL)

add_executable(pujCGAL_dual_graph src/dual_graph.cxx)
target_link_libraries(pujCGAL_dual_graph CGAL::CGAL)

## eof - CMakeLists.txt
