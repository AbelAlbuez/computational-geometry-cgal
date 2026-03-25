## =========================================================================
## @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
## =========================================================================

cmake_minimum_required(VERSION 3.14)
project(DualGraph)

set(CMAKE_CXX_STANDARD 17)

# -- pujCGAL headers: parent of the pujCGAL/ folder
include_directories(
  ${CMAKE_SOURCE_DIR}/lib
  )

# -- Packages
find_package(CGAL REQUIRED)

# -- Executables
add_executable(pujCGAL_triangulate src/triangulate.cxx)
target_link_libraries(pujCGAL_triangulate CGAL::CGAL)

add_executable(dual_graph src/dual_graph.cxx)
target_link_libraries(dual_graph CGAL::CGAL)

## eof - CMakeLists.txt
