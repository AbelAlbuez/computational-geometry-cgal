## =========================================================================
## @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
## =========================================================================

# -- Basic configuration
cmake_minimum_required(VERSION 3.6)
project(SegmentsIntersections)
include("./baseConf.cmake")

include_directories(
  ${CMAKE_BINARY_DIR}/lib
  ${CMAKE_SOURCE_DIR}/lib
  )

# -- Packages
find_package(CGAL REQUIRED COMPONENTS Core)

# -- Executables
add_executable(no_cgal_brute_force src/no_cgal/brute_force.cxx)

add_executable(cgal_brute_force src/cgal/brute_force.cxx)
target_link_libraries(cgal_brute_force CGAL::CGAL)

add_executable(cgal_brute_force_2 src/cgal/brute_force_2.cxx)
target_link_libraries(cgal_brute_force_2 CGAL::CGAL)

add_executable(cgal_compare src/cgal/compare.cxx)
target_link_libraries(cgal_compare CGAL::CGAL)

## eof - CMakeLists.txt
