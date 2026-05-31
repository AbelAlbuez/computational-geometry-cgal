## =========================================================================
## Final project: Geometric interpolation of tumor contours
## Author: Abel Albuez Sanchez
## =========================================================================

cmake_minimum_required(VERSION 3.6)
project(ContourInterpolator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(${PROJECT_BINARY_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
  message(FATAL_ERROR "Build is not allowed in the same base path.")
endif()

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE
    )
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

# -- Packages
find_package(CGAL REQUIRED COMPONENTS Core)

# -- Include path so #include "ContourInterpolator.h" resolves
include_directories(${CMAKE_SOURCE_DIR}/src)

# -- Executable
add_executable(contour_interpolator src/main.cxx)
target_link_libraries(contour_interpolator CGAL::CGAL)

## eof - CMakeLists.txt

# -- pujCGAL header library (SegmentsIntersection, etc.)
include_directories(${CMAKE_SOURCE_DIR}/../../lib)

# pujCGAL/SegmentsIntersection.hxx uses std::iter_value_t (C++20).
target_compile_features(contour_interpolator PRIVATE cxx_std_20)

# -- Register new header-only modules with the existing target so IDEs and
#    dependency tracking pick them up.
target_sources(contour_interpolator PRIVATE
  src/interpolation-lineal/ContourResampler.h
  src/interpolation-lineal/ContourResampler.hxx
  src/interpolation-lineal/LinearInterpolator.h
  src/interpolation-lineal/LinearInterpolator.hxx
  src/interpolation-lineal/SelfIntersectionResolver.h
  src/interpolation-lineal/SelfIntersectionResolver.hxx
)
