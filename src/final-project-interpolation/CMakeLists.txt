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

# Eigen is required by the Poisson reconstruction + jet_estimate_normals
# (Stage 4 only). It is header-only and may be vendored without sudo; pass
# -DEIGEN3_INCLUDE_DIR=/path/to/eigen3 to cmake (see scripts/setup_eigen.sh).
find_package(Eigen3 3.3 REQUIRED)
include(CGAL_Eigen3_support)   # defines the CGAL::Eigen3_support imported target

# -- Include path so #include "ContourInterpolator.h" resolves
include_directories(${CMAKE_SOURCE_DIR}/src)

# -- Executable (main pulls in the Reconstructor -> needs Eigen)
add_executable(contour_interpolator src/main.cxx)
target_link_libraries(contour_interpolator CGAL::CGAL CGAL::Eigen3_support)

# -- Tests
enable_testing()

add_executable(test_linear src/test_linear.cxx)
target_link_libraries(test_linear CGAL::CGAL)
add_test(NAME linear_circle COMMAND test_linear)

add_executable(test_series src/test_series.cxx)
target_link_libraries(test_series CGAL::CGAL)
add_test(NAME series_spline_vs_poly COMMAND test_series)

add_executable(test_sdf src/test_sdf.cxx)
target_link_libraries(test_sdf CGAL::CGAL)
add_test(NAME sdf_circle_and_topology COMMAND test_sdf)

add_executable(test_reconstruct src/test_reconstruct.cxx)
target_link_libraries(test_reconstruct CGAL::CGAL CGAL::Eigen3_support)
add_test(NAME reconstruct_ellipsoid COMMAND test_reconstruct)

add_executable(test_metrics src/test_metrics.cxx)
target_link_libraries(test_metrics CGAL::CGAL CGAL::Eigen3_support)
add_test(NAME metrics_ellipsoid COMMAND test_metrics)

# -- pujCGAL header library (Bentley-Ottmann SegmentsIntersection), used by the
#    `linear` command for self-intersection resolution.
include_directories(${CMAKE_SOURCE_DIR}/../../lib)

# pujCGAL/SegmentsIntersection.hxx uses std::iter_value_t (C++20); the main
# executable therefore compiles as C++20 (the tests stay on C++17).
target_compile_features(contour_interpolator PRIVATE cxx_std_20)

# Register the shared header-only modules with the main target.
target_sources(contour_interpolator PRIVATE
  src/ContourResampler.h
  src/ContourResampler.hxx
  src/LinearInterpolator.h
  src/LinearInterpolator.hxx
  src/SelfIntersectionResolver.h
  src/SelfIntersectionResolver.hxx
)

## eof - CMakeLists.txt
