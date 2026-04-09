## =========================================================================
## First test — conteo de concavidades (polígono simple, pujCGAL)
## =========================================================================

cmake_minimum_required(VERSION 3.14)
project(FirstTest)

set(CMAKE_CXX_STANDARD 17)

# -- pujCGAL headers: parent of the pujCGAL/ folder
include_directories(
  ${CMAKE_SOURCE_DIR}/lib
  )

# -- Packages
find_package(CGAL REQUIRED)

find_package(VTK 9.0 REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  FiltersSources
  InteractionStyle
  IOImage
  RenderingAnnotation
  RenderingCore
  RenderingOpenGL2
  RenderingFreeType
  )

add_executable(concavidades src/concavidades.cxx src/concavidades_viz.cxx)
target_link_libraries(concavidades PRIVATE
  CGAL::CGAL
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::FiltersSources
  VTK::InteractionStyle
  VTK::IOImage
  VTK::RenderingAnnotation
  VTK::RenderingCore
  VTK::RenderingOpenGL2
  VTK::RenderingFreeType
  )
vtk_module_autoinit(
  TARGETS concavidades
  MODULES VTK::RenderingOpenGL2 VTK::RenderingFreeType
  )

add_executable(visualizer src/visualizer.cxx)
target_link_libraries(visualizer PRIVATE
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::FiltersSources
  VTK::InteractionStyle
  VTK::IOImage
  VTK::RenderingAnnotation
  VTK::RenderingCore
  VTK::RenderingOpenGL2
  VTK::RenderingFreeType
  )
vtk_module_autoinit(
  TARGETS visualizer
  MODULES VTK::RenderingOpenGL2 VTK::RenderingFreeType
  )

## eof - CMakeLists.txt
