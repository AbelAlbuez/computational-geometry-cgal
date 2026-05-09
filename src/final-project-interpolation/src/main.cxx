// =============================================================================
// Final project: Geometric Interpolation of Tumor Contours
// Author: Abel Albuez Sanchez
//
// Stub principal: lee dos contornos .obj (slice A y slice B) y confirma que
// fueron cargados. La interpolación geométrica se implementará después.
//
// Usage: ./contour_interpolator slice_A.obj slice_B.obj
// =============================================================================

#include <iostream>
#include <string>

#include "ContourInterpolator.h"

int main( int argc, char** argv )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ]
              << " slice_A.obj slice_B.obj" << std::endl;
    return 1;
  }

  using TInterp = pujCGAL::Final::ContourInterpolator;

  const std::string fa = argv[ 1 ];
  const std::string fb = argv[ 2 ];

  auto A = TInterp::read_obj( fa );
  auto B = TInterp::read_obj( fb );

  std::cout << "Slice A: " << fa << "  (" << A.size( ) << " vertices)\n";
  std::cout << "Slice B: " << fb << "  (" << B.size( ) << " vertices)\n";

  if( A.empty( ) || B.empty( ) )
  {
    std::cerr << "Error: one of the contours is empty." << std::endl;
    return 2;
  }

  std::cout << "[TODO] Interpolacion pendiente" << std::endl;
  return 0;
}

// eof - main.cxx
