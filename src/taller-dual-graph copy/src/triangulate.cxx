// =========================================================================
// @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
// =========================================================================

#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <pujCGAL/IO.h>
#include <pujCGAL/Polygon.h>
#include <pujCGAL/Triangulation.h>
 
int main( int argc, char** argv )
{
  using TKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using TPolygon = pujCGAL::Polygon< TKernel >;
  using TTriangulation = pujCGAL::Triangulation< TKernel >;
  using TPoint = TPolygon::TPoint;

  TPolygon polygon;
  TTriangulation mesh;
  pujCGAL::IO::read( argv[ 1 ], polygon );

  std::cout << "Area = " << polygon.area( ) << std::endl;
  polygon.guarantee_CCW( );
  std::cout << "Corrected area = " << polygon.area( ) << std::endl;

  pujCGAL::triangulate( mesh, polygon );
  pujCGAL::IO::save( argv[ 2 ], mesh );

  return( EXIT_SUCCESS );
}

// eof - triangulate.cxx
