// =========================================================================
// Taller 2 — Grafo Dual de un Polígono Simple
// Paso 1: Cálculo de baricentros de la triangulación
// =========================================================================

#include <iostream>
#include <fstream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <pujCGAL/IO.h>
#include <pujCGAL/Polygon.h>
#include <pujCGAL/Triangulation.h>

int main( int argc, char** argv )
{
  using TKernel        = CGAL::Exact_predicates_inexact_constructions_kernel;
  using TPolygon       = pujCGAL::Polygon< TKernel >;
  using TTriangulation = pujCGAL::Triangulation< TKernel >;
  using TPoint         = TKernel::Point_2;

  // -- Leer polígono
  TPolygon polygon;
  TTriangulation mesh;
  pujCGAL::IO::read( argv[ 1 ], polygon );

  std::cout << "Area = " << polygon.area( ) << std::endl;
  polygon.guarantee_CCW( );
  std::cout << "Area CCW = " << polygon.area( ) << std::endl;

  // -- Triangular
  pujCGAL::triangulate( mesh, polygon );

  // -- Calcular baricentros
  std::vector< TPoint > barycenters;
  std::size_t tri_idx = 0;

  for( auto tIt = mesh.topology_begin( );
            tIt != mesh.topology_end( ); ++tIt, ++tri_idx )
  {
    const TPoint& pa = mesh.point( ( *tIt )[ 0 ] );
    const TPoint& pb = mesh.point( ( *tIt )[ 1 ] );
    const TPoint& pc = mesh.point( ( *tIt )[ 2 ] );

    TPoint bary(
      ( pa[ 0 ] + pb[ 0 ] + pc[ 0 ] ) / 3.0,
      ( pa[ 1 ] + pb[ 1 ] + pc[ 1 ] ) / 3.0
    );

    barycenters.push_back( bary );

    std::cout << "Triangulo " << tri_idx
              << " -> baricentro ("
              << bary[ 0 ] << ", " << bary[ 1 ] << ")"
              << std::endl;
  }

  std::cout << "Total triangulos: " << barycenters.size( ) << std::endl;

  // -- Escribir OBJ de salida con los baricentros como vertices
  std::ofstream ofs( argv[ 2 ] );
  ofs << "# Baricentros de la triangulacion" << std::endl;
  ofs << "# Un vertice por triangulo" << std::endl;
  for( const auto& b : barycenters )
    ofs << "v " << b[ 0 ] << " " << b[ 1 ] << " 0" << std::endl;
  ofs.close( );

  return( EXIT_SUCCESS );
}

// eof - dual_graph.cxx
