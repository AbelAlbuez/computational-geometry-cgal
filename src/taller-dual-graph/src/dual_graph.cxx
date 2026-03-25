// =========================================================================
// Taller 2 — Grafo Dual de un Polígono Simple
// Uso: ./dual_graph input.obj triangulation.obj dual.obj
// @author Santiago Gil Gallego (santiago_gil@javeriana.edu.co)
// @author Abel albueez (aa-albuezs@javeriana.edu.co)
// =========================================================================

#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <pujCGAL/IO.h>
#include <pujCGAL/Polygon.h>
#include <pujCGAL/Triangulation.h>
#include <pujCGAL/DualGraph.h>
#include <pujCGAL/IO_DualGraph.h>

int main( int argc, char** argv )
{
  using TKernel        = CGAL::Exact_predicates_inexact_constructions_kernel;
  using TPolygon       = pujCGAL::Polygon< TKernel >;
  using TTriangulation = pujCGAL::Triangulation< TKernel >;
  using TDualGraph     = pujCGAL::DualGraph< TKernel >;
  using TPoint         = TKernel::Point_2;

  TPolygon       polygon;
  TTriangulation mesh;
  TDualGraph     dual;

  // -- a. Lee polígono
  pujCGAL::IO::read( argv[ 1 ], polygon );

  // -- b. asegura: CCW orientation
  polygon.guarantee_CCW( );

  // -- c. Triangula
  pujCGAL::triangulate( mesh, polygon );

  // -- d. guarda triangulacion
  pujCGAL::IO::save( argv[ 2 ], mesh );

  // -- e. Computa los: barycenters
  std::vector< TPoint > barycenters;
  for(
    auto tIt = mesh.topology_begin( );
    tIt != mesh.topology_end( );
    ++tIt
    )
  {
    const TPoint& pa = mesh.point( ( *tIt )[ 0 ] );
    const TPoint& pb = mesh.point( ( *tIt )[ 1 ] );
    const TPoint& pc = mesh.point( ( *tIt )[ 2 ] );
    barycenters.push_back(
      TPoint(
        ( pa[ 0 ] + pb[ 0 ] + pc[ 0 ] ) / 3.0,
        ( pa[ 1 ] + pb[ 1 ] + pc[ 1 ] ) / 3.0
        )
      );
  } // end for

  // -- f. arma el grafo dual usando la funcion: build dual graph
  pujCGAL::build_dual_graph( dual, mesh, barycenters );

  // -- g. guarda grafo dual
  pujCGAL::IO::save( argv[ 3 ], dual );

  // -- h. síntesis
  std::size_t n_tri      = barycenters.size( );
  std::size_t n_boundary = 0;
  std::size_t n_internal = 0;
  for( std::size_t i = 0; i < n_tri; ++i )
  {
    if( dual.is_boundary( i ) ) n_boundary++;
    else                        n_internal++;
  } // end for

  std::cout << "Triangulos        : " << n_tri      << std::endl;
  std::cout << "  Frontera        : " << n_boundary << std::endl;
  std::cout << "  Internos        : " << n_internal << std::endl;
  std::cout << "Aristas internas  : " << dual.num_internal_edges( ) << std::endl;
  std::cout << "Aristas externas  : " << dual.num_external_edges( ) << std::endl;
  std::cout << "P-infinito        : ( "
            << dual.point_infinity( )[ 0 ] << ", "
            << dual.point_infinity( )[ 1 ] << " )" << std::endl;

  std::size_t n_nodes = dual.num_nodes( );
  std::cout << std::endl
            << "Matriz de adyacencia ("
            << n_nodes << " x " << n_nodes << "):" << std::endl;
  for( std::size_t i = 0; i < n_nodes; ++i )
  {
    for( std::size_t j = 0; j < n_nodes; ++j )
      std::cout << ( dual.adjacent( i, j ) ? 1 : 0 ) << " ";
    std::cout << std::endl;
  } // end for

  return( EXIT_SUCCESS );
}

// eof - dual_graph.cxx
