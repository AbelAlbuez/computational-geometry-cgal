
#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <pujCGAL/Heightmap.h>
#include <pujCGAL/IO.h>

int main( int argc, char** argv )
{
  using TKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using TReal = typename TKernel::RT;
  using TVertices = CGAL::Triangulation_vertex_base_with_info_2< TReal, TKernel >;
  using TTriangulationDS = CGAL::Triangulation_data_structure_2< TVertices >;
  using TDelaunay = CGAL::Delaunay_triangulation_2< TKernel, TTriangulationDS >;
  using TPoint = TKernel::Point_2;

  pujCGAL::Heightmap< TReal > hm;
  hm.read_from_png( argv[ 1 ] );
  hm.set_origin( -10, -10 );
  hm.set_dimensions( 10, 10 );

  std::vector< TPoint > points;
  for( std::size_t h = 0; h < hm.height( ); ++h )
  {
    for( std::size_t w = 0; w < hm.width( ); ++w )
    {
      auto p = hm.point( w, h );
      points.push_back( TPoint( p.first, p.second ) );
    } // end for
  } // end for

  TDelaunay T;
  T.insert( points.begin( ), points.end( ) );
  for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
  {
    auto i = hm.index( v->point( ).x( ), v->point( ).y( ) );
    v->info( ) = hm( i.first, i.second );
  } // end for

  pujCGAL::IO::save( T, argv[ 2 ] );
    
  return( EXIT_SUCCESS );
}

// eof - build_delaunay_from_heightmap.cxx
