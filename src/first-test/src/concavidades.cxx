// =========================================================================
// Parcial 1 — Concavidades en nube de puntos (.obj solo v).
// Cinco pasos según enunciado + VTK (6 PNG + GIF).
// Uso: ./concavidades entrada.obj
// =========================================================================

#include "concavidades_viz.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Segment_2.h>
#include <CGAL/convex_hull_2.h>

#include <pujCGAL/DualGraph.h>
#include <pujCGAL/Polygon.h>
#include <pujCGAL/Triangulation.h>

namespace
{
  using TKernel        = CGAL::Exact_predicates_inexact_constructions_kernel;
  using TPoint         = TKernel::Point_2;
  using TPolygon       = pujCGAL::Polygon< TKernel >;
  using TTriangulation = pujCGAL::Triangulation< TKernel >;
  using TDualGraph     = pujCGAL::DualGraph< TKernel >;
  using Edge           = std::pair< int, int >;

  bool
  read_point_cloud_obj( const char* fname, std::vector< TPoint >& out )
  {
    std::ifstream ifs( fname );
    if( !ifs )
      return( false );
    std::string line;
    while( std::getline( ifs, line ) )
    {
      if( line.empty( ) )
        continue;
      if( line[ 0 ] == '#' || line[ 0 ] == 'o' )
        continue;
      if( line.size( ) < 2 || line[ 0 ] != 'v' || line[ 1 ] != ' ' )
        continue;
      double x, y;
      std::istringstream iss( line.substr( 2 ) );
      iss >> x >> y;
      out.emplace_back( x, y );
    } // end while
    return( !out.empty( ) );
  }

  void
  polygon_from_ring( TPolygon& poly, const std::vector< TPoint >& ring )
  {
    for( const auto& p : ring )
      poly.add_point( p );
    poly.build( );
  }

  std::pair< double, double >
  to_pair( const TPoint& p )
  {
    return(
      std::make_pair(
        CGAL::to_double( p.x( ) ), CGAL::to_double( p.y( ) ) )
      );
  }

  std::vector< TPoint >
  poligono_ajustado_ring( const std::vector< TPoint >& cloud )
  {
    TKernel::FT sx = 0, sy = 0;
    for( const auto& p : cloud )
    {
      sx += p.x( );
      sy += p.y( );
    } // end for
    const TKernel::FT invn = TKernel::FT( 1 ) / TKernel::FT( cloud.size( ) );
    const TKernel::FT cy   = sy * invn;

    std::vector< TPoint > sorted = cloud;
    std::sort(
      sorted.begin( ),
      sorted.end( ),
      []( const TPoint& a, const TPoint& b ) {
        return( a.x( ) < b.x( ) || ( a.x( ) == b.x( ) && a.y( ) < b.y( ) ) );
      }
      );

    std::vector< TPoint > upper;
    std::vector< TPoint > lower;
    for( const auto& p : sorted )
    {
      if( p.y( ) >= cy )
        upper.push_back( p );
      else
        lower.push_back( p );
    } // end for
    std::sort(
      lower.begin( ),
      lower.end( ),
      []( const TPoint& a, const TPoint& b ) {
        return( a.x( ) > b.x( ) || ( a.x( ) == b.x( ) && a.y( ) > b.y( ) ) );
      }
      );

    std::vector< TPoint > ring;
    ring.reserve( upper.size( ) + lower.size( ) );
    ring.insert( ring.end( ), upper.begin( ), upper.end( ) );
    ring.insert( ring.end( ), lower.begin( ), lower.end( ) );
    return( ring );
  }

  std::vector< TPoint >
  puntos_poligono_ccw( const TPolygon& poly )
  {
    std::vector< TPoint > out;
    for(
      auto pIt = poly.topology_begin( );
      pIt != poly.topology_end( );
      ++pIt
      )
      out.push_back( poly.point( *pIt ) );
    return( out );
  }

  bool
  arista_en_ch(
    const TPoint& a,
    const TPoint& b,
    const std::vector< TPoint >& ch_ring
    )
  {
    using TSegment = CGAL::Segment_2< TKernel >;
    const std::size_t n = ch_ring.size( );
    for( std::size_t i = 0; i < n; ++i )
    {
      const TPoint& u = ch_ring[ i ];
      const TPoint& v = ch_ring[ ( i + 1 ) % n ];
      TSegment s( u, v );
      if( s.has_on( a ) && s.has_on( b ) )
        return( true );
    } // end for
    return( false );
  }

  std::map< Edge, int >
  contar_aristas_triangulos( const TTriangulation& mesh )
  {
    std::map< Edge, int > cnt;
    auto bump = [&]( int a, int b ) {
      Edge e = { std::min( a, b ), std::max( a, b ) };
      cnt[ e ]++;
    };
    for(
      auto tIt = mesh.topology_begin( );
      tIt != mesh.topology_end( );
      ++tIt
      )
    {
      int a = static_cast< int >( ( *tIt )[ 0 ] );
      int b = static_cast< int >( ( *tIt )[ 1 ] );
      int c = static_cast< int >( ( *tIt )[ 2 ] );
      bump( a, b );
      bump( b, c );
      bump( a, c );
    } // end for
    return( cnt );
  }

  void
  aristas_borde_por_triangulo(
    const TTriangulation& mesh,
    const std::map< Edge, int >& ecnt,
    std::vector< std::vector< Edge > >& por_tri
    )
  {
    const std::size_t nt = static_cast< std::size_t >(
      std::distance( mesh.topology_begin( ), mesh.topology_end( ) )
      );
    por_tri.assign( nt, {} );
    std::size_t ti = 0;
    for(
      auto tIt = mesh.topology_begin( );
      tIt != mesh.topology_end( );
      ++tIt, ++ti
      )
    {
      int a = static_cast< int >( ( *tIt )[ 0 ] );
      int b = static_cast< int >( ( *tIt )[ 1 ] );
      int c = static_cast< int >( ( *tIt )[ 2 ] );
      Edge e0 = { std::min( a, b ), std::max( a, b ) };
      Edge e1 = { std::min( b, c ), std::max( b, c ) };
      Edge e2 = { std::min( a, c ), std::max( a, c ) };
      if( ecnt.at( e0 ) == 1 )
        por_tri[ ti ].push_back( e0 );
      if( ecnt.at( e1 ) == 1 )
        por_tri[ ti ].push_back( e1 );
      if( ecnt.at( e2 ) == 1 )
        por_tri[ ti ].push_back( e2 );
    } // end for
  }

  std::string
  fmt_visitado( const std::vector< bool >& v )
  {
    std::ostringstream o;
    o << "[";
    for( std::size_t i = 0; i < v.size( ); ++i )
    {
      if( i )
        o << ", ";
      o << ( v[ i ] ? "true" : "false" );
    } // end for
    o << "]";
    return( o.str( ) );
  }

} // namespace

int main( int argc, char** argv )
{
  if( argc != 2 )
  {
    std::cerr << "Uso: " << argv[ 0 ] << "  entrada.obj\n";
    return( EXIT_FAILURE );
  }

  std::vector< TPoint > cloud;
  if( !read_point_cloud_obj( argv[ 1 ], cloud ) )
  {
    std::cerr << "No se pudo leer la nube: " << argv[ 1 ] << std::endl;
    return( EXIT_FAILURE );
  }

  // -- Paso 1: casco convexo
  std::vector< TPoint > ch_ring;
  CGAL::convex_hull_2(
    cloud.begin( ),
    cloud.end( ),
    std::back_inserter( ch_ring )
    );
  TPolygon poly_ch;
  polygon_from_ring( poly_ch, ch_ring );
  poly_ch.guarantee_CCW( );
  std::cout << "[Paso 1] CH: " << ch_ring.size( ) << " vértices" << std::endl;

  // -- Paso 2: polígono ajustado
  std::vector< TPoint > adj_ring = poligono_ajustado_ring( cloud );
  TPolygon poly_adj;
  polygon_from_ring( poly_adj, adj_ring );
  poly_adj.guarantee_CCW( );
  std::cout << "[Paso 2] Polígono ajustado: " << adj_ring.size( )
            << " vértices, CCW" << std::endl;

  // -- Paso 3: triangulación
  TTriangulation mesh;
  pujCGAL::triangulate( mesh, poly_adj );
  const std::size_t n_tri = static_cast< std::size_t >(
    std::distance( mesh.topology_begin( ), mesh.topology_end( ) )
    );
  std::vector< TPoint > barycenters;
  barycenters.reserve( n_tri );
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
  std::cout << "[Paso 3] Triangulación: " << n_tri << " triángulos"
            << std::endl;

  // -- Paso 4: grafo dual y bolsillos
  TDualGraph dual;
  pujCGAL::build_dual_graph( dual, mesh, barycenters );

  std::map< Edge, int > ecnt = contar_aristas_triangulos( mesh );
  std::vector< std::vector< Edge > > borde_por_tri;
  aristas_borde_por_triangulo( mesh, ecnt, borde_por_tri );

  std::vector< bool > is_bolsillo( n_tri, false );
  for( std::size_t i = 0; i < n_tri; ++i )
  {
    if( !dual.is_boundary( i ) )
      continue;
    bool alguna_no_en_ch = false;
    for( const Edge& e : borde_por_tri[ i ] )
    {
      const TPoint& pa = mesh.point( static_cast< std::size_t >( e.first ) );
      const TPoint& pb = mesh.point( static_cast< std::size_t >( e.second ) );
      if( !arista_en_ch( pa, pb, ch_ring ) )
        alguna_no_en_ch = true;
    } // end for
    is_bolsillo[ i ] = alguna_no_en_ch;
  } // end for

  std::cout << "[Paso 4] Grafo dual construido" << std::endl;
  std::cout << "         Bolsillos: [";
  {
    bool first = true;
    for( std::size_t i = 0; i < n_tri; ++i )
    {
      if( !is_bolsillo[ i ] )
        continue;
      if( !first )
        std::cout << ", ";
      first = false;
      std::cout << i;
    } // end for
  }
  std::cout << "]" << std::endl;

  // -- Paso 5: contar concavidades (BFS en dual restringido a bolsillos)
  std::vector< bool > visitado( n_tri, false );
  std::queue< int >    cola;
  int                  contador = 0;
  const int            ntri_i   = static_cast< int >( n_tri );
  std::vector< int >   pocket_component( n_tri, -1 );

  std::cout << "[Paso 5] visitado[] inicial: " << fmt_visitado( visitado )
            << std::endl;

  for( int i = 0; i < ntri_i; ++i )
  {
    if( !is_bolsillo[ static_cast< std::size_t >( i ) ] )
      continue;
    if( visitado[ static_cast< std::size_t >( i ) ] )
      continue;

    contador++;
    visitado[ static_cast< std::size_t >( i ) ] = true;
    pocket_component[ static_cast< std::size_t >( i ) ] = contador - 1;
    cola.push( i );
    std::vector< int > nodos_comp;

    while( !cola.empty( ) )
    {
      int nodo = cola.front( );
      cola.pop( );
      nodos_comp.push_back( nodo );
      for( int j = 0; j < ntri_i; ++j )
      {
        if(
          !dual.adjacent(
            static_cast< std::size_t >( nodo ),
            static_cast< std::size_t >( j )
            )
          )
          continue;
        if( !is_bolsillo[ static_cast< std::size_t >( j ) ] )
          continue;
        if( visitado[ static_cast< std::size_t >( j ) ] )
          continue;
        visitado[ static_cast< std::size_t >( j ) ] = true;
        pocket_component[ static_cast< std::size_t >( j ) ] = contador - 1;
        cola.push( j );
      } // end for j
    } // end while

    std::cout << "Concavidad " << contador << ": nodos B = [";
    for( std::size_t k = 0; k < nodos_comp.size( ); ++k )
    {
      if( k )
        std::cout << ", ";
      std::cout << nodos_comp[ k ];
    } // end for
    std::cout << "]" << std::endl;
    std::cout << "visitado[] ahora: " << fmt_visitado( visitado ) << std::endl;
  } // end for i

  std::cout << "────────────────────────────────────" << std::endl;
  std::cout << "Resultado: " << contador << " concavidad(es)" << std::endl;
  std::cout << "────────────────────────────────────" << std::endl;

  // -- VTK
  const std::vector< TPoint > ch_draw  = puntos_poligono_ccw( poly_ch );
  const std::vector< TPoint > adj_draw = puntos_poligono_ccw( poly_adj );

  ConcavidadVizInput vin;
  for( const auto& p : cloud )
    vin.cloud.push_back( to_pair( p ) );
  for( const auto& p : ch_draw )
    vin.ch_ccw.push_back( to_pair( p ) );
  for( const auto& p : adj_draw )
    vin.adj_ccw.push_back( to_pair( p ) );
  for( auto gIt = mesh.geometry_begin( ); gIt != mesh.geometry_end( ); ++gIt )
    vin.mesh_points.push_back(
      std::make_pair(
        CGAL::to_double( gIt->x( ) ), CGAL::to_double( gIt->y( ) ) )
      );
  for(
    auto tIt = mesh.topology_begin( );
    tIt != mesh.topology_end( );
    ++tIt
    )
    vin.mesh_tris.push_back(
      { ( *tIt )[ 0 ], ( *tIt )[ 1 ], ( *tIt )[ 2 ] }
      );
  for( const auto& p : barycenters )
    vin.barycenters.push_back( to_pair( p ) );
  vin.p_inf = to_pair( dual.point_infinity( ) );

  for( auto it = dual.internal_edges_begin( ); it != dual.internal_edges_end( );
       ++it )
    vin.dual_int_edges.push_back( { ( *it ).first, ( *it ).second } );
  const std::size_t p_inf_idx = n_tri;
  for( auto it = dual.external_edges_begin( ); it != dual.external_edges_end( );
       ++it )
    vin.dual_ext_edges.push_back( { ( *it ).first, p_inf_idx } );

  vin.is_boundary_tri.resize( n_tri );
  for( std::size_t i = 0; i < n_tri; ++i )
    vin.is_boundary_tri[ i ] = dual.is_boundary( i );
  vin.is_bolsillo = is_bolsillo;
  vin.pocket_component = pocket_component;
  vin.num_concavities  = contador;

  run_concavidades_viz( argv[ 0 ], vin );

  return( EXIT_SUCCESS );
}

// eof - concavidades.cxx
