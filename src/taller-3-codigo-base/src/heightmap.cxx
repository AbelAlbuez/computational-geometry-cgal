
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
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
  using TVertexHandle = TDelaunay::Vertex_handle;
  using TFaceCirculator = TDelaunay::Face_circulator;

  if( argc < 2 )
  {
    std::cerr << "Uso: " << argv[ 0 ]
              << " input.png [epsilon] [orden_k]" << std::endl;
    return( EXIT_FAILURE );
  } // end if

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

  // Construir la ruta de output/ un nivel arriba del ejecutable (build/ -> src/)
  std::filesystem::path exe_path = std::filesystem::canonical( argv[ 0 ] );
  std::filesystem::path output_dir = exe_path.parent_path( ).parent_path( ) / "output";
  std::filesystem::create_directories( output_dir );

  // IO::save usa v->info() como indice OBJ. Guardamos y restauramos alturas
  // para no alterar el criterio de planitud del algoritmo.
  std::unordered_map< TVertexHandle, TReal > original_heights;
  for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
    original_heights[ v ] = v->info( );

  pujCGAL::IO::save( T, ( output_dir / "original.obj" ).string( ) );
  for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
    v->info( ) = original_heights[ v ];

  // =========================================================================
  // Taller 3: Simplificación por vecindad de orden k
  //
  // Concepto geométrico:
  //   Un vértice es "redundante" si su altura no difiere significativamente de
  //   sus vecinos en la triangulación (zona plana del heightmap). El error de
  //   planitud mide esa diferencia promedio. Usamos un min-heap (Dijkstra-like)
  //   para procesar primero los vértices con menor error, eliminándolos si el
  //   error < epsilon. CGAL re-triangula el hueco automáticamente (Delaunay).
  //
  //   La vecindad de orden k expande el anillo de vecinos k veces para tener
  //   en cuenta un contexto más amplio al recalcular errores tras una eliminación.
  // =========================================================================

  // Parámetros: epsilon (umbral de planitud) y orden_k (radio de vecindad)
  const double epsilon = ( argc > 2 ) ? std::stod( argv[ 2 ] ) : 10.0;
  const int    orden_k = ( argc > 3 ) ? std::stoi( argv[ 3 ] ) : 2;

  // --- Función lambda: calcula el error de planitud de un vértice v
  //     usando la estrella (Face_circulator ~ half-edge twin→next).
  //     error(v) = promedio de |v->info() - qi->info()| para cada vecino qi
  auto compute_error = [&]( TVertexHandle v ) -> double
  {
    double sum = 0.0;
    int    cnt = 0;
    TFaceCirculator fc = T.incident_faces( v ), done = fc;
    do
    {
      if( !T.is_infinite( fc ) )
      {
        for( int i = 0; i < 3; ++i )
        {
          TVertexHandle q = fc->vertex( i );
          if( q != v && !T.is_infinite( q ) )
          {
            sum += std::abs( static_cast< double >( v->info( ) )
                           - static_cast< double >( q->info( ) ) );
            ++cnt;
          } // end if
        } // end for
      } // end if
      ++fc;
    } while( fc != done );
    return cnt > 0 ? sum / cnt : 0.0;
  };

  // --- Paso 1: construir el min-heap inicial con todos los vértices finitos
  using THeapEntry = std::pair< double, TVertexHandle >;
  std::priority_queue<
    THeapEntry,
    std::vector< THeapEntry >,
    std::greater< THeapEntry >
  > heap;

  std::unordered_map< TVertexHandle, double > error_map;
  for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
  {
    double err = compute_error( v );
    error_map[ v ] = err;
    heap.push( { err, v } );
  } // end for

  std::unordered_set< TVertexHandle > visited;
  const std::size_t before = T.number_of_vertices( );

  // --- Pasos 2–7: procesamiento del heap
  while( !heap.empty( ) )
  {
    auto [err, p] = heap.top( );
    heap.pop( );

    // Caso límite: entrada con error sobre umbral, pero continuar por lazy heap
    if( err >= epsilon ) continue;

    // Saltar si ya fue procesado (lazy deletion del heap)
    if( visited.count( p ) ) continue;

    // Saltar si el error en el mapa ya no coincide (entrada obsoleta)
    auto it = error_map.find( p );
    if( it == error_map.end( ) ) continue;
    if( std::abs( it->second - err ) > 1e-9 ) continue;

    visited.insert( p );

    // --- Paso 3: recorrer la estrella de p con el Face_circulator (DCEL)
    //     Recolectar vecinos directos (anillo de orden 1)
    std::set< TVertexHandle > neighbors_k;
    {
      TFaceCirculator fc = T.incident_faces( p ), done = fc;
      do
      {
        if( !T.is_infinite( fc ) )
        {
          for( int i = 0; i < 3; ++i )
          {
            TVertexHandle q = fc->vertex( i );
            if( q != p && !T.is_infinite( q ) )
              neighbors_k.insert( q );
          } // end for
        } // end if
        ++fc;
      } while( fc != done );
    }

    // --- Paso 4: expandir la vecindad hasta orden k
    //     Cada iteración agrega el siguiente anillo de vecinos
    for( int ring = 1; ring < orden_k; ++ring )
    {
      std::set< TVertexHandle > next_ring;
      for( auto qi : neighbors_k )
      {
        TFaceCirculator fc2 = T.incident_faces( qi ), done2 = fc2;
        do
        {
          if( !T.is_infinite( fc2 ) )
          {
            for( int i = 0; i < 3; ++i )
            {
              TVertexHandle q = fc2->vertex( i );
              if( q != p && !T.is_infinite( q ) && !neighbors_k.count( q ) )
                next_ring.insert( q );
            } // end for
          } // end if
          ++fc2;
        } while( fc2 != done2 );
      } // end for
      neighbors_k.insert( next_ring.begin( ), next_ring.end( ) );
    } // end for

    // --- Paso 5: eliminar el vértice redundante
    //     CGAL re-triangula automáticamente el hueco (Delaunay local).
    //     Caso límite: mantener al menos 4 vértices finitos para preservar
    //     la validez de la triangulación 2D (mínimo una cara finita).
    if( T.number_of_vertices( ) <= 4 )
      break;
    error_map.erase( p );
    T.remove( p );

    // --- Paso 6: recalcular el error de los vecinos afectados
    //     y reinsertarlos en el heap (lazy update)
    for( auto qi : neighbors_k )
    {
      if( visited.count( qi ) ) continue;
      if( error_map.find( qi ) == error_map.end( ) ) continue;
      double new_err = compute_error( qi );
      error_map[ qi ] = new_err;
      heap.push( { new_err, qi } );
    } // end for
  } // end while

  const std::size_t after = T.number_of_vertices( );
  const double reduction = 100.0 * ( 1.0 - static_cast< double >( after )
                                          / static_cast< double >( before ) );
  std::cout << "Vertices antes:   " << before  << std::endl;
  std::cout << "Vertices despues: " << after   << std::endl;
  std::cout << "Reduccion:        " << std::fixed << std::setprecision( 1 )
            << reduction << "%" << std::endl;

  // =========================================================================

  pujCGAL::IO::save( T, ( output_dir / "simplificado.obj" ).string( ) );

  // Llamar al visualizer con ruta absoluta
  std::string visualizer_cmd =
    "python3 \""
    + ( exe_path.parent_path( ).parent_path( ) / "visualizer.py" ).string( )
    + "\" --output \""
    + output_dir.string( )
    + "\"";
  std::system( visualizer_cmd.c_str( ) );
    
  return( EXIT_SUCCESS );
}

// eof - heightmap.cxx
