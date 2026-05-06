#ifndef __pujCGAL__Algoritmos__hxx__
#define __pujCGAL__Algoritmos__hxx__

#include <algorithm>
#include <chrono>
#include <cmath>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <CGAL/number_utils.h>

namespace pujCGAL
{
  namespace detail
  {
    // -----------------------------------------------------------------------
    // Helpers privados de los algoritmos.
    // -----------------------------------------------------------------------

    /**
     * Interpolacion lineal (baricentrica) de altura sobre una cara
     * de la triangulacion. Asume que cada vertice de la cara tiene
     * info() = altura. Si la cara es degenerada (area ~ 0) devuelve
     * el promedio simple.
     */
    template< class TFaceHandle >
    inline double interpolar_altura_baricentrica(
        TFaceHandle f, double x, double y
        )
    {
      const auto& p0 = f->vertex( 0 )->point( );
      const auto& p1 = f->vertex( 1 )->point( );
      const auto& p2 = f->vertex( 2 )->point( );
      const double x0 = CGAL::to_double( p0.x( ) );
      const double y0 = CGAL::to_double( p0.y( ) );
      const double x1 = CGAL::to_double( p1.x( ) );
      const double y1 = CGAL::to_double( p1.y( ) );
      const double x2 = CGAL::to_double( p2.x( ) );
      const double y2 = CGAL::to_double( p2.y( ) );
      const double h0 = static_cast< double >( f->vertex( 0 )->info( ) );
      const double h1 = static_cast< double >( f->vertex( 1 )->info( ) );
      const double h2 = static_cast< double >( f->vertex( 2 )->info( ) );

      const double denom = ( y1 - y2 ) * ( x0 - x2 )
                         + ( x2 - x1 ) * ( y0 - y2 );
      if( std::abs( denom ) < 1e-15 )
        return( ( h0 + h1 + h2 ) / 3.0 );

      const double a = ( ( y1 - y2 ) * ( x  - x2 )
                       + ( x2 - x1 ) * ( y  - y2 ) ) / denom;
      const double b = ( ( y2 - y0 ) * ( x  - x2 )
                       + ( x0 - x2 ) * ( y  - y2 ) ) / denom;
      const double c = 1.0 - a - b;
      return( a * h0 + b * h1 + c * h2 );
    }

    /**
     * Devuelve true si los tres vertices de la cara coinciden con la
     * posicion del pixel (w, h) del heightmap (en coordenadas reales).
     * Util para evitar registrar pixeles que ya estan representados
     * exactamente por un vertice (caso de las esquinas).
     */
    template< class TFaceHandle, class THeightmap >
    inline bool pixel_es_vertice(
        TFaceHandle f, std::size_t w, std::size_t h, const THeightmap& hm
        )
    {
      const auto p = hm.point( w, h );
      for( int i = 0; i < 3; ++i )
      {
        const auto& vp = f->vertex( i )->point( );
        if( CGAL::to_double( vp.x( ) ) == p.first &&
            CGAL::to_double( vp.y( ) ) == p.second )
          return( true );
      } // end for
      return( false );
    }

  } // end namespace detail


  // =========================================================================
  // 1. Submuestreo uniforme
  // =========================================================================
  template< class TDelaunay, class THeightmap >
  ResultadoSimplificacion simplificar_submuestreo_uniforme(
      TDelaunay& T,
      const THeightmap& hm,
      int paso
      )
  {
    using TPoint = typename TDelaunay::Point;

    ResultadoSimplificacion r;
    r.nombre_algoritmo
      = "submuestreo_uniforme(paso=" + std::to_string( paso ) + ")";
    r.vertices_antes = hm.width( ) * hm.height( );

    if( paso < 1 ) paso = 1;

    const auto t0 = std::chrono::steady_clock::now( );

    // Submuestrear: incluir las 4 esquinas explicitamente para garantizar
    // un dominio bien formado, luego una rejilla regular cada `paso` pixeles.
    std::vector< TPoint > puntos;
    const std::size_t W = hm.width( );
    const std::size_t H = hm.height( );
    puntos.reserve( ( W / paso + 2 ) * ( H / paso + 2 ) );

    auto agregar = [&] ( std::size_t w, std::size_t h )
    {
      const auto p = hm.point( w, h );
      puntos.push_back( TPoint( p.first, p.second ) );
    };
    agregar( 0,     0     );
    agregar( W - 1, 0     );
    agregar( 0,     H - 1 );
    agregar( W - 1, H - 1 );
    for( std::size_t h = 0; h < H; h += paso )
      for( std::size_t w = 0; w < W; w += paso )
        agregar( w, h );

    T.insert( puntos.begin( ), puntos.end( ) );

    // Asignar info() = altura a cada vertice.
    for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
    {
      const auto idx = hm.index( v->point( ).x( ), v->point( ).y( ) );
      v->info( ) = hm( idx.first, idx.second );
    } // end for

    const auto t1 = std::chrono::steady_clock::now( );
    r.tiempo_milisegundos
      = std::chrono::duration< double, std::milli >( t1 - t0 ).count( );
    r.vertices_despues   = T.number_of_vertices( );
    r.triangulos_despues = T.number_of_faces( );
    r.reduccion_porcentaje
      = 100.0 * ( 1.0 - static_cast< double >( r.vertices_despues )
                       / static_cast< double >( r.vertices_antes  ) );
    return( r );
  }


  // =========================================================================
  // 2. Decimacion con error L1-promedio (algoritmo del taller)
  // =========================================================================
  template< class TDelaunay, class THeightmap >
  ResultadoSimplificacion simplificar_decimacion_L1(
      TDelaunay& T,
      const THeightmap& hm,
      double epsilon,
      int orden_k
      )
  {
    using TPoint          = typename TDelaunay::Point;
    using TVertexHandle   = typename TDelaunay::Vertex_handle;
    using TFaceCirculator = typename TDelaunay::Face_circulator;

    ResultadoSimplificacion r;
    r.nombre_algoritmo
      = "decimacion_L1(eps=" + std::to_string( epsilon )
      + ",k="                + std::to_string( orden_k ) + ")";
    r.vertices_antes = hm.width( ) * hm.height( );

    const auto t0 = std::chrono::steady_clock::now( );

    // --- Fase 1: insertar todos los pixeles y asignar alturas.
    {
      std::vector< TPoint > puntos;
      puntos.reserve( hm.width( ) * hm.height( ) );
      for( std::size_t h = 0; h < hm.height( ); ++h )
        for( std::size_t w = 0; w < hm.width( ); ++w )
        {
          const auto p = hm.point( w, h );
          puntos.push_back( TPoint( p.first, p.second ) );
        } // end for
      T.insert( puntos.begin( ), puntos.end( ) );
      for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
      {
        const auto idx = hm.index( v->point( ).x( ), v->point( ).y( ) );
        v->info( ) = hm( idx.first, idx.second );
      } // end for
    }

    // Lambda: error de planitud L1 (promedio de |dz| con vecinos en la estrella).
    auto compute_error = [ & ] ( TVertexHandle v ) -> double
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
      return( cnt > 0 ? sum / cnt : 0.0 );
    };

    // --- Fase 2: construir min-heap inicial con todos los vertices.
    using THeapEntry = std::pair< double, TVertexHandle >;
    std::priority_queue<
      THeapEntry,
      std::vector< THeapEntry >,
      std::greater< THeapEntry >
    > heap;

    std::unordered_map< TVertexHandle, double > error_map;
    for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
    {
      const double err = compute_error( v );
      error_map[ v ] = err;
      heap.push( { err, v } );
    } // end for

    std::unordered_set< TVertexHandle > visited;

    // --- Fase 3: bucle principal.
    while( !heap.empty( ) )
    {
      auto [ err, p ] = heap.top( );
      heap.pop( );

      if( err >= epsilon ) continue;          // por encima del umbral
      if( visited.count( p ) ) continue;      // ya procesado (lazy delete)
      auto it = error_map.find( p );
      if( it == error_map.end( ) ) continue;  // entrada obsoleta
      if( std::abs( it->second - err ) > 1e-9 ) continue;

      visited.insert( p );

      // Recolectar vecindad (anillo 1) por face circulator (DCEL).
      std::set< TVertexHandle > neighbors_k;
      {
        TFaceCirculator fc = T.incident_faces( p ), done = fc;
        do
        {
          if( !T.is_infinite( fc ) )
            for( int i = 0; i < 3; ++i )
            {
              TVertexHandle q = fc->vertex( i );
              if( q != p && !T.is_infinite( q ) )
                neighbors_k.insert( q );
            } // end for
          ++fc;
        } while( fc != done );
      }

      // Expandir hasta orden k.
      for( int ring = 1; ring < orden_k; ++ring )
      {
        std::set< TVertexHandle > next_ring;
        for( auto qi : neighbors_k )
        {
          TFaceCirculator fc2 = T.incident_faces( qi ), done2 = fc2;
          do
          {
            if( !T.is_infinite( fc2 ) )
              for( int i = 0; i < 3; ++i )
              {
                TVertexHandle q = fc2->vertex( i );
                if( q != p && !T.is_infinite( q ) && !neighbors_k.count( q ) )
                  next_ring.insert( q );
              } // end for
            ++fc2;
          } while( fc2 != done2 );
        } // end for
        neighbors_k.insert( next_ring.begin( ), next_ring.end( ) );
      } // end for

      // Eliminar p; CGAL retriangula (Delaunay local).
      // Mantener al menos 4 vertices para que la triangulacion 2D sea valida.
      if( T.number_of_vertices( ) <= 4 ) break;
      error_map.erase( p );
      T.remove( p );

      // Recalcular errores de los vecinos afectados y reinsertar.
      for( auto qi : neighbors_k )
      {
        if( visited.count( qi ) ) continue;
        if( error_map.find( qi ) == error_map.end( ) ) continue;
        const double new_err = compute_error( qi );
        error_map[ qi ] = new_err;
        heap.push( { new_err, qi } );
      } // end for
    } // end while

    const auto t1 = std::chrono::steady_clock::now( );
    r.tiempo_milisegundos
      = std::chrono::duration< double, std::milli >( t1 - t0 ).count( );
    r.vertices_despues   = T.number_of_vertices( );
    r.triangulos_despues = T.number_of_faces( );
    r.reduccion_porcentaje
      = 100.0 * ( 1.0 - static_cast< double >( r.vertices_despues )
                       / static_cast< double >( r.vertices_antes  ) );
    return( r );
  }


  // =========================================================================
  // 3. Insercion golosa de Garland-Heckbert (1995)
  // =========================================================================
  template< class TDelaunay, class THeightmap >
  ResultadoSimplificacion simplificar_insercion_greedy(
      TDelaunay& T,
      const THeightmap& hm,
      double epsilon_max,
      std::size_t max_vertices
      )
  {
    using TPoint        = typename TDelaunay::Point;
    using TFaceHandle   = typename TDelaunay::Face_handle;
    using TVertexHandle = typename TDelaunay::Vertex_handle;

    ResultadoSimplificacion r;
    r.nombre_algoritmo
      = "insercion_greedy(eps=" + std::to_string( epsilon_max ) + ")";
    r.vertices_antes = hm.width( ) * hm.height( );

    const auto t0 = std::chrono::steady_clock::now( );

    // Pixel candidato dentro del bucket de una cara.
    struct PixelCandidato
    {
      std::size_t w;
      std::size_t h;
      double      altura_real;
    };

    // --- Fase 1: insertar las 4 esquinas como semilla.
    auto insertar_esquina = [ & ] ( std::size_t w, std::size_t h ) -> TVertexHandle
    {
      const auto p = hm.point( w, h );
      auto v = T.insert( TPoint( p.first, p.second ) );
      v->info( ) = hm( w, h );
      return( v );
    };
    insertar_esquina( 0,                0                );
    insertar_esquina( hm.width( ) - 1,  0                );
    insertar_esquina( 0,                hm.height( ) - 1 );
    insertar_esquina( hm.width( ) - 1,  hm.height( ) - 1 );

    // Por cada cara: lista de pixeles candidatos y peor (error, pixel).
    std::unordered_map< TFaceHandle, std::vector< PixelCandidato > > face_buckets;
    std::unordered_map<
      TFaceHandle, std::pair< double, PixelCandidato >
    > face_worst;

    // Recalcula el peor pixel de una cara dada y actualiza face_worst.
    auto recompute_face_worst = [ & ] ( TFaceHandle f )
    {
      auto bit = face_buckets.find( f );
      if( bit == face_buckets.end( ) || bit->second.empty( ) )
      {
        face_worst.erase( f );
        return;
      } // end if
      double best_err = -1.0;
      PixelCandidato best_px { 0, 0, 0.0 };
      for( const auto& px : bit->second )
      {
        const auto p = hm.point( px.w, px.h );
        const double interp
          = detail::interpolar_altura_baricentrica( f, p.first, p.second );
        const double err = std::abs( px.altura_real - interp );
        if( err > best_err )
        {
          best_err = err;
          best_px  = px;
        } // end if
      } // end for
      face_worst[ f ] = { best_err, best_px };
    };

    // --- Fase 2: bucketing inicial. Para cada pixel: localizar la cara
    //     finita que lo contiene y agregarlo a su bucket.
    {
      TFaceHandle hint;
      for( std::size_t hh = 0; hh < hm.height( ); ++hh )
        for( std::size_t ww = 0; ww < hm.width( ); ++ww )
        {
          const auto p = hm.point( ww, hh );
          TPoint qp( p.first, p.second );
          TFaceHandle f = T.locate( qp, hint );
          hint = f;
          if( T.is_infinite( f ) )                    continue;
          if( detail::pixel_es_vertice( f, ww, hh, hm ) ) continue;
          face_buckets[ f ].push_back(
            { ww, hh, static_cast< double >( hm( ww, hh ) ) }
            );
        } // end for
      for( const auto& kv : face_buckets )
        recompute_face_worst( kv.first );
    }

    // --- Fase 3: bucle goloso con max-heap (default priority_queue).
    using THeapEntry = std::pair< double, TFaceHandle >;
    std::priority_queue< THeapEntry > heap;
    for( const auto& kv : face_worst )
      heap.push( { kv.second.first, kv.first } );

    while( !heap.empty( ) )
    {
      // Cota dura por numero de vertices.
      if( max_vertices > 0 && T.number_of_vertices( ) >= max_vertices )
        break;

      auto [ err, f ] = heap.top( );
      heap.pop( );

      // Lazy delete: la entrada puede estar obsoleta.
      auto it = face_worst.find( f );
      if( it == face_worst.end( ) )                         continue;
      if( std::abs( it->second.first - err ) > 1e-9 )       continue;
      if( err <= epsilon_max )                              break;

      const PixelCandidato peor = it->second.second;
      const auto p = hm.point( peor.w, peor.h );
      const TPoint qp( p.first, p.second );

      // Antes de insertar: recolectar caras del conflict region (Bowyer-Watson).
      std::vector< TFaceHandle > conflict_faces;
      T.get_conflicts( qp, std::back_inserter( conflict_faces ), f );

      // Recolectar pixeles huerfanos y limpiar tablas.
      std::vector< PixelCandidato > orfanos;
      for( auto cf : conflict_faces )
      {
        auto bit = face_buckets.find( cf );
        if( bit != face_buckets.end( ) )
        {
          for( const auto& px : bit->second )
            if( !( px.w == peor.w && px.h == peor.h ) )
              orfanos.push_back( px );
          face_buckets.erase( bit );
        } // end if
        face_worst.erase( cf );
      } // end for

      // Insertar el peor pixel como nuevo vertice.
      auto v = T.insert( qp );
      v->info( ) = peor.altura_real;

      // Re-bucketing de huerfanos en las caras nuevas.
      TFaceHandle hint = v->face( );
      for( const auto& px : orfanos )
      {
        const auto pp = hm.point( px.w, px.h );
        TFaceHandle nf = T.locate( TPoint( pp.first, pp.second ), hint );
        hint = nf;
        if( T.is_infinite( nf ) )                            continue;
        if( detail::pixel_es_vertice( nf, px.w, px.h, hm ) ) continue;
        face_buckets[ nf ].push_back( px );
      } // end for

      // Recomputar peor en las caras incidentes a v y empujarlas al heap.
      auto fc = T.incident_faces( v ), fc_done = fc;
      do
      {
        if( !T.is_infinite( fc ) )
        {
          recompute_face_worst( fc );
          auto fwit = face_worst.find( fc );
          if( fwit != face_worst.end( ) )
            heap.push( { fwit->second.first, fc } );
        } // end if
        ++fc;
      } while( fc != fc_done );
    } // end while

    const auto t1 = std::chrono::steady_clock::now( );
    r.tiempo_milisegundos
      = std::chrono::duration< double, std::milli >( t1 - t0 ).count( );
    r.vertices_despues   = T.number_of_vertices( );
    r.triangulos_despues = T.number_of_faces( );
    r.reduccion_porcentaje
      = 100.0 * ( 1.0 - static_cast< double >( r.vertices_despues )
                       / static_cast< double >( r.vertices_antes  ) );
    return( r );
  }


  // =========================================================================
  // Calculo de metricas de calidad
  // =========================================================================
  template< class TDelaunay, class THeightmap >
  void calcular_metricas(
      ResultadoSimplificacion& resultado,
      const TDelaunay& T,
      const THeightmap& hm
      )
  {
    using TPoint      = typename TDelaunay::Point;
    using TFaceHandle = typename TDelaunay::Face_handle;

    // Errores verticales: recorrer todos los pixeles y comparar con interp.
    double      max_err  = 0.0;
    double      sum_sq   = 0.0;
    double      sum_abs  = 0.0;
    std::size_t n        = 0;

    TFaceHandle hint;
    for( std::size_t h = 0; h < hm.height( ); ++h )
      for( std::size_t w = 0; w < hm.width( ); ++w )
      {
        const auto p = hm.point( w, h );
        TFaceHandle f = T.locate( TPoint( p.first, p.second ), hint );
        hint = f;
        if( T.is_infinite( f ) ) continue;

        const double interp
          = detail::interpolar_altura_baricentrica( f, p.first, p.second );
        const double err
          = std::abs( static_cast< double >( hm( w, h ) ) - interp );
        max_err = std::max( max_err, err );
        sum_sq  += err * err;
        sum_abs += err;
        ++n;
      } // end for

    resultado.error_L_infinito = max_err;
    resultado.error_L2_RMSE
      = ( n > 0 ) ? std::sqrt( sum_sq / static_cast< double >( n ) ) : 0.0;
    resultado.error_promedio
      = ( n > 0 ) ? sum_abs / static_cast< double >( n )             : 0.0;

    // Angulos: recorrer todas las caras finitas.
    constexpr double PI = 3.14159265358979323846;
    double ang_min = 180.0;
    double ang_max = 0.0;

    auto angulo_en_b = [ ] (
        double ax, double ay,
        double bx, double by,
        double cx, double cy
        ) -> double
    {
      const double ux = ax - bx, uy = ay - by;
      const double vx = cx - bx, vy = cy - by;
      const double dot = ux * vx + uy * vy;
      const double nu  = std::sqrt( ux * ux + uy * uy );
      const double nv  = std::sqrt( vx * vx + vy * vy );
      if( nu < 1e-15 || nv < 1e-15 ) return( 0.0 );
      double cos_t = dot / ( nu * nv );
      cos_t = std::clamp( cos_t, -1.0, 1.0 );
      return( std::acos( cos_t ) * 180.0 / PI );
    };

    for( auto f = T.finite_faces_begin( ); f != T.finite_faces_end( ); ++f )
    {
      const auto& p0 = f->vertex( 0 )->point( );
      const auto& p1 = f->vertex( 1 )->point( );
      const auto& p2 = f->vertex( 2 )->point( );
      const double x0 = CGAL::to_double( p0.x( ) );
      const double y0 = CGAL::to_double( p0.y( ) );
      const double x1 = CGAL::to_double( p1.x( ) );
      const double y1 = CGAL::to_double( p1.y( ) );
      const double x2 = CGAL::to_double( p2.x( ) );
      const double y2 = CGAL::to_double( p2.y( ) );

      const double a0 = angulo_en_b( x1, y1, x0, y0, x2, y2 );
      const double a1 = angulo_en_b( x0, y0, x1, y1, x2, y2 );
      const double a2 = std::max( 0.0, 180.0 - a0 - a1 );
      ang_min = std::min( { ang_min, a0, a1, a2 } );
      ang_max = std::max( { ang_max, a0, a1, a2 } );
    } // end for

    resultado.angulo_minimo_grados = ang_min;
    resultado.angulo_maximo_grados = ang_max;
  }

} // end namespace pujCGAL

#endif // __pujCGAL__Algoritmos__hxx__

// eof - Algoritmos.hxx
