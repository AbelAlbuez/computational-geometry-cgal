// =========================================================================
// @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__Polygon_Functions__hxx__
#define __pujCGAL__Polygon_Functions__hxx__

#include <list>
#include <stack>
#include <tuple>
#include <vector>

#include <CGAL/partition_2.h>

// -------------------------------------------------------------------------
template< class _TTriangulation, class _TPolygon >
void pujCGAL::
triangulate( _TTriangulation& mesh, const _TPolygon& polygon )
{
  // Compute a monotone partition
  std::list< typename _TPolygon::TTopology > mono_polys;
  CGAL::y_monotone_partition_2(
    polygon.topology_begin( ), polygon.topology_end( ),
    std::back_inserter( mono_polys ), polygon.traits( )
    );

  // Copy points
  mesh.clear( );
  for(
    auto pIt = polygon.geometry_begin( );
    pIt != polygon.geometry_end( );
    ++pIt
    )
    mesh.add_point( *pIt );

  // Compute triangles
  for( const auto& m: mono_polys )
    pujCGAL::triangulate( mesh, polygon, m );
}

// -------------------------------------------------------------------------
template< class _TTriangulation, class _TPolygon >
void pujCGAL::
triangulate(
  _TTriangulation& T,
  const _TPolygon& P, const typename _TPolygon::TTopology& M
  )
{
  if( M.size( ) > 3 )
  {
    using TPoint = typename _TPolygon::TPoint;
    using TIndex = typename _TPolygon::TIndex;
    using TSweepKey = std::tuple< TPoint, TIndex, bool >;
    using TSweep = std::vector< TSweepKey >;

    // Sort input
    TSweep sweep;
    for( auto mIt = M.vertices_begin( ); mIt != M.vertices_end( ); ++mIt )
      sweep.push_back( TSweepKey( P.point( *mIt ), *mIt, false ) );
    std::sort(
      sweep.begin( ), sweep.end( ),
      []( const TSweepKey& a, const TSweepKey& b ) -> bool
      {
        const TPoint& p = std::get< 0 >( a );
        const TPoint& q = std::get< 0 >( b );
        if( p[ 1 ] == q[ 1 ] ) return( p[ 0 ] < q[ 0 ] );
        else                   return( q[ 1 ] < p[ 1 ] );
      }
      );

    // Mark chains
    auto mark_left = []( TSweep& s, const TIndex& i ) -> void
      {
        auto sIt = s.begin( );
        bool stop = false;
        while( sIt != s.end( ) && !stop )
        {
          if( std::get< 1 >( *sIt ) == i )
          {
            std::get< 2 >( *sIt ) = true;
            stop = true;
          } // end if
          sIt++;
        } // end while
      };
    auto vIt = sweep.begin( );
    auto bIdx = std::find( M.begin( ), M.end( ), std::get< 1 >( *vIt ) );
    auto eIdx
      =
      std::find( M.begin( ), M.end( ), std::get< 1 >( *sweep.rbegin( ) ) );
    auto i = bIdx;
    while( i != eIdx )
    {
      mark_left( sweep, *i );
      i++;
      if( i == M.end( ) )
        i = M.begin( );
    } // end while
    mark_left( sweep, *i );

    // Initialize stack
    std::stack< TSweepKey > s;
    s.push( *( vIt++ ) );
    s.push( *( vIt++ ) );

    // Main loop
    for( ; std::next( vIt ) != sweep.end( ); ++vIt )
    {
      auto sv = s.top( );
      if( std::get< 2 >( *vIt ) == std::get< 2 >( sv ) ) // Same chain
      {
        bool valid = true;
        auto psv = sv;
        s.pop( );
        while( valid && s.size( ) > 0 )
        {
          auto ssv = s.top( );
          auto ori
            =
            CGAL::orientation(
              std::get< 0 >( *vIt ),
              std::get< 0 >( psv ),
              std::get< 0 >( ssv )
              );
          if(
            ( std::get< 2 >( *vIt ) && ( ori != CGAL::LEFT_TURN ) )
            ||
            ( !( std::get< 2 >( *vIt ) ) && ( ori != CGAL::RIGHT_TURN ) )
            )
          {
            T.add_triangle(
              std::get< 1 >( *vIt ), std::get< 1 >( ssv ), std::get< 1 >( psv )
              );
            s.pop( );
            psv = ssv;
          }
          else
            valid = false;
        } // end while

        s.push( psv );
        s.push( *vIt );
      }
      else // Different chain
      {
        while( s.size( ) > 1 )
        {
          auto st = s.top( );
          s.pop( );
          auto sb = s.top( );

          T.add_triangle(
            std::get< 1 >( *vIt ), std::get< 1 >( st ), std::get< 1 >( sb )
            );
        } // end while
        s.pop( );

        s.push( sv );
        s.push( *vIt );
      } // end if
    } // end for

    // Last vertex
      while( s.size( ) > 1 )
      {
        auto st = s.top( );
        s.pop( );
        auto sb = s.top( );

        T.add_triangle(
          std::get< 1 >( *vIt ), std::get< 1 >( st ), std::get< 1 >( sb )
          );
      } // end while
  }
  else if( M.size( ) == 3 )
    T.add_triangle( M[ 0 ], M[ 1 ], M[ 2 ] );
  else
    // do nothing
    ;
}

#endif // __pujCGAL__Polygon_Functions__hxx__

// eof - Polygon_Functions.hxx
