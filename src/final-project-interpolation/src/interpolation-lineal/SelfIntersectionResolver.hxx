#ifndef __SelfIntersectionResolver__hxx__
#define __SelfIntersectionResolver__hxx__

#include <algorithm>
#include <cmath>
#include <iterator>
#include <unordered_map>
#include <vector>

#include <CGAL/squared_distance_2.h>

#include <pujCGAL/SegmentsIntersection.h>

// -------------------------------------------------------------------------
namespace pujCGAL
{
  namespace Final
  {
    namespace _sir_detail
    {
      template< class TPoint >
      inline double signed_area( const std::vector< TPoint >& L )
      {
        const std::size_t n = L.size( );
        if( n < 3 ) return 0.0;
        double a = 0.0;
        for( std::size_t i = 0; i < n; ++i )
        {
          const TPoint& p = L[ i ];
          const TPoint& q = L[ ( i + 1 ) % n ];
          a += CGAL::to_double( p.x( ) ) * CGAL::to_double( q.y( ) )
             - CGAL::to_double( q.x( ) ) * CGAL::to_double( p.y( ) );
        }
        return 0.5 * a;
      }
    } // end namespace
  } // end namespace
} // end namespace

// -------------------------------------------------------------------------
inline bool
pujCGAL::Final::SelfIntersectionResolver::
has_self_intersections( const TContour& contour )
{
  const std::size_t n = contour.size( );
  if( n < 4 ) return false;

  std::vector< TSegment > segs;
  segs.reserve( n );
  for( std::size_t i = 0; i < n; ++i )
    segs.emplace_back( contour[ i ], contour[ ( i + 1 ) % n ] );

  std::vector< TPoint > inters;
  pujCGAL::SegmentsIntersection::BentleyOttmann(
    segs.begin( ), segs.end( ), std::back_inserter( inters )
    );

  // Filter out trivial "intersections" that are just the shared endpoint
  // between consecutive segments.
  for( const TPoint& p : inters )
  {
    bool shared_endpoint = false;
    for( std::size_t i = 0; i < n; ++i )
    {
      const TPoint& v = contour[ i ];
      const double dx = CGAL::to_double( v.x( ) ) - CGAL::to_double( p.x( ) );
      const double dy = CGAL::to_double( v.y( ) ) - CGAL::to_double( p.y( ) );
      if( dx * dx + dy * dy < 1e-18 )
      { shared_endpoint = true; break; }
    }
    if( !shared_endpoint ) return true;
  }
  return false;
}

// -------------------------------------------------------------------------
inline pujCGAL::Final::SelfIntersectionResolver::TContour
pujCGAL::Final::SelfIntersectionResolver::
resolve( const TContour& contour )
{
  const std::size_t n = contour.size( );
  if( n < 4 ) return contour;

  // -- 1) Build segments of the closed contour.
  std::vector< TSegment > segs;
  segs.reserve( n );
  for( std::size_t i = 0; i < n; ++i )
    segs.emplace_back( contour[ i ], contour[ ( i + 1 ) % n ] );

  // -- 2) Detect self-intersections.
  std::vector< TPoint > inters;
  pujCGAL::SegmentsIntersection::BentleyOttmann(
    segs.begin( ), segs.end( ), std::back_inserter( inters )
    );

  // Drop intersections that coincide with an original vertex (shared
  // endpoint between consecutive edges) -- these are not real crossings.
  auto is_original_vertex = [ & ]( const TPoint& p ) -> bool
  {
    for( std::size_t i = 0; i < n; ++i )
    {
      const double dx = CGAL::to_double( contour[ i ].x( ) )
                      - CGAL::to_double( p.x( ) );
      const double dy = CGAL::to_double( contour[ i ].y( ) )
                      - CGAL::to_double( p.y( ) );
      if( dx * dx + dy * dy < 1e-18 ) return true;
    }
    return false;
  };

  std::vector< TPoint > real_inters;
  real_inters.reserve( inters.size( ) );
  for( const TPoint& p : inters )
    if( !is_original_vertex( p ) )
      real_inters.push_back( p );

  if( real_inters.empty( ) ) return contour;

  // -- 3) For each segment, list the intersection points it carries with
  //       their parameter t (so we can sort along the segment direction).
  std::vector< std::vector< std::pair< double, int > > > per_seg( n );
  for( int j = 0; j < (int) real_inters.size( ); ++j )
  {
    const TPoint& p = real_inters[ j ];
    for( std::size_t i = 0; i < n; ++i )
    {
      const TPoint& a = segs[ i ].source( );
      const TPoint& b = segs[ i ].target( );
      const double ax = CGAL::to_double( a.x( ) );
      const double ay = CGAL::to_double( a.y( ) );
      const double bx = CGAL::to_double( b.x( ) );
      const double by = CGAL::to_double( b.y( ) );
      const double px = CGAL::to_double( p.x( ) );
      const double py = CGAL::to_double( p.y( ) );
      const double vx = bx - ax;
      const double vy = by - ay;
      const double len2 = vx * vx + vy * vy;
      if( len2 <= 0.0 ) continue;

      const double t = ( ( px - ax ) * vx + ( py - ay ) * vy ) / len2;
      if( t < -1e-9 || t > 1.0 + 1e-9 ) continue;

      const double tc = std::clamp( t, 0.0, 1.0 );
      const double cx = ax + tc * vx - px;
      const double cy = ay + tc * vy - py;
      if( cx * cx + cy * cy < 1e-12 )
        per_seg[ i ].emplace_back( tc, j );
    }
  }
  for( auto& v : per_seg )
    std::sort( v.begin( ), v.end( ) );

  // -- 4) Build the subdivided traversal: original vertex, then the
  //       intersection points along that edge, in arc-length order.
  struct Node { TPoint p; int tag; }; // tag = -1 if original, else index in real_inters
  std::vector< Node > path;
  path.reserve( n + 2 * real_inters.size( ) );
  for( std::size_t i = 0; i < n; ++i )
  {
    path.push_back( { contour[ i ], -1 } );
    for( const auto& tp : per_seg[ i ] )
      path.push_back( { real_inters[ tp.second ], tp.second } );
  }

  // -- 5) Stack-based loop extraction. When the same intersection tag is
  //       seen for the second time we close a loop with everything pushed
  //       since the first occurrence.
  std::vector< Node > stack;
  std::unordered_map< int, int > tag_pos;
  std::vector< std::vector< TPoint > > loops;

  for( const Node& nd : path )
  {
    if( nd.tag >= 0 )
    {
      auto it = tag_pos.find( nd.tag );
      if( it != tag_pos.end( ) )
      {
        const int idx = it->second;
        std::vector< TPoint > loop;
        loop.reserve( stack.size( ) - idx );
        for( std::size_t k = idx; k < stack.size( ); ++k )
          loop.push_back( stack[ k ].p );
        loops.push_back( std::move( loop ) );

        // Pop the loop; keep the meeting point so traversal continues.
        for( std::size_t k = idx + 1; k < stack.size( ); ++k )
          if( stack[ k ].tag >= 0 )
            tag_pos.erase( stack[ k ].tag );
        stack.resize( idx + 1 );
        continue;
      }
      tag_pos[ nd.tag ] = (int) stack.size( );
    }
    stack.push_back( nd );
  }

  // The remainder of the stack closes back through the start of the path.
  if( stack.size( ) >= 3 )
  {
    std::vector< TPoint > loop;
    loop.reserve( stack.size( ) );
    for( const Node& s : stack )
      loop.push_back( s.p );
    loops.push_back( std::move( loop ) );
  }

  if( loops.empty( ) ) return contour;

  // -- 6) Return the loop with the largest absolute area.
  std::size_t best = 0;
  double best_a = std::abs( _sir_detail::signed_area( loops[ 0 ] ) );
  for( std::size_t i = 1; i < loops.size( ); ++i )
  {
    const double a = std::abs( _sir_detail::signed_area( loops[ i ] ) );
    if( a > best_a ) { best_a = a; best = i; }
  }
  return loops[ best ];
}

#endif // __SelfIntersectionResolver__hxx__

// eof - SelfIntersectionResolver.hxx
