#ifndef __DistanceField__hxx__
#define __DistanceField__hxx__

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/Polygon_2.h>

namespace pujCGAL
{
  namespace Final
  {
    namespace df_detail
    {
      using K        = DistanceField::TKernel;
      using TPoint   = DistanceField::TPoint;       // Point_2
      using TContour = DistanceField::TContour;
      using Point3   = K::Point_3;
      using Segment3 = K::Segment_3;
      using Polygon2 = CGAL::Polygon_2< K >;

      // AABB tree over 2D segments embedded at z = 0 (CGAL's segment
      // primitive is 3D; with z=0 the 3D closest distance equals the 2D one).
      using SegIt      = std::vector< Segment3 >::const_iterator;
      using Primitive  = CGAL::AABB_segment_primitive< K, SegIt >;
      using AABBTraits = CGAL::AABB_traits< K, Primitive >;
      using Tree       = CGAL::AABB_tree< AABBTraits >;

      // -- small vector helpers -------------------------------------------
      inline TPoint rotate( const TPoint& p, double ang )
      {
        const double c = std::cos( ang ), s = std::sin( ang );
        return TPoint( c * p.x( ) - s * p.y( ), s * p.x( ) + c * p.y( ) );
      }

      inline TPoint centroid_loops( const std::vector< TContour >& loops )
      {
        double sx = 0, sy = 0; std::size_t n = 0;
        for( const auto& L : loops )
          for( const auto& p : L ) { sx += p.x( ); sy += p.y( ); ++n; }
        if( n == 0 ) return TPoint( 0, 0 );
        return TPoint( sx / n, sy / n );
      }

      // Principal-axis angle from the 2x2 point covariance (PCA). For an
      // isotropic shape (e.g. a circle) this returns a stable 0.
      inline double principal_angle( const std::vector< TContour >& loops,
                                     const TPoint& c )
      {
        double sxx = 0, syy = 0, sxy = 0; std::size_t n = 0;
        for( const auto& L : loops )
          for( const auto& p : L )
          {
            const double dx = p.x( ) - c.x( ), dy = p.y( ) - c.y( );
            sxx += dx * dx; syy += dy * dy; sxy += dx * dy; ++n;
          }
        if( n == 0 ) return 0.0;
        if( std::fabs( sxy ) < 1e-12 && std::fabs( sxx - syy ) < 1e-12 )
          return 0.0;                                  // isotropic -> canonical 0
        return 0.5 * std::atan2( 2.0 * sxy, sxx - syy );
      }

      inline double loop_signed_area( const TContour& c )
      {
        const std::size_t n = c.size( );
        if( n < 3 ) return 0.0;
        double a = 0.0;
        for( std::size_t i = 0; i < n; ++i )
          a += c[ i ].x( ) * c[ ( i + 1 ) % n ].y( )
             - c[ ( i + 1 ) % n ].x( ) * c[ i ].y( );
        return 0.5 * a;
      }

      // -- regular grid holding a sampled scalar field --------------------
      struct Grid
      {
        double x0 = 0, y0 = 0, h = 1.0;
        int    nx = 0, ny = 0;                 // number of cells per axis
        std::vector< double > phi;             // (nx+1)*(ny+1) node values
        double  val( int i, int j ) const { return phi[ (std::size_t) i * ( ny + 1 ) + j ]; }
        double& val( int i, int j )       { return phi[ (std::size_t) i * ( ny + 1 ) + j ]; }
        TPoint  pos( int i, int j ) const { return TPoint( x0 + i * h, y0 + j * h ); }
      };

      // Signed distance field of a set of loops sampled on grid g
      // (positive inside any loop). Distance from a CGAL::AABB_tree of the
      // contour segments; sign from CGAL::Polygon_2::bounded_side.
      inline void build_field( const std::vector< TContour >& loops, Grid& g )
      {
        std::vector< Segment3 > segs;
        for( const auto& L : loops )
        {
          const std::size_t n = L.size( );
          for( std::size_t i = 0; i < n; ++i )
            segs.emplace_back( Point3( L[ i ].x( ), L[ i ].y( ), 0 ),
                               Point3( L[ ( i + 1 ) % n ].x( ),
                                       L[ ( i + 1 ) % n ].y( ), 0 ) );
        }
        Tree tree( segs.begin( ), segs.end( ) );
        tree.accelerate_distance_queries( );

        std::vector< Polygon2 > polys;
        for( const auto& L : loops )
        {
          Polygon2 p;
          for( const auto& pt : L ) p.push_back( pt );
          polys.push_back( p );
        }

        g.phi.assign( (std::size_t)( g.nx + 1 ) * ( g.ny + 1 ), 0.0 );
        for( int i = 0; i <= g.nx; ++i )
          for( int j = 0; j <= g.ny; ++j )
          {
            const TPoint q = g.pos( i, j );
            const double d = std::sqrt(
                tree.squared_distance( Point3( q.x( ), q.y( ), 0 ) ) );
            bool inside = false;
            for( const auto& p : polys )
              if( p.bounded_side( q ) == CGAL::ON_BOUNDED_SIDE ) { inside = true; break; }
            g.val( i, j ) = inside ? d : -d;
          }
      }

      // Marching squares: extract every closed loop of the zero level set of
      // g. Crossings are indexed by the grid edge they lie on so neighboring
      // cells share endpoints exactly -> clean closed polylines.
      inline std::vector< TContour > extract_loops( const Grid& g )
      {
        const int nx = g.nx, ny = g.ny;
        const long long HBASE = (long long) nx * ( ny + 1 );
        auto Hid = [&]( int i, int j ) -> long long { return (long long) j * nx + i; };
        auto Vid = [&]( int i, int j ) -> long long { return HBASE + (long long) i * ny + j; };

        std::unordered_map< long long, TPoint > cross;
        auto getH = [&]( int i, int j ) -> long long {
          const long long id = Hid( i, j );
          auto it = cross.find( id );
          if( it == cross.end( ) )
          {
            const double v0 = g.val( i, j ), v1 = g.val( i + 1, j );
            const double tt = ( v0 == v1 ) ? 0.5 : v0 / ( v0 - v1 );
            const TPoint a = g.pos( i, j ), b = g.pos( i + 1, j );
            cross[ id ] = TPoint( a.x( ) + tt * ( b.x( ) - a.x( ) ),
                                  a.y( ) + tt * ( b.y( ) - a.y( ) ) );
          }
          return id;
        };
        auto getV = [&]( int i, int j ) -> long long {
          const long long id = Vid( i, j );
          auto it = cross.find( id );
          if( it == cross.end( ) )
          {
            const double v0 = g.val( i, j ), v1 = g.val( i, j + 1 );
            const double tt = ( v0 == v1 ) ? 0.5 : v0 / ( v0 - v1 );
            const TPoint a = g.pos( i, j ), b = g.pos( i, j + 1 );
            cross[ id ] = TPoint( a.x( ) + tt * ( b.x( ) - a.x( ) ),
                                  a.y( ) + tt * ( b.y( ) - a.y( ) ) );
          }
          return id;
        };

        // For each case, which cell edges to connect. Edge ids within a cell:
        // 0 = bottom, 1 = right, 2 = top, 3 = left. -1 terminates the list.
        static const int TBL[ 16 ][ 4 ] = {
          { -1, -1, -1, -1 }, // 0
          {  3,  0, -1, -1 }, // 1  c00
          {  0,  1, -1, -1 }, // 2  c10
          {  3,  1, -1, -1 }, // 3
          {  1,  2, -1, -1 }, // 4  c11
          {  3,  0,  1,  2 }, // 5  saddle
          {  0,  2, -1, -1 }, // 6
          {  2,  3, -1, -1 }, // 7
          {  2,  3, -1, -1 }, // 8  c01
          {  0,  2, -1, -1 }, // 9
          {  0,  1,  2,  3 }, // 10 saddle
          {  1,  2, -1, -1 }, // 11
          {  3,  1, -1, -1 }, // 12
          {  0,  1, -1, -1 }, // 13
          {  3,  0, -1, -1 }, // 14
          { -1, -1, -1, -1 }  // 15
        };

        std::unordered_map< long long, std::vector< long long > > adj;
        auto connect = [&]( long long a, long long b ) {
          adj[ a ].push_back( b );
          adj[ b ].push_back( a );
        };

        for( int i = 0; i < nx; ++i )
          for( int j = 0; j < ny; ++j )
          {
            const double v00 = g.val( i, j ),     v10 = g.val( i + 1, j );
            const double v11 = g.val( i + 1, j + 1 ), v01 = g.val( i, j + 1 );
            const int cs = ( v00 >= 0 ? 1 : 0 ) | ( v10 >= 0 ? 2 : 0 )
                         | ( v11 >= 0 ? 4 : 0 ) | ( v01 >= 0 ? 8 : 0 );
            auto edge = [&]( int e ) -> long long {
              switch( e )
              {
                case 0: return getH( i, j );       // bottom
                case 1: return getV( i + 1, j );   // right
                case 2: return getH( i, j + 1 );   // top
                default: return getV( i, j );      // left
              }
            };
            for( int k = 0; k < 4; k += 2 )
            {
              if( TBL[ cs ][ k ] < 0 ) break;
              connect( edge( TBL[ cs ][ k ] ), edge( TBL[ cs ][ k + 1 ] ) );
            }
          }

        // Trace each connected component into a closed polyline.
        std::vector< TContour > loops;
        std::unordered_set< long long > used;
        for( const auto& kv : adj )
        {
          const long long start = kv.first;
          if( used.count( start ) ) continue;
          TContour loop;
          long long prev = -1, cur = start;
          while( cur != -1 && !used.count( cur ) )
          {
            used.insert( cur );
            loop.push_back( cross[ cur ] );
            long long nxt = -1;
            for( long long cand : adj[ cur ] )
              if( cand != prev ) { nxt = cand; break; }
            prev = cur;
            cur  = nxt;
          }
          if( loop.size( ) >= 3 ) loops.push_back( std::move( loop ) );
        }

        // Largest |area| first.
        std::sort( loops.begin( ), loops.end( ),
                   []( const TContour& a, const TContour& b ) {
                     return std::fabs( loop_signed_area( a ) )
                          > std::fabs( loop_signed_area( b ) );
                   } );
        return loops;
      }
    } // namespace df_detail

    // =========================================================================
    // helpers
    // =========================================================================
    inline double
    DistanceField::signed_area( const TContour& c )
    { return df_detail::loop_signed_area( c ); }

    inline DistanceField::TPoint
    DistanceField::centroid( const TContour& c )
    { return df_detail::centroid_loops( { c } ); }

    // =========================================================================
    // multi-loop SDF interpolation
    // =========================================================================
    inline std::vector< DistanceField::TContour >
    DistanceField::interpolate_loops( const std::vector< TContour >& A_in,
                                      const std::vector< TContour >& B_in,
                                      double t, bool align, double spacing )
    {
      using namespace df_detail;
      if( A_in.empty( ) || B_in.empty( ) ) return {};
      const double h = ( spacing > 0 ) ? spacing : 1.0;

      // Rigid frames: blend the gross translation/rotation separately so the
      // SDF blend only sees the residual shape (avoids the fade in/out
      // artifact). With align == false this is the identity.
      const TPoint cA = centroid_loops( A_in );
      const TPoint cB = centroid_loops( B_in );
      const double aA = align ? principal_angle( A_in, cA ) : 0.0;
      const double aB = align ? principal_angle( B_in, cB ) : 0.0;

      auto to_frame = [&]( const std::vector< TContour >& loops,
                           const TPoint& c, double ang )
      {
        std::vector< TContour > out;
        for( const auto& L : loops )
        {
          TContour o; o.reserve( L.size( ) );
          for( const auto& p : L )
            o.push_back( align
              ? rotate( TPoint( p.x( ) - c.x( ), p.y( ) - c.y( ) ), -ang )
              : p );
          out.push_back( std::move( o ) );
        }
        return out;
      };

      const std::vector< TContour > A = to_frame( A_in, cA, aA );
      const std::vector< TContour > B = to_frame( B_in, cB, aB );

      // Common grid covering both (canonical) shapes plus a margin so the
      // boundary nodes are safely outside (negative) and loops close.
      double minx = 1e300, miny = 1e300, maxx = -1e300, maxy = -1e300;
      for( const auto* S : { &A, &B } )
        for( const auto& L : *S )
          for( const auto& p : L )
          {
            minx = std::min( minx, p.x( ) ); maxx = std::max( maxx, p.x( ) );
            miny = std::min( miny, p.y( ) ); maxy = std::max( maxy, p.y( ) );
          }
      const double pad = 4.0 * h;
      Grid g;
      g.h  = h;
      g.x0 = minx - pad;  g.y0 = miny - pad;
      g.nx = std::max( 2, (int) std::ceil( ( maxx - minx + 2 * pad ) / h ) );
      g.ny = std::max( 2, (int) std::ceil( ( maxy - miny + 2 * pad ) / h ) );
      // Keep the grid bounded if a contour is huge: coarsen instead of exploding.
      const int CAP = 1024;
      if( g.nx > CAP || g.ny > CAP )
      {
        const double sx = ( maxx - minx + 2 * pad ) / CAP;
        const double sy = ( maxy - miny + 2 * pad ) / CAP;
        g.h  = std::max( h, std::max( sx, sy ) );
        g.nx = std::max( 2, (int) std::ceil( ( maxx - minx + 2 * pad ) / g.h ) );
        g.ny = std::max( 2, (int) std::ceil( ( maxy - miny + 2 * pad ) / g.h ) );
      }

      Grid gA = g, gB = g;
      build_field( A, gA );
      build_field( B, gB );

      // Phi_t = (1-t) Phi_A + t Phi_B
      Grid gT = g;
      gT.phi.resize( gA.phi.size( ) );
      for( std::size_t k = 0; k < gA.phi.size( ); ++k )
        gT.phi[ k ] = ( 1.0 - t ) * gA.phi[ k ] + t * gB.phi[ k ];

      std::vector< TContour > loops = extract_loops( gT );

      // Map the canonical result back through the interpolated rigid frame.
      if( align )
      {
        const double ang_t = ( 1.0 - t ) * aA + t * aB;
        const TPoint c_t( ( 1.0 - t ) * cA.x( ) + t * cB.x( ),
                          ( 1.0 - t ) * cA.y( ) + t * cB.y( ) );
        for( auto& L : loops )
          for( auto& p : L )
          {
            const TPoint r = rotate( p, ang_t );
            p = TPoint( r.x( ) + c_t.x( ), r.y( ) + c_t.y( ) );
          }
      }

      // Normalize every loop to CCW so orientation is consistent with the rest
      // of the pipeline (correspondence / outward normals in later stages).
      for( auto& L : loops )
        if( loop_signed_area( L ) < 0.0 )
          std::reverse( L.begin( ), L.end( ) );
      return loops;
    }

    // =========================================================================
    // pairwise SDF interpolation (returns the largest component)
    // =========================================================================
    inline DistanceField::TContour
    DistanceField::interpolate( const TContour& A, const TContour& B,
                                double t, bool align, double spacing )
    {
      auto loops = interpolate_loops( { A }, { B }, t, align, spacing );
      if( loops.empty( ) ) return {};
      return loops.front( );   // already sorted largest-area first
    }
  } // namespace Final
} // namespace pujCGAL

#endif // __DistanceField__hxx__

// eof - DistanceField.hxx
