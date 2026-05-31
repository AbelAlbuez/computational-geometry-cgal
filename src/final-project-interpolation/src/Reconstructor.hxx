#ifndef __Reconstructor__hxx__
#define __Reconstructor__hxx__

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <tuple>
#include <vector>

#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

namespace pujCGAL
{
  namespace Final
  {
    namespace rc_detail
    {
      // Trailing integer of a stem ("slice_0070" -> 70, else -1).
      inline long stem_index( const std::string& stem )
      {
        std::size_t i = stem.size( );
        while( i > 0 && std::isdigit( (unsigned char) stem[ i - 1 ] ) ) --i;
        if( i == stem.size( ) ) return -1;
        return std::stol( stem.substr( i ) );
      }
    } // namespace rc_detail

    // =========================================================================
    // load_stack
    // =========================================================================
    inline std::vector< Reconstructor::Slice >
    Reconstructor::load_stack( const std::string& dir, double dz )
    {
      namespace fs = std::filesystem;
      std::vector< std::pair< long, std::string > > files;
      for( const auto& e : fs::directory_iterator( dir ) )
      {
        if( e.path( ).extension( ) != ".obj" ) continue;
        files.emplace_back( rc_detail::stem_index( e.path( ).stem( ).string( ) ),
                            e.path( ).string( ) );
      }
      std::sort( files.begin( ), files.end( ) );

      std::vector< Slice > stack;
      for( std::size_t k = 0; k < files.size( ); ++k )
      {
        auto c = ContourInterpolator::read_obj( files[ k ].second );
        if( c.size( ) < 3 ) continue;
        const long idx = files[ k ].first;
        const double z = ( idx >= 0 ? (double) idx : (double) k ) * dz;
        stack.push_back( Slice{ std::move( c ), z } );
      }
      return stack;
    }

    // =========================================================================
    // densify
    // =========================================================================
    inline std::vector< Reconstructor::Slice >
    Reconstructor::densify( const std::vector< Slice >& orig, InterpKind kind, int K )
    {
      std::vector< Slice > out;
      const std::size_t n = orig.size( );
      if( n == 0 ) return out;
      if( n == 1 || K < 1 ) return orig;

      if( kind == InterpKind::Spline || kind == InterpKind::Polynomial )
      {
        std::vector< TContour > slices;
        std::vector< double >   zs;
        for( const auto& s : orig ) { slices.push_back( s.contour ); zs.push_back( s.z ); }
        std::vector< double > q;
        for( std::size_t k = 0; k + 1 < n; ++k )
          for( int j = 0; j < K; ++j )
            q.push_back( zs[ k ] + ( zs[ k + 1 ] - zs[ k ] ) * ( double( j ) / K ) );
        q.push_back( zs.back( ) );
        const auto cs = ContourInterpolator::interpolate_series( slices, zs, q, kind );
        for( std::size_t i = 0; i < cs.size( ); ++i )
          out.push_back( Slice{ cs[ i ], q[ i ] } );
        return out;
      }

      // Pairwise (Linear | Sdf) between consecutive originals.
      for( std::size_t k = 0; k + 1 < n; ++k )
      {
        out.push_back( orig[ k ] );
        for( int j = 1; j < K; ++j )
        {
          const double tau = double( j ) / K;
          const double z = orig[ k ].z + ( orig[ k + 1 ].z - orig[ k ].z ) * tau;
          TContour c = ContourInterpolator::interpolate(
              orig[ k ].contour, orig[ k + 1 ].contour, tau, kind, /*align=*/true );
          out.push_back( Slice{ std::move( c ), z } );
        }
      }
      out.push_back( orig.back( ) );
      return out;
    }

    // =========================================================================
    // poisson_reconstruct
    // =========================================================================
    inline bool
    Reconstructor::poisson_reconstruct( const std::vector< Slice >& stack, Mesh& out )
    {
      // Tuple = (point, normal, seed-in-plane-outward-normal). The seed travels
      // with the point through mst_orient_normals (which may reorder), so we
      // can repair the sign afterwards regardless of reordering.
      using Tuple = std::tuple< Point3, Vector3, Vector3 >;
      using PMap  = CGAL::Nth_of_tuple_property_map< 0, Tuple >;
      using NMap  = CGAL::Nth_of_tuple_property_map< 1, Tuple >;

      std::vector< Tuple > pts;
      for( const auto& s : stack )
      {
        TContour c = s.contour;
        ContourInterpolator::ensure_ccw( c );        // outward normal needs CCW
        const std::size_t m = c.size( );
        if( m < 3 ) continue;
        for( std::size_t i = 0; i < m; ++i )
        {
          const Point2& pm = c[ ( i + m - 1 ) % m ];
          const Point2& pp = c[ ( i + 1 ) % m ];
          const double tx = pp.x( ) - pm.x( ), ty = pp.y( ) - pm.y( );
          const double len = std::hypot( tx, ty );
          if( len <= 0.0 ) continue;
          // Outward normal of a CCW contour = tangent rotated by -90deg.
          const double nx = ty / len, ny = -tx / len;
          pts.emplace_back( Point3( c[ i ].x( ), c[ i ].y( ), s.z ),
                            Vector3( nx, ny, 0.0 ),
                            Vector3( nx, ny, 0.0 ) );
        }
      }
      if( pts.size( ) < 20 ) return false;

      const int kk = std::min< int >( 18, (int) pts.size( ) - 1 );
      const double spacing = CGAL::compute_average_spacing< CGAL::Sequential_tag >(
          pts, 6, CGAL::parameters::point_map( PMap( ) ) );

      CGAL::jet_estimate_normals< CGAL::Sequential_tag >(
          pts, kk, CGAL::parameters::point_map( PMap( ) ).normal_map( NMap( ) ) );
      auto un = CGAL::mst_orient_normals(
          pts, kk, CGAL::parameters::point_map( PMap( ) ).normal_map( NMap( ) ) );
      pts.erase( un, pts.end( ) );
      if( pts.size( ) < 20 ) return false;

      // Repair: on the lateral surface (substantially horizontal normal) force
      // the sign to agree with the in-plane outward seed; leave near-vertical
      // cap normals to the MST result.
      for( auto& tp : pts )
      {
        Vector3& nrm = std::get< 1 >( tp );
        const Vector3& seed = std::get< 2 >( tp );
        const double horiz = nrm.x( ) * nrm.x( ) + nrm.y( ) * nrm.y( );
        if( horiz > 0.25 )
        {
          const double d = nrm.x( ) * seed.x( ) + nrm.y( ) * seed.y( );
          if( d < 0.0 ) nrm = Vector3( -nrm.x( ), -nrm.y( ), -nrm.z( ) );
        }
      }

      out.clear( );
      const bool ok = CGAL::poisson_surface_reconstruction_delaunay(
          pts.begin( ), pts.end( ), PMap( ), NMap( ), out, spacing );
      return ok && out.number_of_faces( ) > 0;
    }

    // =========================================================================
    // I/O + verification
    // =========================================================================
    inline bool
    Reconstructor::write_mesh( const std::string& fname, const Mesh& m )
    {
      return CGAL::IO::write_polygon_mesh(
          fname, m, CGAL::parameters::stream_precision( 17 ) );
    }

    inline bool
    Reconstructor::verify_closed_orientable( Mesh& m, bool& closed,
                                             bool& bounds_volume )
    {
      closed = CGAL::is_closed( m );
      bounds_volume = false;
      if( closed )
      {
        CGAL::Polygon_mesh_processing::orient_to_bound_a_volume( m );
        bounds_volume = CGAL::Polygon_mesh_processing::does_bound_a_volume( m );
      }
      return closed && bounds_volume;
    }

    inline double
    Reconstructor::volume( const Mesh& m )
    {
      return CGAL::to_double( CGAL::Polygon_mesh_processing::volume( m ) );
    }

    // =========================================================================
    // OPTIONAL Boissonnat-style capped band (comparison baseline)
    // =========================================================================
    inline bool
    Reconstructor::boissonnat_band( const std::vector< Slice >& orig, Mesh& out )
    {
      const std::size_t n = orig.size( );
      if( n < 2 ) return false;
      const std::size_t N = 64;                       // common ring resolution

      // Resample + chain the cyclic correspondence so ring index i is coherent.
      std::vector< TContour > al( n );
      {
        TContour c0 = orig[ 0 ].contour;
        ContourInterpolator::ensure_ccw( c0 );
        al[ 0 ] = ContourInterpolator::resample_uniform( c0, N );
      }
      for( std::size_t k = 1; k < n; ++k )
      {
        TContour ck = orig[ k ].contour;
        ContourInterpolator::ensure_ccw( ck );
        TContour rk = ContourInterpolator::resample_uniform( ck, N );
        std::size_t sigma = 0; bool rev = false;
        ContourInterpolator::cyclic_align( al[ k - 1 ], rk, sigma, rev );
        TContour aligned( N );
        for( std::size_t i = 0; i < N; ++i )
        {
          const std::size_t kk = ( i + sigma ) % N;
          aligned[ i ] = rk[ rev ? ( N - 1 - kk ) : kk ];
        }
        al[ k ] = aligned;
      }

      out.clear( );
      std::vector< std::vector< Mesh::Vertex_index > > V( n,
          std::vector< Mesh::Vertex_index >( N ) );
      for( std::size_t k = 0; k < n; ++k )
        for( std::size_t i = 0; i < N; ++i )
          V[ k ][ i ] = out.add_vertex(
              Point3( al[ k ][ i ].x( ), al[ k ][ i ].y( ), orig[ k ].z ) );

      // Lateral quads (two CCW triangles each) between consecutive rings.
      for( std::size_t k = 0; k + 1 < n; ++k )
        for( std::size_t i = 0; i < N; ++i )
        {
          const auto a = V[ k ][ i ],         b = V[ k ][ ( i + 1 ) % N ];
          const auto cc = V[ k + 1 ][ ( i + 1 ) % N ], d = V[ k + 1 ][ i ];
          out.add_face( a, b, cc );
          out.add_face( a, cc, d );
        }

      // Centroid caps (wound so the shared ring edges oppose the lateral faces).
      auto cap = [&]( std::size_t k, bool bottom )
      {
        double sx = 0, sy = 0;
        for( std::size_t i = 0; i < N; ++i ) { sx += al[ k ][ i ].x( ); sy += al[ k ][ i ].y( ); }
        const auto cv = out.add_vertex( Point3( sx / N, sy / N, orig[ k ].z ) );
        for( std::size_t i = 0; i < N; ++i )
        {
          const auto p0 = V[ k ][ i ], p1 = V[ k ][ ( i + 1 ) % N ];
          if( bottom ) out.add_face( cv, p1, p0 );
          else         out.add_face( cv, p0, p1 );
        }
      };
      cap( 0, true );
      cap( n - 1, false );

      return out.number_of_faces( ) > 0;
    }
  } // namespace Final
} // namespace pujCGAL

#endif // __Reconstructor__hxx__

// eof - Reconstructor.hxx
