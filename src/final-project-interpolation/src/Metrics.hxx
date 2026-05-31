#ifndef __Metrics__hxx__
#define __Metrics__hxx__

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Side_of_triangle_mesh.h>

namespace pujCGAL
{
  namespace Final
  {
    namespace mt_detail
    {
      using K      = ContourInterpolator::TKernel;   // EPIC, same as the mesh
      using Point3 = K::Point_3;
      using Vec3   = K::Vector_3;
    }

    // =========================================================================
    // volume / area
    // =========================================================================
    inline Metrics::VolArea
    Metrics::volume_area( Mesh& m )
    {
      namespace PMP = CGAL::Polygon_mesh_processing;
      PMP::orient_to_bound_a_volume( m );
      VolArea r;
      r.volume_mm3 = std::fabs( CGAL::to_double( PMP::volume( m ) ) );
      r.area_mm2   = CGAL::to_double( PMP::area( m ) );
      return r;
    }

    // =========================================================================
    // stack volume (trapezoidal rule on shoelace areas)
    // =========================================================================
    inline double
    Metrics::stack_volume( const std::vector< Slice >& stack )
    {
      std::vector< Slice > s = stack;
      std::sort( s.begin( ), s.end( ),
                 []( const Slice& a, const Slice& b ) { return a.z < b.z; } );
      double V = 0.0;
      for( std::size_t k = 0; k + 1 < s.size( ); ++k )
      {
        const double a0 = std::fabs( ContourInterpolator::signed_area( s[ k ].contour ) );
        const double a1 = std::fabs( ContourInterpolator::signed_area( s[ k + 1 ].contour ) );
        const double dz = s[ k + 1 ].z - s[ k ].z;
        V += 0.5 * ( a0 + a1 ) * dz;
      }
      return V;
    }

    // =========================================================================
    // reflective symmetry score
    // =========================================================================
    inline Metrics::Symmetry
    Metrics::symmetry_score( const Mesh& m, int grid_res )
    {
      using namespace mt_detail;
      namespace PMP = CGAL::Polygon_mesh_processing;
      Symmetry best;
      if( m.number_of_faces( ) == 0 ) return best;

      // Area-weighted centroid and inertia tensor over the triangle faces.
      Eigen::Vector3d C = Eigen::Vector3d::Zero( );
      double mass = 0.0;
      std::vector< std::pair< Eigen::Vector3d, double > > fg; // (centroid, area)
      for( auto f : m.faces( ) )
      {
        std::vector< Point3 > p;
        for( auto v : CGAL::vertices_around_face( m.halfedge( f ), m ) )
          p.push_back( m.point( v ) );
        if( p.size( ) != 3 ) continue;
        const Eigen::Vector3d a( p[0].x( ), p[0].y( ), p[0].z( ) );
        const Eigen::Vector3d b( p[1].x( ), p[1].y( ), p[1].z( ) );
        const Eigen::Vector3d c( p[2].x( ), p[2].y( ), p[2].z( ) );
        const double area = 0.5 * ( b - a ).cross( c - a ).norm( );
        const Eigen::Vector3d g = ( a + b + c ) / 3.0;
        fg.emplace_back( g, area );
        C += area * g; mass += area;
      }
      if( mass <= 0.0 ) return best;
      C /= mass;

      Eigen::Matrix3d I = Eigen::Matrix3d::Zero( );
      for( const auto& it : fg )
      {
        const Eigen::Vector3d r = it.first - C;
        I += it.second * ( r.squaredNorm( ) * Eigen::Matrix3d::Identity( )
                           - r * r.transpose( ) );
      }
      Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > es( I );
      const Eigen::Matrix3d axes = es.eigenvectors( );

      // Bounding box (padded) shared by the inside tests.
      double mnx = 1e300, mny = 1e300, mnz = 1e300;
      double mxx = -1e300, mxy = -1e300, mxz = -1e300;
      for( auto v : m.vertices( ) )
      {
        const Point3& p = m.point( v );
        mnx = std::min( mnx, p.x( ) ); mxx = std::max( mxx, p.x( ) );
        mny = std::min( mny, p.y( ) ); mxy = std::max( mxy, p.y( ) );
        mnz = std::min( mnz, p.z( ) ); mxz = std::max( mxz, p.z( ) );
      }
      const double ext = std::max( { mxx - mnx, mxy - mny, mxz - mnz } );
      const double pad = 0.05 * ext;
      mnx -= pad; mny -= pad; mnz -= pad;
      mxx += pad; mxy += pad; mxz += pad;
      const double h = ext / std::max( 4, grid_res );

      CGAL::Side_of_triangle_mesh< Mesh, K > inside_m( m );

      for( int a = 0; a < 3; ++a )
      {
        const Eigen::Vector3d u = axes.col( a ).normalized( );

        // Reflect a copy of the mesh about the plane (C, u).
        Mesh mp = m;
        for( auto v : mp.vertices( ) )
        {
          const Point3& p = mp.point( v );
          const Eigen::Vector3d r( p.x( ), p.y( ), p.z( ) );
          const double d = ( r - C ).dot( u );
          const Eigen::Vector3d q = r - 2.0 * d * u;
          mp.point( v ) = Point3( q.x( ), q.y( ), q.z( ) );
        }
        CGAL::Side_of_triangle_mesh< Mesh, K > inside_mp( mp );

        std::size_t in_m = 0, in_mp = 0, in_both = 0;
        for( double x = mnx + 0.5 * h; x < mxx; x += h )
          for( double y = mny + 0.5 * h; y < mxy; y += h )
            for( double z = mnz + 0.5 * h; z < mxz; z += h )
            {
              const Point3 q( x, y, z );
              const bool a_in = inside_m( q )  == CGAL::ON_BOUNDED_SIDE;
              const bool b_in = inside_mp( q ) == CGAL::ON_BOUNDED_SIDE;
              if( a_in )          ++in_m;
              if( b_in )          ++in_mp;
              if( a_in && b_in )  ++in_both;
            }
        const double dice = ( in_m + in_mp > 0 )
            ? ( 2.0 * in_both / double( in_m + in_mp ) ) : 0.0;
        const double hd = CGAL::Polygon_mesh_processing::
            approximate_symmetric_Hausdorff_distance< CGAL::Sequential_tag >( m, mp );

        if( best.axis < 0 || dice > best.dice )
        { best.axis = a; best.dice = dice; best.hausdorff = hd; }
      }
      return best;
    }

    // =========================================================================
    // discrete mean curvature (cotangent Laplacian) statistics
    // =========================================================================
    inline Metrics::CurvStats
    Metrics::mean_curvature_stats( const Mesh& m )
    {
      using namespace mt_detail;
      CurvStats st;
      const std::size_t nv = m.number_of_vertices( );
      if( nv == 0 ) return st;

      std::vector< Eigen::Vector3d > L( nv, Eigen::Vector3d::Zero( ) );
      std::vector< double >          A( nv, 0.0 );

      auto idx = []( Mesh::Vertex_index v ) { return (std::size_t) v.idx( ); };
      auto ev  = []( const Point3& p ) { return Eigen::Vector3d( p.x( ), p.y( ), p.z( ) ); };

      for( auto f : m.faces( ) )
      {
        std::vector< Mesh::Vertex_index > vs;
        for( auto v : CGAL::vertices_around_face( m.halfedge( f ), m ) ) vs.push_back( v );
        if( vs.size( ) != 3 ) continue;
        const Eigen::Vector3d p0 = ev( m.point( vs[0] ) );
        const Eigen::Vector3d p1 = ev( m.point( vs[1] ) );
        const Eigen::Vector3d p2 = ev( m.point( vs[2] ) );
        const Eigen::Vector3d cr = ( p1 - p0 ).cross( p2 - p0 );
        const double twoArea = cr.norm( );
        if( twoArea <= 0.0 ) continue;
        const double area = 0.5 * twoArea;

        // cot at each vertex = dot(edge1,edge2) / |cross(edge1,edge2)|
        auto cot = []( const Eigen::Vector3d& u, const Eigen::Vector3d& w ) {
          const double c = u.cross( w ).norm( );
          return ( c > 0.0 ) ? ( u.dot( w ) / c ) : 0.0;
        };
        const double cot0 = cot( p1 - p0, p2 - p0 ); // angle at vs[0], edge (1,2)
        const double cot1 = cot( p0 - p1, p2 - p1 ); // angle at vs[1], edge (0,2)
        const double cot2 = cot( p0 - p2, p1 - p2 ); // angle at vs[2], edge (0,1)

        const std::size_t i0 = idx( vs[0] ), i1 = idx( vs[1] ), i2 = idx( vs[2] );
        L[ i1 ] += cot0 * ( p2 - p1 );  L[ i2 ] += cot0 * ( p1 - p2 );
        L[ i0 ] += cot1 * ( p2 - p0 );  L[ i2 ] += cot1 * ( p0 - p2 );
        L[ i0 ] += cot2 * ( p1 - p0 );  L[ i1 ] += cot2 * ( p0 - p1 );
        A[ i0 ] += area / 3.0;  A[ i1 ] += area / 3.0;  A[ i2 ] += area / 3.0;
      }

      std::vector< double > H;
      H.reserve( nv );
      for( std::size_t i = 0; i < nv; ++i )
        if( A[ i ] > 0.0 )
          H.push_back( L[ i ].norm( ) / ( 4.0 * A[ i ] ) );  // |Delta p| / 2, A mixed
      if( H.empty( ) ) return st;

      std::sort( H.begin( ), H.end( ) );
      double sum = 0.0; for( double h : H ) sum += h;
      const double mean = sum / H.size( );
      double var = 0.0; for( double h : H ) var += ( h - mean ) * ( h - mean );
      st.n      = H.size( );
      st.mean   = mean;
      st.stddev = std::sqrt( var / H.size( ) );
      st.min    = H.front( );
      st.max    = H.back( );
      st.median = H[ H.size( ) / 2 ];
      st.p90    = H[ std::min( H.size( ) - 1, (std::size_t)( 0.9 * H.size( ) ) ) ];
      return st;
    }
  } // namespace Final
} // namespace pujCGAL

#endif // __Metrics__hxx__

// eof - Metrics.hxx
