// =============================================================================
// Phase 3 test: shape-based (signed distance field) interpolation.
//
// Test 1 (geometry): two offset circles of different radius, blended at t=0.5
//   with rigid pre-alignment, must give a circle at the midpoint center with
//   the mean radius.  A: c(0,0) r8   B: c(20,4) r12   ->  c(10,2) r10.
//
// Test 2 (topology robustness): a single circle vs. TWO separated circles.
//   Plain vertex correspondence cannot represent this split, but the SDF zero
//   set naturally has the right number of components. Near the two-circle end
//   the blended field must yield two loops.
// Returns 0 on PASS, 1 otherwise.
// =============================================================================
#include <cmath>
#include <iostream>
#include <vector>

#include "DistanceField.h"

using DF       = pujCGAL::Final::DistanceField;
using TContour = DF::TContour;

namespace
{
  TContour circle( double cx, double cy, double r, int m )
  {
    TContour c; c.reserve( m );
    for( int i = 0; i < m; ++i )
    {
      const double a = 2.0 * M_PI * i / m;
      c.emplace_back( cx + r * std::cos( a ), cy + r * std::sin( a ) );
    }
    return c;
  }

  void radius_stats( const TContour& c, const DF::TPoint& ctr,
                     double& rmean, double& rmin, double& rmax )
  {
    rmean = 0; rmin = 1e18; rmax = 0;
    for( const auto& p : c )
    {
      const double r = std::hypot( p.x( ) - ctr.x( ), p.y( ) - ctr.y( ) );
      rmean += r; rmin = std::min( rmin, r ); rmax = std::max( rmax, r );
    }
    rmean /= c.size( );
  }
}

int main( )
{
  bool ok = true;

  // -- Test 1: offset circles, different radius -----------------------------
  {
    const auto A = circle(  0.0, 0.0,  8.0, 80 );
    const auto B = circle( 20.0, 4.0, 12.0, 90 );
    const auto C = DF::interpolate( A, B, 0.5, /*align=*/true, /*spacing=*/1.0 );

    const auto ctr = DF::centroid( C );
    double rmean, rmin, rmax; radius_stats( C, ctr, rmean, rmin, rmax );
    const double cx_err   = std::hypot( ctr.x( ) - 10.0, ctr.y( ) - 2.0 );
    const double roundness = ( rmax - rmin ) / rmean;

    std::cout << "[T1] vertices=" << C.size( )
              << "  centroid=(" << ctr.x( ) << ", " << ctr.y( ) << ") exp (10,2)"
              << "  err=" << cx_err << "\n"
              << "     radius mean=" << rmean << " exp 10  roundness=" << roundness << "\n";
    if( cx_err > 0.8 )                    { std::cout << "FAIL T1: centroid\n"; ok = false; }
    if( std::fabs( rmean - 10.0 ) > 0.8 ) { std::cout << "FAIL T1: radius\n";   ok = false; }
    if( roundness > 0.12 )                { std::cout << "FAIL T1: roundness\n"; ok = false; }
  }

  // -- Test 2: topology split (1 contour vs 2) ------------------------------
  {
    const std::vector< TContour > A = { circle(   0.0, 0.0, 8.0, 80 ) };
    const std::vector< TContour > B = { circle( -13.0, 0.0, 6.0, 60 ),
                                        circle(  13.0, 0.0, 6.0, 60 ) };
    // Near the two-circle end, with absolute positions kept (no alignment).
    const auto loops = DF::interpolate_loops( A, B, 0.8, /*align=*/false, 1.0 );

    int big = 0;
    for( const auto& L : loops )
      if( std::fabs( DF::signed_area( L ) ) > 20.0 ) ++big;

    std::cout << "[T2] loops=" << loops.size( ) << "  significant components=" << big
              << "  (expect 2)\n";
    if( big != 2 ) { std::cout << "FAIL T2: topology not split into two\n"; ok = false; }
  }

  std::cout << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok ? 0 : 1;
}

// eof - test_sdf.cxx
