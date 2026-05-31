// =============================================================================
// Phase 1 test: linear interpolation of a circle and a translated/scaled
// circle at t=0.5. The LERP of two aligned circles is a circle whose center
// is the midpoint of the centers and whose radius is the mean of the radii.
//   A: center (0,0)  r=10  (64 pts)
//   B: center (5,5)  r=20  (80 pts)
//   expected C at t=0.5: center (2.5,2.5), radius 15.
// Returns 0 on PASS, 1 on FAIL.
// =============================================================================
#include <cmath>
#include <iostream>

#include "ContourInterpolator.h"

using TInterp  = pujCGAL::Final::ContourInterpolator;
using TContour = TInterp::TContour;

namespace
{
  TContour circle( double cx, double cy, double r, int m )
  {
    TContour c;
    c.reserve( m );
    for( int i = 0; i < m; ++i )
    {
      const double a = 2.0 * M_PI * i / m;
      c.emplace_back( cx + r * std::cos( a ), cy + r * std::sin( a ) );
    }
    return c;
  }
}

int main( )
{
  const auto A = circle( 0.0, 0.0, 10.0, 64 );
  const auto B = circle( 5.0, 5.0, 20.0, 80 );
  const auto C = TInterp::interpolate( A, B, 0.5 );

  const auto   ctr  = TInterp::centroid( C );
  const double area = std::fabs( TInterp::signed_area( C ) );

  double rmean = 0.0, rmin = 1e18, rmax = 0.0;
  for( const auto& p : C )
  {
    const double r = std::hypot( p.x( ) - ctr.x( ), p.y( ) - ctr.y( ) );
    rmean += r;
    rmin = std::min( rmin, r );
    rmax = std::max( rmax, r );
  }
  rmean /= static_cast< double >( C.size( ) );

  const double expected_area = M_PI * 15.0 * 15.0;          // ~706.86
  const double cx_err   = std::hypot( ctr.x( ) - 2.5, ctr.y( ) - 2.5 );
  const double area_err = std::fabs( area - expected_area ) / expected_area;
  const double roundness = ( rmax - rmin ) / rmean;          // 0 = perfect circle

  std::cout << "result vertices = " << C.size( ) << "\n";
  std::cout << "centroid  = (" << ctr.x( ) << ", " << ctr.y( )
            << ")   expected (2.5, 2.5)   err = " << cx_err << "\n";
  std::cout << "area      = " << area << "   expected ~ " << expected_area
            << "   relerr = " << area_err << "\n";
  std::cout << "radius    mean = " << rmean << "  min = " << rmin
            << "  max = " << rmax << "  roundness = " << roundness << "\n";

  bool ok = true;
  if( cx_err   > 0.30 ) { std::cout << "FAIL: centroid off\n";  ok = false; }
  if( area_err > 0.05 ) { std::cout << "FAIL: area off\n";      ok = false; }
  if( std::fabs( rmean - 15.0 ) > 0.50 ) { std::cout << "FAIL: radius off\n"; ok = false; }
  if( roundness > 0.06 ) { std::cout << "FAIL: not circle-like\n"; ok = false; }

  std::cout << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok ? 0 : 1;
}

// eof - test_linear.cxx
