// =============================================================================
// Phase 2 test: series interpolation along z.
//
// Slices are circles whose center and radius follow smooth, wiggly analytic
// functions of z, sampled on a long window. We hide the true midpoints and
// reconstruct them from neighbours with the natural cubic SPLINE and with the
// (degree-clamped) global POLYNOMIAL. The spline, being a local C2 fit, must
// recover the midpoints with smaller error than the global low-degree
// polynomial -- the textbook spline-vs-Runge/underfit contrast.
// Returns 0 on PASS (spline beats polynomial), 1 otherwise.
// =============================================================================
#include <cmath>
#include <iostream>
#include <vector>

#include "ContourInterpolator.h"

using TInterp   = pujCGAL::Final::ContourInterpolator;
using TContour  = TInterp::TContour;
using InterpKind = TInterp::InterpKind;

namespace
{
  double cx_of( double z ) { return 10.0 * std::sin( 0.6 * z ); }
  double cy_of( double )   { return 0.0; }
  double r_of ( double z ) { return 20.0 + 5.0 * std::sin( 0.9 * z + 1.0 ); }

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

  double mean_radius( const TContour& c, const TInterp::TPoint& ctr )
  {
    double s = 0.0;
    for( const auto& p : c ) s += std::hypot( p.x( ) - ctr.x( ), p.y( ) - ctr.y( ) );
    return s / c.size( );
  }
}

int main( )
{
  const int M = 9;                                   // long window: z = 0..8
  std::vector< TContour > slices;
  std::vector< double >   zs;
  for( int z = 0; z < M; ++z )
  {
    slices.push_back( circle( cx_of( z ), cy_of( z ), r_of( z ), 60 ) );
    zs.push_back( (double) z );
  }

  std::vector< double > q;                           // hidden midpoints
  for( int z = 0; z + 1 < M; ++z ) q.push_back( z + 0.5 );

  const auto spl  = TInterp::interpolate_series( slices, zs, q, InterpKind::Spline );
  const auto poly = TInterp::interpolate_series( slices, zs, q, InterpKind::Polynomial );

  double err_spl = 0.0, err_poly = 0.0;
  for( std::size_t i = 0; i < q.size( ); ++i )
  {
    const double gx = cx_of( q[ i ] ), gy = cy_of( q[ i ] ), gr = r_of( q[ i ] );

    const auto cs = TInterp::centroid( spl[ i ] );
    err_spl += std::hypot( cs.x( ) - gx, cs.y( ) - gy )
             + std::fabs( mean_radius( spl[ i ], cs ) - gr );

    const auto cp = TInterp::centroid( poly[ i ] );
    err_poly += std::hypot( cp.x( ) - gx, cp.y( ) - gy )
              + std::fabs( mean_radius( poly[ i ], cp ) - gr );
  }

  std::cout << "window slices   = " << M << "\n";
  std::cout << "midpoints       = " << q.size( ) << "\n";
  std::cout << "spline   total error = " << err_spl  << "\n";
  std::cout << "polynom. total error = " << err_poly << "  (deg=min(3,m-1))\n";

  bool ok = true;
  if( !( err_spl < err_poly ) )      { std::cout << "FAIL: spline not better than polynomial\n"; ok = false; }
  if( err_spl > 2.0 )                { std::cout << "FAIL: spline error too large\n"; ok = false; }
  std::cout << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok ? 0 : 1;
}

// eof - test_series.cxx
