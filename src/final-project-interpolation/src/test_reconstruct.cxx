// =============================================================================
// Phase 4 test: Poisson reconstruction from a stacked contour cloud.
//
// Build an ellipsoid of revolution (a=b=18, c=28) as a stack of axial circles
// of radius r(z) = a*sqrt(1-((z-c)/c)^2), reconstruct it, and require that the
// result is a CLOSED, volume-bounding triangle mesh whose volume is within a
// generous band of the analytic ellipsoid volume (4/3*pi*a^2*c).
// A small radial jitter avoids a perfectly degenerate (cospherical) cloud.
// Returns 0 on PASS, 1 otherwise.
// =============================================================================
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "Reconstructor.h"

using TRecon   = pujCGAL::Final::Reconstructor;
using TContour = TRecon::TContour;

int main( )
{
  const double a = 18.0, c = 28.0;
  std::mt19937 rng( 11 );
  std::uniform_real_distribution< double > jit( -0.01, 0.01 );

  std::vector< TRecon::Slice > stack;
  for( int z = 2; z <= int( 2 * c ) - 2; ++z )             // skip the very tips
  {
    const double u = ( z - c ) / c;
    const double r = a * std::sqrt( std::max( 0.0, 1.0 - u * u ) );
    if( r < 3.0 ) continue;
    const int m = std::max( 16, (int) std::round( 2.0 * M_PI * r / 2.0 ) );
    TContour ring;
    for( int i = 0; i < m; ++i )
    {
      const double ang = 2.0 * M_PI * i / m;
      const double rr = r * ( 1.0 + jit( rng ) );
      ring.emplace_back( rr * std::cos( ang ), rr * std::sin( ang ) );
    }
    stack.push_back( TRecon::Slice{ ring, (double) z } );
  }

  std::cout << "stack slices = " << stack.size( ) << "\n";

  TRecon::Mesh mesh;
  if( !TRecon::poisson_reconstruct( stack, mesh ) )
  { std::cout << "FAIL: poisson_reconstruct returned false\n"; return 1; }

  bool closed = false, bounds = false;
  TRecon::verify_closed_orientable( mesh, closed, bounds );
  const double vol      = bounds ? TRecon::volume( mesh ) : 0.0;
  const double analytic = 4.0 / 3.0 * M_PI * a * a * c;

  std::cout << "mesh V=" << mesh.number_of_vertices( )
            << " F=" << mesh.number_of_faces( ) << "\n";
  std::cout << "is_closed=" << ( closed ? "yes" : "no" )
            << "  bounds_volume=" << ( bounds ? "yes" : "no" ) << "\n";
  std::cout << "volume=" << vol << " mm^3   analytic ~ " << analytic
            << "   ratio=" << ( analytic > 0 ? vol / analytic : 0.0 ) << "\n";

  bool ok = true;
  if( !closed )            { std::cout << "FAIL: mesh not closed\n";        ok = false; }
  if( !bounds )            { std::cout << "FAIL: does not bound a volume\n"; ok = false; }
  if( mesh.number_of_faces( ) < 50 ) { std::cout << "FAIL: too few faces\n"; ok = false; }
  if( bounds )
  {
    const double ratio = vol / analytic;
    if( ratio < 0.5 || ratio > 1.5 ) { std::cout << "FAIL: volume off\n"; ok = false; }
  }

  std::cout << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok ? 0 : 1;
}

// eof - test_reconstruct.cxx
