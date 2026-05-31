// =============================================================================
// Phase 5 test: metrics on a reconstructed ellipsoid.
//
// Reconstruct an a=b=18, c=28 ellipsoid of revolution (symmetric about all
// three principal axes), then check that the metrics are sane:
//   - volume > 0 and area > 0;
//   - stack-volume cross-check within ~25% of the mesh volume;
//   - best reflective symmetry Dice is high (the ellipsoid is symmetric);
//   - mean curvature stats are finite and positive (a convex blob).
// Returns 0 on PASS, 1 otherwise.
// =============================================================================
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "Reconstructor.h"
#include "Metrics.h"

using TRecon   = pujCGAL::Final::Reconstructor;
using TMetrics = pujCGAL::Final::Metrics;
using TContour = TRecon::TContour;

int main( )
{
  const double a = 18.0, c = 28.0;
  std::mt19937 rng( 5 );
  std::uniform_real_distribution< double > jit( -0.01, 0.01 );

  std::vector< TRecon::Slice > stack;
  for( int z = 2; z <= int( 2 * c ) - 2; ++z )
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

  TRecon::Mesh mesh;
  if( !TRecon::poisson_reconstruct( stack, mesh ) )
  { std::cout << "FAIL: reconstruction\n"; return 1; }

  const auto va  = TMetrics::volume_area( mesh );
  const double sv = TMetrics::stack_volume( stack );
  const auto sym = TMetrics::symmetry_score( mesh, 36 );
  const auto cur = TMetrics::mean_curvature_stats( mesh );
  const double stack_diff = std::fabs( va.volume_mm3 - sv ) / va.volume_mm3;

  std::cout << "volume=" << va.volume_mm3 << " mm^3  area=" << va.area_mm2 << " mm^2\n";
  std::cout << "stack volume=" << sv << " mm^3  diff=" << 100*stack_diff << " %\n";
  std::cout << "symmetry: axis " << sym.axis << "  Dice=" << sym.dice
            << "  Hausdorff=" << sym.hausdorff << " mm\n";
  std::cout << "mean curvature: mean=" << cur.mean << " median=" << cur.median
            << " max=" << cur.max << " n=" << cur.n << "\n";

  bool ok = true;
  if( va.volume_mm3 <= 0 || va.area_mm2 <= 0 ) { std::cout << "FAIL: vol/area\n"; ok = false; }
  if( stack_diff > 0.25 )                      { std::cout << "FAIL: stack cross-check\n"; ok = false; }
  if( sym.dice < 0.75 )                        { std::cout << "FAIL: symmetry Dice low\n"; ok = false; }
  if( !( cur.mean > 0 ) || !std::isfinite( cur.mean ) || cur.n == 0 )
                                               { std::cout << "FAIL: curvature\n"; ok = false; }

  std::cout << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok ? 0 : 1;
}

// eof - test_metrics.cxx
