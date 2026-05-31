// =============================================================================
// Final project: Geometric Interpolation of Tumor Contours
// Author: Abel Albuez Sanchez
//
// Reads two .obj contours (slice A and slice B), aligns them (orientation,
// centroid, optimal start vertex), resamples them to a common number of
// vertices, linearly interpolates them at t=0.5, resolves any
// self-intersections of the resulting contour and writes it back to
// output/contour_interpolated.obj (or to the path given as 3rd argument).
//
// Usage: ./interpolation_lineal slice_A.obj slice_B.obj [output.obj]
// =============================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

// Sibling headers (interpolation-lineal/).
#include "ContourResampler.h"
#include "LinearInterpolator.h"
#include "SelfIntersectionResolver.h"

// Parent module header.
#include "../ContourInterpolator.h"

using TInterp    = pujCGAL::Final::ContourInterpolator;
using TResampler = pujCGAL::Final::ContourResampler;
using TLinear    = pujCGAL::Final::LinearInterpolator;
using TResolver  = pujCGAL::Final::SelfIntersectionResolver;
using TContour   = TInterp::TContour;
using TPoint     = TInterp::TPoint;

// -----------------------------------------------------------------------------
// Helpers for pre-interpolation alignment.
// -----------------------------------------------------------------------------

// Centroid of a contour (average of vertices).
static TPoint centroid( const TContour& c )
{
  double cx = 0.0, cy = 0.0;
  for( const auto& p : c )
  {
    cx += CGAL::to_double( p.x( ) );
    cy += CGAL::to_double( p.y( ) );
  }
  const double n = static_cast< double >( c.size( ) );
  return TPoint( cx / n, cy / n );
}

// Translate a contour by (dx, dy).
static TContour translate( const TContour& c, double dx, double dy )
{
  TContour out;
  out.reserve( c.size( ) );
  for( const auto& p : c )
    out.emplace_back(
      CGAL::to_double( p.x( ) ) + dx,
      CGAL::to_double( p.y( ) ) + dy
      );
  return out;
}

// Signed area via shoelace; positive means CCW, negative means CW.
static double signed_area( const TContour& c )
{
  const std::size_t n = c.size( );
  double s = 0.0;
  for( std::size_t i = 0; i < n; ++i )
  {
    const auto& a = c[ i ];
    const auto& b = c[ ( i + 1 ) % n ];
    s += CGAL::to_double( a.x( ) ) * CGAL::to_double( b.y( ) )
       - CGAL::to_double( b.x( ) ) * CGAL::to_double( a.y( ) );
  }
  return 0.5 * s;
}

// Force counter-clockwise orientation.
static TContour to_ccw( const TContour& c )
{
  if( signed_area( c ) < 0.0 )
    return TContour( c.rbegin( ), c.rend( ) );
  return c;
}

// Find the cyclic rotation k of B that minimises
// sum_i |A[i] - B[(i+k) % n]|^2.
static int best_rotation( const TContour& A, const TContour& B )
{
  const int n = static_cast< int >( A.size( ) );
  int best_k = 0;
  double best_dist = std::numeric_limits< double >::max( );
  for( int k = 0; k < n; ++k )
  {
    double dist = 0.0;
    for( int i = 0; i < n; ++i )
      dist += CGAL::to_double(
        CGAL::squared_distance( A[ i ], B[ ( i + k ) % n ] )
        );
    if( dist < best_dist )
    {
      best_dist = dist;
      best_k = k;
    }
  }
  return best_k;
}

// Rotate B by k positions: result[i] = B[(i+k) % n].
static TContour rotate( const TContour& B, int k )
{
  const int n = static_cast< int >( B.size( ) );
  TContour out;
  out.reserve( n );
  for( int i = 0; i < n; ++i )
    out.push_back( B[ ( i + k ) % n ] );
  return out;
}

// Writes a 2D contour as a 3D .obj (v x y 0.0) so that 3D viewers such as
// ParaView, Online 3D Viewer or ImageToSTL accept it. ContourInterpolator
// cannot be modified, so we provide this helper locally.
static void write_obj_3d( const std::string& fname, const TContour& contour )
{
  std::ofstream ofs( fname );
  if( !ofs )
  {
    std::cerr << "Error: could not write " << fname << std::endl;
    return;
  }

  ofs << "# Contorno interpolado 3D - " << contour.size( ) << " vertices\n";
  for( const auto& p : contour )
    ofs << "v " << p.x( ) << " " << p.y( ) << " 0.0\n";

  const std::size_t n = contour.size( );
  for( std::size_t i = 1; i < n; ++i )
    ofs << "l " << i << " " << ( i + 1 ) << "\n";
  if( n >= 2 )
    ofs << "l " << n << " 1\n";
}

int main( int argc, char** argv )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ]
              << " slice_A.obj slice_B.obj [output.obj]" << std::endl;
    return 1;
  }

  const std::string fa = argv[ 1 ];
  const std::string fb = argv[ 2 ];
  const std::string output_path =
    ( argc >= 4 ) ? argv[ 3 ] : std::string( "output/contour_interpolated.obj" );

  auto A = TInterp::read_obj( fa );
  auto B = TInterp::read_obj( fb );

  std::cout << "Contour A: " << A.size( ) << " vertices\n";
  std::cout << "Contour B: " << B.size( ) << " vertices\n";

  if( A.empty( ) || B.empty( ) )
  {
    std::cerr << "Error: one of the contours is empty." << std::endl;
    return 2;
  }

  // -- (1) Normalise orientation: force both contours to be CCW.
  A = to_ccw( A );
  B = to_ccw( B );

  // -- (2) Centroid-based alignment.
  const TPoint cA = centroid( A );
  const TPoint cB = centroid( B );
  const double cAx = CGAL::to_double( cA.x( ) );
  const double cAy = CGAL::to_double( cA.y( ) );
  const double cBx = CGAL::to_double( cB.x( ) );
  const double cBy = CGAL::to_double( cB.y( ) );
  auto A0 = translate( A, -cAx, -cAy );
  auto B0 = translate( B, -cBx, -cBy );

  // -- (3) Resample both contours to the same number of vertices.
  const int n = static_cast< int >( std::max( A0.size( ), B0.size( ) ) );
  auto Ar = TResampler::resample( A0, n );
  auto Br = TResampler::resample( B0, n );
  std::cout << "Resampled to: " << n << " vertices\n";

  // -- (4) Optimal cyclic rotation of B.
  const int k = best_rotation( Ar, Br );
  if( k != 0 )
  {
    Br = rotate( Br, k );
    std::cout << "Optimal rotation of B: " << k << "\n";
  }
  else
  {
    std::cout << "Optimal rotation of B: 0\n";
  }

  // -- (5) Linear interpolation at t = 0.5 in the centred frame.
  auto C0 = TLinear::interpolate( Ar, Br, 0.5 );

  // -- (6) Translate the result to the average centroid.
  const double mx = 0.5 * ( cAx + cBx );
  const double my = 0.5 * ( cAy + cBy );
  auto C = translate( C0, mx, my );
  std::cout << "Interpolated contour: " << C.size( ) << " vertices\n";

  // -- (7) Self-intersection resolution.
  const bool had_inter = TResolver::has_self_intersections( C );
  std::cout << "Self-intersections detected: "
            << ( had_inter ? "yes" : "no" ) << "\n";
  if( had_inter )
    C = TResolver::resolve( C );

  write_obj_3d( output_path, C );
  std::cout << "Result written to " << output_path << std::endl;
  return 0;
}

// eof - main.cxx
