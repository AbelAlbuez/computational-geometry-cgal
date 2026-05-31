// =============================================================================
// Final project: Geometric Interpolation of Tumor Contours
// Author: Abel Albuez Sanchez
//
// Reads two .obj contours (slice A and slice B), resamples them to a common
// number of vertices, linearly interpolates them at t=0.5, resolves any
// self-intersections of the resulting contour and writes it back to
// output/contour_interpolated.obj.
//
// Usage: ./contour_interpolator slice_A.obj slice_B.obj
// =============================================================================

#include <fstream>
#include <iostream>
#include <string>

#include "ContourInterpolator.h"
#include "ContourResampler.h"
#include "LinearInterpolator.h"
#include "SelfIntersectionResolver.h"

// Writes a 2D contour as a 3D .obj (v x y 0.0) so that 3D viewers such as
// ParaView, Online 3D Viewer or ImageToSTL accept it. ContourInterpolator
// cannot be modified, so we provide this helper locally.
static void write_obj_3d(
  const std::string& fname,
  const pujCGAL::Final::ContourInterpolator::TContour& contour
  )
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
              << " slice_A.obj slice_B.obj" << std::endl;
    return 1;
  }

  using TInterp    = pujCGAL::Final::ContourInterpolator;
  using TResampler = pujCGAL::Final::ContourResampler;
  using TLinear    = pujCGAL::Final::LinearInterpolator;
  using TResolver  = pujCGAL::Final::SelfIntersectionResolver;

  const std::string fa = argv[ 1 ];
  const std::string fb = argv[ 2 ];

  auto A = TInterp::read_obj( fa );
  auto B = TInterp::read_obj( fb );

  std::cout << "Contour A: " << A.size( ) << " vertices\n";
  std::cout << "Contour B: " << B.size( ) << " vertices\n";

  if( A.empty( ) || B.empty( ) )
  {
    std::cerr << "Error: one of the contours is empty." << std::endl;
    return 2;
  }

  // -- Align both contours to the same vertex count by arc-length resampling.
  const int n = static_cast< int >( std::max( A.size( ), B.size( ) ) );
  auto Ar = TResampler::resample( A, n );
  auto Br = TResampler::resample( B, n );
  std::cout << "Resampled to: " << n << " vertices\n";

  // -- Linear interpolation at t = 0.5.
  auto C = TLinear::interpolate( Ar, Br, 0.5 );
  std::cout << "Interpolated contour: " << C.size( ) << " vertices\n";

  // -- Self-intersection resolution.
  const bool had_inter = TResolver::has_self_intersections( C );
  std::cout << "Self-intersections detected: "
            << ( had_inter ? "yes" : "no" ) << "\n";
  if( had_inter )
    C = TResolver::resolve( C );

  const std::string out = "output/contour_interpolated.obj";
  write_obj_3d( out, C );
  std::cout << "Result written to " << out << std::endl;
  return 0;
}

// eof - main.cxx
