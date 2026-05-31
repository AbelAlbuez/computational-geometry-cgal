// =============================================================================
// Final project: Geometric Interpolation of Tumor Contours
// Author: Abel Albuez Sanchez
//
// CLI entry point. Modes:
//   linear A.obj B.obj t out.obj [--adaptive [lambda]]
//   sdf    A.obj B.obj t out.obj [--no-align] [--spacing s]
//   series --kind {spline|poly|linear} --window dir --upsample K out_dir
// =============================================================================

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "ContourInterpolator.h"
#include "Reconstructor.h"
#include "Metrics.h"

using TInterp   = pujCGAL::Final::ContourInterpolator;
using TRecon    = pujCGAL::Final::Reconstructor;
using TMetrics  = pujCGAL::Final::Metrics;
using InterpKind = TInterp::InterpKind;
namespace fs = std::filesystem;

namespace
{
  void usage( const char* prog )
  {
    std::cerr
      << "Usage:\n"
      << "  " << prog << " linear A.obj B.obj t out.obj [--adaptive [lambda]]\n"
      << "  " << prog << " sdf    A.obj B.obj t out.obj [--no-align] [--spacing s]\n"
      << "  " << prog << " series --kind {spline|poly|linear} --window dir "
         "--upsample K out_dir\n"
      << "  " << prog << " reconstruct --contours dir --method {linear|spline|sdf} "
         "--out mesh.off [--dz mm] [--upsample K] [--boissonnat band.off]\n"
      << "  " << prog << " metrics --mesh mesh.off [--contours dir --dz mm] [--grid N]\n";
  }

  InterpKind parse_method( const std::string& v, bool& okk )
  {
    okk = true;
    if( v == "linear" ) return InterpKind::Linear;
    if( v == "spline" ) return InterpKind::Spline;
    if( v == "sdf" )    return InterpKind::Sdf;
    if( v == "poly" || v == "polynomial" ) return InterpKind::Polynomial;
    okk = false; return InterpKind::Linear;
  }

  std::string strip_ext( const std::string& p )
  {
    const std::size_t dot = p.find_last_of( '.' );
    const std::size_t sep = p.find_last_of( "/\\" );
    if( dot == std::string::npos || ( sep != std::string::npos && dot < sep ) ) return p;
    return p.substr( 0, dot );
  }

  // ----- linear (M1) -----------------------------------------------------
  int run_linear( int argc, char** argv )
  {
    if( argc < 6 ) { usage( argv[ 0 ] ); return 1; }
    const std::string fa = argv[ 2 ], fb = argv[ 3 ], out = argv[ 5 ];
    const double t = std::atof( argv[ 4 ] );

    bool adaptive = false; double lambda = 2.0;
    for( int i = 6; i < argc; ++i )
    {
      const std::string a = argv[ i ];
      if( a == "--adaptive" )
      {
        adaptive = true;
        if( i + 1 < argc )
        {
          char* end = nullptr;
          const double v = std::strtod( argv[ i + 1 ], &end );
          if( end != argv[ i + 1 ] ) { lambda = v; ++i; }
        }
      }
      else { std::cerr << "Unknown argument: " << a << "\n"; usage( argv[ 0 ] ); return 1; }
    }

    auto A = TInterp::read_obj( fa );
    auto B = TInterp::read_obj( fb );
    if( A.size( ) < 3 || B.size( ) < 3 )
    { std::cerr << "Error: a contour has < 3 vertices.\n"; return 2; }

    const auto C = TInterp::interpolate( A, B, t, adaptive, lambda );
    TInterp::write_obj( out, C );
    const auto cc = TInterp::centroid( C );
    std::cout << "linear t=" << t;
    if( adaptive ) std::cout << "  [curvature-adaptive lambda=" << lambda << "]";
    std::cout << "\n  A: " << A.size( ) << " v  area=" << TInterp::signed_area( A )
              << "\n  B: " << B.size( ) << " v  area=" << TInterp::signed_area( B )
              << "\n  C: " << C.size( ) << " v  area=" << TInterp::signed_area( C )
              << "  centroid=(" << cc.x( ) << ", " << cc.y( ) << ")\n"
              << "  wrote " << out << "\n";
    return 0;
  }

  // ----- sdf (M3) --------------------------------------------------------
  int run_sdf( int argc, char** argv )
  {
    if( argc < 6 ) { usage( argv[ 0 ] ); return 1; }
    const std::string fa = argv[ 2 ], fb = argv[ 3 ], out = argv[ 5 ];
    const double t = std::atof( argv[ 4 ] );

    bool align = true; double spacing = 1.0;
    for( int i = 6; i < argc; ++i )
    {
      const std::string a = argv[ i ];
      if( a == "--no-align" ) align = false;
      else if( a == "--spacing" && i + 1 < argc ) spacing = std::atof( argv[ ++i ] );
      else { std::cerr << "Unknown argument: " << a << "\n"; usage( argv[ 0 ] ); return 1; }
    }

    auto A = TInterp::read_obj( fa );
    auto B = TInterp::read_obj( fb );
    if( A.size( ) < 3 || B.size( ) < 3 )
    { std::cerr << "Error: a contour has < 3 vertices.\n"; return 2; }

    const auto C = TInterp::interpolate( A, B, t, InterpKind::Sdf, align );
    TInterp::write_obj( out, C );
    const auto cc = TInterp::centroid( C );
    std::cout << "sdf t=" << t << ( align ? "" : "  [no-align]" )
              << "  spacing=" << spacing
              << "\n  A: " << A.size( ) << " v  area=" << TInterp::signed_area( A )
              << "\n  B: " << B.size( ) << " v  area=" << TInterp::signed_area( B )
              << "\n  C: " << C.size( ) << " v  area=" << TInterp::signed_area( C )
              << "  centroid=(" << cc.x( ) << ", " << cc.y( ) << ")\n"
              << "  wrote " << out << "\n";
    return 0;
  }

  // ----- series (M2) -----------------------------------------------------
  // Parse a trailing integer from a stem like "slice_0070" -> 70 (else -1).
  long stem_index( const std::string& stem )
  {
    std::size_t i = stem.size( );
    while( i > 0 && std::isdigit( (unsigned char) stem[ i - 1 ] ) ) --i;
    if( i == stem.size( ) ) return -1;
    return std::stol( stem.substr( i ) );
  }

  int run_series( int argc, char** argv )
  {
    InterpKind kind = InterpKind::Spline;
    std::string window, out_dir;
    int K = 4;
    for( int i = 2; i < argc; ++i )
    {
      const std::string a = argv[ i ];
      if( a == "--kind" && i + 1 < argc )
      {
        const std::string v = argv[ ++i ];
        if( v == "spline" )      kind = InterpKind::Spline;
        else if( v == "poly" || v == "polynomial" ) kind = InterpKind::Polynomial;
        else if( v == "linear" ) kind = InterpKind::Linear;
        else { std::cerr << "Unknown --kind: " << v << "\n"; return 1; }
      }
      else if( a == "--window"   && i + 1 < argc ) window  = argv[ ++i ];
      else if( a == "--upsample" && i + 1 < argc ) K       = std::atoi( argv[ ++i ] );
      else if( a.rfind( "--", 0 ) == 0 ) { std::cerr << "Unknown flag: " << a << "\n"; return 1; }
      else out_dir = a;                              // positional out_dir
    }
    if( window.empty( ) || out_dir.empty( ) || K < 1 )
    { usage( argv[ 0 ] ); return 1; }

    // Collect slice_*.obj, sorted, with their z indices.
    std::vector< std::pair< long, std::string > > files;
    for( const auto& e : fs::directory_iterator( window ) )
    {
      if( e.path( ).extension( ) != ".obj" ) continue;
      const long z = stem_index( e.path( ).stem( ).string( ) );
      files.emplace_back( z, e.path( ).string( ) );
    }
    std::sort( files.begin( ), files.end( ) );
    if( files.size( ) < 2 )
    { std::cerr << "Need >= 2 .obj slices in " << window << "\n"; return 2; }

    std::vector< TInterp::TContour > slices;
    std::vector< double > zs;
    for( std::size_t i = 0; i < files.size( ); ++i )
    {
      auto c = TInterp::read_obj( files[ i ].second );
      if( c.size( ) < 3 ) continue;
      slices.push_back( std::move( c ) );
      zs.push_back( files[ i ].first >= 0 ? (double) files[ i ].first : (double) i );
    }

    // Query heights: originals plus K-1 intermediates in every gap.
    std::vector< double > query;
    for( std::size_t k = 0; k + 1 < zs.size( ); ++k )
      for( int j = 0; j < K; ++j )
        query.push_back( zs[ k ] + ( zs[ k + 1 ] - zs[ k ] ) * ( double( j ) / K ) );
    query.push_back( zs.back( ) );

    const auto result = TInterp::interpolate_series( slices, zs, query, kind );

    fs::create_directories( out_dir );
    for( std::size_t i = 0; i < result.size( ); ++i )
    {
      char name[ 64 ];
      std::snprintf( name, sizeof( name ), "interp_%04zu.obj", i );
      TInterp::write_obj( ( fs::path( out_dir ) / name ).string( ), result[ i ] );
    }
    std::cout << "series kind="
              << ( kind == InterpKind::Spline ? "spline"
                 : kind == InterpKind::Polynomial ? "poly" : "linear" )
              << "  slices=" << slices.size( ) << "  upsample=" << K
              << "  -> " << result.size( ) << " contours in " << out_dir << "\n";
    return 0;
  }
  // ----- reconstruct (Stage 3-4: Poisson) --------------------------------
  int run_reconstruct( int argc, char** argv )
  {
    std::string contours, out;
    InterpKind kind = InterpKind::Linear;
    double dz = 1.0;            // BraTS GLI is 1 mm isotropic by default
    int K = 4;
    std::string boissonnat;
    for( int i = 2; i < argc; ++i )
    {
      const std::string a = argv[ i ];
      if( a == "--contours" && i + 1 < argc ) contours = argv[ ++i ];
      else if( a == "--method" && i + 1 < argc )
      {
        bool okk; kind = parse_method( argv[ ++i ], okk );
        if( !okk ) { std::cerr << "Unknown --method\n"; return 1; }
      }
      else if( a == "--out"        && i + 1 < argc ) out        = argv[ ++i ];
      else if( a == "--dz"         && i + 1 < argc ) dz         = std::atof( argv[ ++i ] );
      else if( a == "--upsample"   && i + 1 < argc ) K          = std::atoi( argv[ ++i ] );
      else if( a == "--boissonnat" && i + 1 < argc ) boissonnat = argv[ ++i ];
      else { std::cerr << "Unknown argument: " << a << "\n"; usage( argv[ 0 ] ); return 1; }
    }
    if( contours.empty( ) || out.empty( ) ) { usage( argv[ 0 ] ); return 1; }

    auto orig = TRecon::load_stack( contours, dz );
    if( orig.size( ) < 2 )
    { std::cerr << "Need >= 2 contour slices in " << contours << "\n"; return 2; }

    const auto stack = TRecon::densify( orig, kind, K );
    std::cout << "reconstruct method="
              << ( kind == InterpKind::Spline ? "spline"
                 : kind == InterpKind::Sdf ? "sdf"
                 : kind == InterpKind::Polynomial ? "poly" : "linear" )
              << "  originals=" << orig.size( )
              << "  densified=" << stack.size( )
              << "  dz=" << dz << " mm  upsample=" << K << "\n";

    TRecon::Mesh mesh;
    if( !TRecon::poisson_reconstruct( stack, mesh ) )
    { std::cerr << "Poisson reconstruction failed.\n"; return 3; }

    bool closed = false, bounds = false;
    TRecon::verify_closed_orientable( mesh, closed, bounds );

    const std::string stem = strip_ext( out );
    TRecon::write_mesh( stem + ".off", mesh );
    TRecon::write_mesh( stem + ".ply", mesh );

    std::cout << "  mesh: V=" << mesh.number_of_vertices( )
              << " F=" << mesh.number_of_faces( )
              << "  is_closed=" << ( closed ? "yes" : "no" )
              << "  bounds_volume=" << ( bounds ? "yes" : "no" );
    if( bounds ) std::cout << "  volume=" << TRecon::volume( mesh ) << " mm^3";
    std::cout << "\n  wrote " << stem << ".off / " << stem << ".ply\n";

    if( !boissonnat.empty( ) )
    {
      TRecon::Mesh band;
      if( TRecon::boissonnat_band( orig, band ) )
      {
        bool bc = false, bb = false;
        TRecon::verify_closed_orientable( band, bc, bb );
        TRecon::write_mesh( boissonnat, band );
        std::cout << "  boissonnat band: V=" << band.number_of_vertices( )
                  << " F=" << band.number_of_faces( )
                  << "  is_closed=" << ( bc ? "yes" : "no" )
                  << "  -> " << boissonnat << "\n";
      }
    }
    return ( closed && bounds ) ? 0 : 4;
  }
  // ----- metrics (Stage 5) -----------------------------------------------
  int run_metrics( int argc, char** argv )
  {
    std::string mesh_file, contours;
    double dz = 1.0;
    int grid = 40;
    for( int i = 2; i < argc; ++i )
    {
      const std::string a = argv[ i ];
      if( a == "--mesh"     && i + 1 < argc ) mesh_file = argv[ ++i ];
      else if( a == "--contours" && i + 1 < argc ) contours = argv[ ++i ];
      else if( a == "--dz"       && i + 1 < argc ) dz       = std::atof( argv[ ++i ] );
      else if( a == "--grid"     && i + 1 < argc ) grid     = std::atoi( argv[ ++i ] );
      else { std::cerr << "Unknown argument: " << a << "\n"; usage( argv[ 0 ] ); return 1; }
    }
    if( mesh_file.empty( ) ) { usage( argv[ 0 ] ); return 1; }

    TRecon::Mesh m;
    if( !CGAL::IO::read_polygon_mesh( mesh_file, m ) || m.number_of_faces( ) == 0 )
    { std::cerr << "Could not read a triangle mesh from " << mesh_file << "\n"; return 2; }

    const auto va  = TMetrics::volume_area( m );
    const auto sym = TMetrics::symmetry_score( m, grid );
    const auto cur = TMetrics::mean_curvature_stats( m );

    std::cout << "=== Metrics: " << mesh_file << " ===\n";
    std::cout << "  volume         : " << va.volume_mm3 << " mm^3\n";
    std::cout << "  surface area   : " << va.area_mm2   << " mm^2\n";
    if( !contours.empty( ) )
    {
      const auto stack = TRecon::load_stack( contours, dz );
      const double sv = TMetrics::stack_volume( stack );
      const double pd = ( va.volume_mm3 > 0 )
          ? 100.0 * ( va.volume_mm3 - sv ) / va.volume_mm3 : 0.0;
      std::cout << "  stack volume   : " << sv << " mm^3   (mesh vs stack diff "
                << pd << " %)\n";
    }
    std::cout << "  symmetry       : axis " << sym.axis << "  Dice=" << sym.dice
              << "  Hausdorff=" << sym.hausdorff << " mm\n";
    std::cout << "  mean curvature : mean=" << cur.mean << "  std=" << cur.stddev
              << "  min=" << cur.min << "  max=" << cur.max
              << "  median=" << cur.median << "  p90=" << cur.p90
              << "  (1/mm, n=" << cur.n << ")\n";
    return 0;
  }
} // anonymous namespace

int main( int argc, char** argv )
{
  if( argc >= 2 )
  {
    const std::string mode = argv[ 1 ];
    if( mode == "linear" )      return run_linear( argc, argv );
    if( mode == "sdf" )         return run_sdf( argc, argv );
    if( mode == "series" )      return run_series( argc, argv );
    if( mode == "reconstruct" ) return run_reconstruct( argc, argv );
    if( mode == "metrics" )     return run_metrics( argc, argv );
  }
  usage( argv[ 0 ] );
  return 1;
}

// eof - main.cxx
