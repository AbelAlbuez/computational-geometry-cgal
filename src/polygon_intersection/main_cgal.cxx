// =============================================================================
// Workshop: Convex Polygon Intersection — CGAL native version
// Author: Abel Albuez Sanchez
//
// Uses CGAL::Polygon_2 and CGAL::Polygon_with_holes_2 + boolean intersection.
//
// Usage: ./polygon_intersection_cgal polygon_P.obj polygon_Q.obj output.obj
// =============================================================================

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Exact arithmetic kernel — required for Polygon_set_2 boolean operations.
// Cartesian<long double> is not exact enough for boolean set operations.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>

using TKernel           = CGAL::Exact_predicates_exact_constructions_kernel;
using TPoint            = TKernel::Point_2;
using TPolygon          = CGAL::Polygon_2< TKernel >;
using TPolygonWithHoles = CGAL::Polygon_with_holes_2< TKernel >;

// =============================================================================
// BLOCK 1 — I/O  (same logic as main.cxx, only the TPoint type changes)
// =============================================================================

TPolygon read_polygon( const std::string& fname )
{
  std::ifstream ifs( fname );
  if( !ifs )
  {
    std::cerr << "Error: could not open " << fname << std::endl;
    return {};
  }

  std::vector< TPoint > P;
  std::vector< std::pair< std::size_t, std::size_t > > L;

  std::string line;
  while( std::getline( ifs, line ) )
  {
    if( line.empty( ) || line[0] == '#' ) continue;

    if( line[0] == 'v' )
    {
      double x, y;
      std::istringstream( line.substr( 1 ) ) >> x >> y;
      P.push_back( TPoint( x, y ) );
    }
    else if( line[0] == 'l' )
    {
      std::size_t a, b;
      std::istringstream( line.substr( 1 ) ) >> a >> b;
      // .obj uses 1-based indices
      L.push_back( { a - 1, b - 1 } );
    }
  }

  // Build Polygon_2 by pushing vertices in edge order
  TPolygon poly;
  for( const auto& edge : L )
    poly.push_back( P[ edge.first ] );

  return poly;
}

void save_polygon( const TPolygon& poly, const std::string& fname )
{
  std::ofstream ofs( fname );
  if( !ofs )
  {
    std::cerr << "Error: could not write " << fname << std::endl;
    return;
  }

  // Polygon_2 exposes vertex iterators directly — cleaner than the manual vector
  for( auto it = poly.vertices_begin( ); it != poly.vertices_end( ); ++it )
    ofs << "v " << CGAL::to_double( it->x( ) )
        << " "  << CGAL::to_double( it->y( ) ) << " 0\n";

  const std::size_t n = poly.size( );
  for( std::size_t i = 0; i < n; ++i )
    ofs << "l " << (i + 1) << " " << ((i + 1) % n + 1) << "\n";
}

// =============================================================================
// MAIN
// With arguments:    ./polygon_intersection_cgal P.obj Q.obj output.obj
// Without arguments: reads case 1 from data/samples/ by default
// =============================================================================

int main( int argc, char** argv )
{
  std::string fP, fQ, fOut;

  if( argc >= 4 )
  {
    fP   = argv[1];
    fQ   = argv[2];
    fOut = argv[3];
  }
  else
  {
    // Default path relative to the execution directory (build/)
    const std::string base = "../data/samples/";
    fP   = base + "caso1_parcial_P.obj";
    fQ   = base + "caso1_parcial_Q.obj";
    fOut = base + "result_cgal.obj";
    std::cout << "No arguments provided — using default sample:\n"
              << "  P:      " << fP   << "\n"
              << "  Q:      " << fQ   << "\n"
              << "  Output: " << fOut << "\n\n";
  }

  TPolygon P = read_polygon( fP );
  TPolygon Q = read_polygon( fQ );

  std::cout << "P has " << P.size( ) << " vertices" << std::endl;
  std::cout << "Q has " << Q.size( ) << " vertices" << std::endl;

  // CGAL::intersection returns a list of Polygon_with_holes_2.
  // For convex polygons the result is always a single polygon with no holes.
  std::vector< TPolygonWithHoles > result;
  CGAL::intersection( P, Q, std::back_inserter( result ) );

  std::cout << "Intersection has "
            << ( result.empty( ) ? 0 : result[0].outer_boundary( ).size( ) )
            << " vertices" << std::endl;

  if( !result.empty( ) )
    save_polygon( result[0].outer_boundary( ), fOut );
  else
    save_polygon( TPolygon{ }, fOut ); // no intersection -> empty file

  return EXIT_SUCCESS;
}