// =============================================================================
// Workshop: Convex Polygon Intersection with CGAL
// Author: Abel Albuez Sanchez
//
// Algorithm: Sutherland-Hodgman
//   Idea: clip polygon P against each half-plane defined by the edges of Q.
//         The result is the intersection of both polygons.
//
// Usage: ./polygon_intersection polygon_P.obj polygon_Q.obj output.obj
// =============================================================================

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include <CGAL/Cartesian.h>

// Kernel types — same as the rest of the project
using TKernel  = CGAL::Cartesian< long double >;
using TReal    = TKernel::RT;
using TPoint   = TKernel::Point_2;
using TSegment = TKernel::Segment_2;
using TLine    = TKernel::Line_2;

// A polygon is an ordered list of vertices (CCW)
using TPolygon = std::vector< TPoint >;

// =============================================================================
// BLOCK 1 — I/O
// Reuses the logic from IO.hxx but builds an ordered vertex list instead of
// loose segments, since Sutherland-Hodgman requires vertex order.
// =============================================================================

TPolygon read_polygon( const std::string& fname )
{
  // Same parsing as IO.hxx, but we build a vector of points in the ORDER
  // the 'l' entries appear, preserving the CCW/CW orientation of the polygon.

  std::ifstream ifs( fname );
  if( !ifs )
  {
    std::cerr << "Error: could not open " << fname << std::endl;
    return {};
  }

  std::vector< TPoint > P;                                  // all vertices in the file
  std::vector< std::pair< std::size_t, std::size_t > > L;  // edges (i, j)

  std::string line;
  while( std::getline( ifs, line ) )
  {
    if( line.empty( ) || line[0] == '#' ) continue;

    if( line[0] == 'v' )
    {
      TReal x, y;
      std::istringstream ss( line.substr( 1 ) );
      ss >> x >> y;
      P.push_back( TPoint( x, y ) );
    }
    else if( line[0] == 'l' )
    {
      std::size_t a, b;
      std::istringstream ss( line.substr( 1 ) );
      ss >> a >> b;
      // .obj uses 1-based indices
      L.push_back( { a - 1, b - 1 } );
    }
  }

  // Rebuild the polygon following edge order.
  // Each edge l[i] = (a, b) gives point 'a' as vertex i of the polygon.
  TPolygon poly;
  for( const auto& edge : L )
    poly.push_back( P[ edge.first ] );

  return poly;
}

void save_polygon( const TPolygon& poly, const std::string& fname )
{
  // Same .obj format used by the project: 'v x y 0' and 'l i j' lines
  std::ofstream ofs( fname );
  if( !ofs )
  {
    std::cerr << "Error: could not write " << fname << std::endl;
    return;
  }

  // Vertices
  for( const auto& p : poly )
    ofs << "v " << p.x( ) << " " << p.y( ) << " 0\n";

  // Edges (close the polygon with the last edge n->1)
  const std::size_t n = poly.size( );
  for( std::size_t i = 0; i < n; ++i )
    ofs << "l " << (i + 1) << " " << ((i + 1) % n + 1) << "\n";
}

// =============================================================================
// BLOCK 2 — Geometric predicates
// These are the building blocks of the algorithm, as covered in Stage 2.
// =============================================================================

// Predicate: is P inside the half-plane defined by edge A->B?
// In a CCW polygon, the interior lies to the LEFT of each directed edge.
// TLine( A, B ) defines the directed line from A to B.
// has_on_positive_side(P) == true when P is to the LEFT (cross product > 0).
// Points ON the edge are also accepted (has_on_negative_side == false).
bool inside( const TPoint& P, const TPoint& A, const TPoint& B )
{
  TLine edge( A, B );
  return !edge.has_on_negative_side( P );
  // includes the boundary (COLLINEAR and LEFT_TURN)
}

// Predicate: intersection between segment P1->P2 and line A->B.
// Precondition: the segment is known to cross the line (guaranteed by Sutherland-Hodgman).
TPoint intersect_edge( const TPoint& P1, const TPoint& P2,
                       const TPoint& A,  const TPoint& B )
{
  TSegment seg( P1, P2 );
  TLine    line( A, B );

  // CGAL returns std::optional<std::variant<Point_2, Segment_2>>
  auto result = CGAL::intersection( seg, line );
  if( result )
    if( const TPoint* p = std::get_if< TPoint >( &*result ) )
      return *p;

  // Should never reach here when called correctly
  // (only invoked when a point intersection is guaranteed)
  std::cerr << "WARNING: intersect_edge — no intersection point found\n";
  return P1;
}

// =============================================================================
// BLOCK 3 — Sutherland-Hodgman
// Clips polygon 'subject' against ONE half-plane defined by edge A->B.
// Called once per edge of the clip polygon Q.
// =============================================================================

TPolygon clip_by_edge( const TPolygon& subject,
                       const TPoint& A, const TPoint& B )
{
  TPolygon output;
  if( subject.empty( ) ) return output;

  const std::size_t n = subject.size( );

  for( std::size_t i = 0; i < n; ++i )
  {
    // Current edge of subject: current -> next (wraps around with % n)
    const TPoint& current = subject[ i ];
    const TPoint& next    = subject[ (i + 1) % n ];

    bool current_inside = inside( current, A, B );
    bool next_inside    = inside( next,    A, B );

    // ┌──────────────────┬──────────────────────────────────────────────────┐
    // │       CASE        │  Action                                          │
    // ├──────────────────┼──────────────────────────────────────────────────┤
    // │ inside -> inside  │ add 'next'                                       │
    // │ inside -> outside │ add intersection with clip edge                  │
    // │ outside -> inside │ add intersection + 'next'                        │
    // │ outside -> outside│ add nothing                                      │
    // └──────────────────┴──────────────────────────────────────────────────┘

    if( current_inside && next_inside )
    {
      // Case 1: both inside -> only add 'next'
      // ('current' was already added in the previous iteration)
      output.push_back( next );
    }
    else if( current_inside && !next_inside )
    {
      // Case 2: exiting clip polygon -> add the crossing point
      output.push_back( intersect_edge( current, next, A, B ) );
    }
    else if( !current_inside && next_inside )
    {
      // Case 3: entering clip polygon -> add crossing point + 'next'
      output.push_back( intersect_edge( current, next, A, B ) );
      output.push_back( next );
    }
    // Case 4: both outside -> add nothing
  }

  return output;
}

// =============================================================================
// BLOCK 4 — Main algorithm
// Applies clip_by_edge once per edge of polygon Q.
// =============================================================================

TPolygon polygon_intersection( const TPolygon& P, const TPolygon& Q )
{
  // Start with P as the subject polygon to be clipped
  TPolygon result = P;

  const std::size_t m = Q.size( );

  // For each edge A->B of Q, clip the accumulated result
  for( std::size_t i = 0; i < m; ++i )
  {
    if( result.empty( ) ) break; // no intersection — early exit

    const TPoint& A = Q[ i ];
    const TPoint& B = Q[ (i + 1) % m ];

    result = clip_by_edge( result, A, B );
  }

  return result;
}

// =============================================================================
// MAIN
// With arguments:    ./polygon_intersection P.obj Q.obj output.obj
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
    fOut = base + "result_sh.obj";
    std::cout << "No arguments provided — using default sample:\n"
              << "  P:      " << fP   << "\n"
              << "  Q:      " << fQ   << "\n"
              << "  Output: " << fOut << "\n\n";
  }

  TPolygon P = read_polygon( fP );
  TPolygon Q = read_polygon( fQ );

  std::cout << "P has " << P.size( ) << " vertices" << std::endl;
  std::cout << "Q has " << Q.size( ) << " vertices" << std::endl;

  TPolygon result = polygon_intersection( P, Q );

  std::cout << "Intersection has " << result.size( ) << " vertices" << std::endl;

  save_polygon( result, fOut );

  return EXIT_SUCCESS;
}