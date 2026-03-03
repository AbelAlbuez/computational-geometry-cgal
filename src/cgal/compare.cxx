#include <iostream>
#include <iterator>
#include <vector>

#include <CGAL/Cartesian.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <pujCGAL/IO.h>
#include <pujCGAL/SegmentsIntersection.h>

int main( int argc, char** argv )
{
  using TKernel = CGAL::Cartesian< long double >;
  // using TKernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using TPoint = TKernel::Point_2;
  using TSegment = TKernel::Segment_2;

  std::vector< TSegment > S;
  std::vector< TPoint > BF, BO;

  pujCGAL::IO::read( argv[ 1 ], std::back_inserter( S ) );
  pujCGAL::SegmentsIntersection::BruteForce(
    S.begin( ), S.end( ), std::back_inserter( BF )
    );
  pujCGAL::SegmentsIntersection::BentleyOttmann(
    S.begin( ), S.end( ), std::back_inserter( BO )
    );

  pujCGAL::IO::save( argv[ 2 ], BF.begin( ), BF.end( ) );
  pujCGAL::IO::save( argv[ 3 ], BO.begin( ), BO.end( ) );

  return( EXIT_SUCCESS );
}

// eof - compare.cxx
