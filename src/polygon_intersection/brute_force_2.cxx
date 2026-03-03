#include <iostream>
#include <iterator>
#include <vector>

#include <CGAL/Cartesian.h>
#include <CGAL/intersections.h>

#include <pujCGAL/IO.h>

template< class TSegmentsIt, class TPointsIt >
void intersections( TSegmentsIt sB, TSegmentsIt sE, TPointsIt pIt )
{
  using TPoint = std::iter_value_t< typename TPointsIt::container_type >;

  for( auto i = sB; i != sE; ++i )
    for( auto j = i; j != sE; ++j )
      if( auto r = CGAL::intersection( *i, *j ) )
        if( const TPoint* p = std::get_if< TPoint >( &*r ) )
          *pIt = *p;
}

int main( int argc, char** argv )
{
  using TKernel = CGAL::Cartesian< long double >;
  using TPoint = TKernel::Point_2;
  using TSegment = TKernel::Segment_2;

  std::vector< TSegment > S;
  std::vector< TPoint > I;

  pujCGAL::IO::read( argv[ 1 ], std::back_inserter( S ) );
  intersections( S.begin( ), S.end( ), std::back_inserter( I ) );
  pujCGAL::IO::save( argv[ 2 ], I.begin( ), I.end( ) );

  return( EXIT_SUCCESS );
}

// eof - brute_force_2.cxx
