
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

#include <CGAL/Cartesian.h>
#include <CGAL/intersections.h>

using TKernel = CGAL::Cartesian< long double >;
using TReal = TKernel::RT;
using TPoint = TKernel::Point_2;
using TSegment = TKernel::Segment_2;
using TSegments = std::vector< TSegment >;
using TPoints = std::vector< TPoint >;

bool read( TSegments& segments, const std::string& fname )
{
  std::ifstream ifs( fname.c_str( ) );
  if( ifs )
  {
    TPoints P;
    std::vector< std::pair< std::size_t, std::size_t > > L;

    std::string line;
    while( std::getline( ifs, line ) )
    {
      if( line[ 0 ] == 'v' )
      {
        TReal x, y;
        std::istringstream( line.substr( 1 ) ) >> x >> y;
        P.push_back( TPoint( x, y ) );
      }
      else if( line[ 0 ] == 'l' )
      {
        std::pair< std::size_t, std::size_t > l;
        std::istringstream( line.substr( 1 ) ) >> l.first >> l.second;
        L.push_back( l );
      } // end if
    } // end while
    ifs.close( );

    for( const auto& l: L )
      segments.push_back( TSegment( P[ l.first - 1 ], P[ l.second - 1 ] ) );

    return( true );
  }
  else
    return( false );
}

bool save( const TPoints& points, const std::string& fname )
{
  std::ofstream ofs( fname.c_str( ) );
  if( ofs )
  {
    for( const auto& p: points )
      ofs << "v " << p[ 0 ] << " " << p[ 1 ] << " 0" << std::endl;
    return( true );
  }
  else
    return( false );
}

std::pair< bool, TPoint > intersect( const TSegment& a, const TSegment& b )
{
  auto r = CGAL::intersection( a, b );
  if( r )
  {
    if( const TPoint* p = std::get_if< TPoint >( &*r ) )
      return( std::make_pair( true, *p ) );
    else
      return( std::make_pair( false, TPoint( ) ) );
  }
  else
    return( std::make_pair( false, TPoint( ) ) );
}

void intersections( TPoints& I, const TSegments& S )
{
  for( auto i = S.begin( ); i != S.end( ); ++i )
  {
    for( auto j = i; j != S.end( ); ++j )
    {
      auto id = intersect( *i, *j );
      if( id.first )
        I.push_back( id.second );
    } // end for
  } // end for
}

int main( int argc, char** argv )
{
  TSegments S;
  TPoints I;

  read( S, argv[ 1 ] );
  intersections( I, S );
  save( I, argv[ 2 ] );

  return( EXIT_SUCCESS );
}

// eof - brute_force.cxx
