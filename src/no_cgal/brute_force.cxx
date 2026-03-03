#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>
#include <tuple>

using TReal = long double;
using TPoint = std::pair< TReal, TReal >;
using TSegment = std::pair< TPoint, TPoint >;
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
        TPoint p;
        std::istringstream( line.substr( 1 ) ) >> p.first >> p.second;
        P.push_back( p );
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
      ofs << "v " << p.first << " " << p.second << " 0" << std::endl;
    return( true );
  }
  else
    return( false );
}

std::tuple< TReal, TReal, TPoint > intersect(
  const TSegment& a, const TSegment& b
  )
{
  TReal x1 = a.first.first;
  TReal y1 = a.first.second;
  TReal x2 = a.second.first;
  TReal y2 = a.second.second;
  TReal x3 = b.first.first;
  TReal y3 = b.first.second;
  TReal x4 = b.second.first;
  TReal y4 = b.second.second;

  TReal t = std::numeric_limits< TReal >::quiet_NaN( );
  TReal u = t;
  TPoint i( t, t );

  TReal d = ( ( x1 - x2 ) * ( y3 - y4 ) ) - ( ( y1 - y2 ) * ( x3 - x4 ) );
  if( d != TReal( 0 ) )
  {
    t = ( ( ( x1 - x3 ) * ( y3 - y4 ) ) - ( ( y1 - y3 ) * ( x3 - x4 ) ) ) / d;
    u = ( ( ( y1 - y2 ) * ( x1 - x3 ) ) - ( ( x1 - x2 ) * ( y1 - y3 ) ) ) / d;
    i.first  = x1 + ( t * ( x2 - x1 ) );
    i.second = y1 + ( t * ( y2 - y1 ) );
  } // end if

  return( std::make_tuple( t, u, i ) );
}

void intersections( TPoints& I, const TSegments& S )
{
  for( auto i = S.begin( ); i != S.end( ); ++i )
  {
    for( auto j = i; j != S.end( ); ++j )
    {
      auto id = intersect( *i, *j );
      if(
        !(
          std::isnan( std::get< 0 >( id ) )
          ||
          std::isnan( std::get< 1 >( id ) )
          )
        &&
        TReal( 0 ) <= std::get< 0 >( id ) && std::get< 0 >( id ) <= TReal( 1 )
        &&
        TReal( 0 ) <= std::get< 1 >( id ) && std::get< 1 >( id ) <= TReal( 1 )
        )
        I.push_back( std::get< 2 >( id ) );
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
