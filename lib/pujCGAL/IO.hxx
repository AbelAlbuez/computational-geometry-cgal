#ifndef __pujCGAL__IO__hxx__
#define __pujCGAL__IO__hxx__

#include <fstream>
#include <iterator>
#include <sstream>
#include <utility>
#include <vector>

// -------------------------------------------------------------------------
template< class TSegmentsIt >
bool pujCGAL::IO::
read( const std::string& fname, TSegmentsIt sIt )
{
  using TSegment = std::iter_value_t< typename TSegmentsIt::container_type >;
  using TKernel = typename TSegment::R;
  using TReal = typename TKernel::RT;
  using TPoint = typename TKernel::Point_2;
  using TPoints = std::vector< TPoint >;
  using TLine = std::pair< std::size_t, std::size_t >;
  using TLines = std::vector< TLine >;

  std::ifstream ifs( fname.c_str( ) );
  if( ifs )
  {
    TPoints P;
    TLines L;

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
        std::size_t a, b;
        std::istringstream( line.substr( 1 ) ) >> a >> b;
        L.push_back( TLine( a - 1, b - 1 ) );
      } // end if
    } // end while
    ifs.close( );

    for( const auto& l: L )
      *sIt = TSegment( P[ l.first ], P[ l.second ] );

    return( true );
  }
  else
    return( false );
}

// -------------------------------------------------------------------------
template< class TPointsIt >
bool pujCGAL::IO::
save( const std::string& fname, TPointsIt pB, TPointsIt pE )
{
  std::ofstream ofs( fname.c_str( ) );
  if( ofs )
  {
    for( auto pIt = pB; pIt != pE; ++pIt )
      ofs << "v " << pIt->x( ) << " " << pIt->y( ) << " 0" << std::endl;
    return( true );
  }
  else
    return( false );
}

#endif // __pujCGAL__IO__hxx__

// eof - IO.hxx
