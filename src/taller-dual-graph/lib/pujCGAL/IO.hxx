#ifndef __pujCGAL__IO__hxx__
#define __pujCGAL__IO__hxx__

#include <fstream>
#include <iterator>
#include <sstream>
#include <utility>
#include <vector>

#include <pujCGAL/Polygon.h>
#include <pujCGAL/Triangulation.h>

// -------------------------------------------------------------------------
template< class TKernel >
bool pujCGAL::IO::
read( const std::string& fname, pujCGAL::Polygon< TKernel >& polygon )
{
  using TReal    = typename TKernel::RT;
  using TPolygon = pujCGAL::Polygon< TKernel >;
  using TPoint   = typename TPolygon::TPoint;
  using TPoints  = std::vector< TPoint >;
  using TLines   = std::vector< std::size_t >;
  
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
      else if( line[ 0 ] == 'l' || line[ 0 ] == 'f' )
      {
        std::size_t a;
        std::istringstream line_str( line.substr( 1 ) );
        line_str >> a;
        while( line_str )
        {
          L.push_back( a - 1 );
          line_str >> a;
        } // end if
      } // end if
    } // end while
    ifs.close( );

    polygon.add_point( P[ L.front( ) ] );
    for( const auto& l: L )
      if( l != L.front( ) )
        polygon.add_point( P[ l ] );
    polygon.build( );

    return( true );
   }
   else
   return( false );
}

/* TODO
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
*/

// -------------------------------------------------------------------------
template< class TKernel >
bool pujCGAL::IO::
save( const std::string& fname, const pujCGAL::Triangulation< TKernel >& mesh )
{
  std::ofstream ofs( fname.c_str( ) );
  if( ofs )
  {
    for(
      auto pIt = mesh.geometry_begin( ); pIt != mesh.geometry_end( ); ++pIt
      )
      ofs << "v "
          << ( *pIt )[ 0 ] << " "
          << ( *pIt )[ 1 ] << " 0" << std::endl;

    ofs << std::endl;

    for(
      auto tIt = mesh.topology_begin( ); tIt != mesh.topology_end( ); ++tIt
      )
      ofs << "f "
          << ( ( *tIt )[ 0 ] + 1 ) << " "
          << ( ( *tIt )[ 1 ] + 1 ) << " "
          << ( ( *tIt )[ 2 ] + 1 ) << std::endl;

    ofs << std::endl;

    ofs.close( );
    return( true );
  }
  else
    return( false );
}


#endif // __pujCGAL__IO__hxx__

// eof - IO.hxx
