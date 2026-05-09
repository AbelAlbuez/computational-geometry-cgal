#ifndef __ContourInterpolator__hxx__
#define __ContourInterpolator__hxx__

#include <fstream>
#include <iostream>
#include <sstream>

// -------------------------------------------------------------------------
inline pujCGAL::Final::ContourInterpolator::TContour
pujCGAL::Final::ContourInterpolator::
read_obj( const std::string& fname )
{
  TContour contour;
  std::ifstream ifs( fname );
  if( !ifs )
  {
    std::cerr << "Error: could not open " << fname << std::endl;
    return contour;
  }

  std::string line;
  while( std::getline( ifs, line ) )
  {
    if( line.empty( ) || line[ 0 ] == '#' ) continue;
    if( line[ 0 ] != 'v' ) continue;

    std::istringstream ss( line.substr( 1 ) );
    double x, y;
    if( ss >> x >> y )
      contour.emplace_back( x, y );
  } // end while
  return contour;
}

// -------------------------------------------------------------------------
inline void
pujCGAL::Final::ContourInterpolator::
write_obj( const std::string& fname, const TContour& contour )
{
  std::ofstream ofs( fname );
  if( !ofs )
  {
    std::cerr << "Error: could not write " << fname << std::endl;
    return;
  }

  ofs << "# Contorno interpolado - " << contour.size( ) << " vertices\n";
  for( const auto& p : contour )
    ofs << "v " << p.x( ) << " " << p.y( ) << "\n";

  const std::size_t n = contour.size( );
  for( std::size_t i = 1; i < n; ++i )
    ofs << "l " << i << " " << ( i + 1 ) << "\n";
  if( n >= 2 )
    ofs << "l " << n << " 1\n";
}

// -------------------------------------------------------------------------
inline pujCGAL::Final::ContourInterpolator::TContour
pujCGAL::Final::ContourInterpolator::
interpolate( const TContour& A, const TContour& B, double /*t*/ )
{
  // TODO: implementación geométrica (correspondencia + interpolación).
  // Por ahora regresa A para no romper el flujo.
  ( void ) B;
  return A;
}

#endif // __ContourInterpolator__hxx__

// eof - ContourInterpolator.hxx
