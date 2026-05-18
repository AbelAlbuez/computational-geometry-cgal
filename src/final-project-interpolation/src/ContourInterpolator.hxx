#ifndef __ContourInterpolator__hxx__
#define __ContourInterpolator__hxx__

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <vector>

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
    ofs << "v " << p.x( ) << " " << p.y( ) << " 0.000000\n";

  const std::size_t n = contour.size( );
  for( std::size_t i = 1; i < n; ++i )
    ofs << "l " << i << " " << ( i + 1 ) << "\n";
  if( n >= 2 )
    ofs << "l " << n << " 1\n";
}

// -------------------------------------------------------------------------
// Funciones auxiliares para la correspondencia de vértices por ángulo polar
// -------------------------------------------------------------------------
inline pujCGAL::Final::ContourInterpolator::TPoint
calcular_centroide( const pujCGAL::Final::ContourInterpolator::TContour& contorno )
{
  double cx = 0.0, cy = 0.0;
  for( const auto& p : contorno )
  {
    cx += p.x( );
    cy += p.y( );
  }
  const double n = static_cast< double >( contorno.size( ) );
  return pujCGAL::Final::ContourInterpolator::TPoint( cx / n, cy / n );
}

inline pujCGAL::Final::ContourInterpolator::TContour
ordenar_por_angulo( const pujCGAL::Final::ContourInterpolator::TContour& contorno )
{
  if( contorno.empty( ) ) return contorno;

  auto centroide = calcular_centroide( contorno );
  pujCGAL::Final::ContourInterpolator::TContour ordenado = contorno;

  std::sort( ordenado.begin( ), ordenado.end( ),
             [ &centroide ]( const pujCGAL::Final::ContourInterpolator::TPoint& p1,
                             const pujCGAL::Final::ContourInterpolator::TPoint& p2 )
             {
               double ang1 = std::atan2( p1.y( ) - centroide.y( ), p1.x( ) - centroide.x( ) );
               double ang2 = std::atan2( p2.y( ) - centroide.y( ), p2.x( ) - centroide.x( ) );
               return ang1 < ang2;
             } );

  return ordenado;
}

// -------------------------------------------------------------------------
inline pujCGAL::Final::ContourInterpolator::TContour
pujCGAL::Final::ContourInterpolator::
interpolate( const TContour& A, const TContour& B, double t )
{
  if( A.empty( ) || B.empty( ) ) return A;

  // Paso 1: Ordenar ambos contornos por ángulo polar respecto a sus centroides
  TContour ordenadoA = ordenar_por_angulo( A );
  TContour ordenadoB = ordenar_por_angulo( B );

  TContour resultado;
  resultado.reserve( ordenadoA.size( ) );

  const std::size_t nA = ordenadoA.size( );
  const std::size_t nB = ordenadoB.size( );

  // Paso 2: Generar la correspondencia e interpolación lineal entre pares
  for( std::size_t i = 0; i < nA; ++i )
  {
    // Mapeo indexado proporcional por si los contornos poseen diferente número de vértices
    std::size_t idxB = ( i * nB ) / nA;
    if( idxB >= nB ) idxB = nB - 1;

    const auto& pA = ordenadoA[ i ];
    const auto& pB = ordenadoB[ idxB ];

    // Interpolación lineal paramétrica
    double x = pA.x( ) + t * ( pB.x( ) - pA.x( ) );
    double y = pA.y( ) + t * ( pB.y( ) - pA.y( ) );

    resultado.emplace_back( x, y );
  }

  return resultado;
}

#endif // __ContourInterpolator__hxx__