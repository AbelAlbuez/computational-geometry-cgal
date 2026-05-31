#ifndef __ContourResampler__hxx__
#define __ContourResampler__hxx__

#include <cmath>
#include <vector>

// -------------------------------------------------------------------------
inline pujCGAL::Final::ContourResampler::TContour
pujCGAL::Final::ContourResampler::
resample( const TContour& contour, int n )
{
  TContour out;
  const std::size_t m = contour.size( );
  if( m == 0 || n <= 0 )
    return out;
  if( m == 1 )
  {
    out.assign( static_cast< std::size_t >( n ), contour[ 0 ] );
    return out;
  }

  // -- Per-vertex cumulative arc length (circular: index m corresponds
  //    to the closing edge back to vertex 0).
  std::vector< double > acc( m + 1, 0.0 );
  for( std::size_t i = 0; i < m; ++i )
  {
    const TPoint& a = contour[ i ];
    const TPoint& b = contour[ ( i + 1 ) % m ];
    const double dx = CGAL::to_double( b.x( ) ) - CGAL::to_double( a.x( ) );
    const double dy = CGAL::to_double( b.y( ) ) - CGAL::to_double( a.y( ) );
    acc[ i + 1 ] = acc[ i ] + std::sqrt( dx * dx + dy * dy );
  }
  const double total = acc[ m ];
  if( total <= 0.0 )
  {
    out.assign( static_cast< std::size_t >( n ), contour[ 0 ] );
    return out;
  }

  out.reserve( static_cast< std::size_t >( n ) );
  const double step = total / static_cast< double >( n );
  std::size_t seg = 0;
  for( int k = 0; k < n; ++k )
  {
    const double s = static_cast< double >( k ) * step;

    // Advance to the segment containing arc length `s`.
    while( seg < m && acc[ seg + 1 ] < s )
      ++seg;
    if( seg >= m ) seg = m - 1;

    const double seg_len = acc[ seg + 1 ] - acc[ seg ];
    const double u = ( seg_len > 0.0 ) ? ( s - acc[ seg ] ) / seg_len : 0.0;

    const TPoint& a = contour[ seg ];
    const TPoint& b = contour[ ( seg + 1 ) % m ];
    const double x = ( 1.0 - u ) * CGAL::to_double( a.x( ) )
                   + u * CGAL::to_double( b.x( ) );
    const double y = ( 1.0 - u ) * CGAL::to_double( a.y( ) )
                   + u * CGAL::to_double( b.y( ) );
    out.emplace_back( x, y );
  }
  return out;
}

#endif // __ContourResampler__hxx__

// eof - ContourResampler.hxx
