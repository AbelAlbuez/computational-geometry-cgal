#ifndef __LinearInterpolator__hxx__
#define __LinearInterpolator__hxx__

// -------------------------------------------------------------------------
inline pujCGAL::Final::LinearInterpolator::TContour
pujCGAL::Final::LinearInterpolator::
interpolate( const TContour& A, const TContour& B, double t )
{
  TContour out;
  const std::size_t n = A.size( );
  out.reserve( n );

  const double s = 1.0 - t;
  for( std::size_t i = 0; i < n; ++i )
  {
    const double ax = CGAL::to_double( A[ i ].x( ) );
    const double ay = CGAL::to_double( A[ i ].y( ) );
    const double bx = CGAL::to_double( B[ i ].x( ) );
    const double by = CGAL::to_double( B[ i ].y( ) );
    out.emplace_back( s * ax + t * bx, s * ay + t * by );
  }
  return out;
}

#endif // __LinearInterpolator__hxx__

// eof - LinearInterpolator.hxx
