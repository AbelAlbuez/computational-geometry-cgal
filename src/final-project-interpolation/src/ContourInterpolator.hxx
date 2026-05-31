#ifndef __ContourInterpolator__hxx__
#define __ContourInterpolator__hxx__

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// =============================================================================
// I/O
// =============================================================================
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

// =============================================================================
// Geometry helpers
// =============================================================================
inline double
pujCGAL::Final::ContourInterpolator::
signed_area( const TContour& c )
{
  const std::size_t n = c.size( );
  if( n < 3 ) return 0.0;
  double a = 0.0;
  for( std::size_t i = 0; i < n; ++i )
  {
    const TPoint& p = c[ i ];
    const TPoint& q = c[ ( i + 1 ) % n ];
    a += p.x( ) * q.y( ) - q.x( ) * p.y( );
  }
  return 0.5 * a;
}

// -------------------------------------------------------------------------
inline pujCGAL::Final::ContourInterpolator::TPoint
pujCGAL::Final::ContourInterpolator::
centroid( const TContour& c )
{
  if( c.empty( ) ) return TPoint( 0, 0 );
  double sx = 0, sy = 0;
  for( const auto& p : c ) { sx += p.x( ); sy += p.y( ); }
  const double n = static_cast< double >( c.size( ) );
  return TPoint( sx / n, sy / n );
}

// -------------------------------------------------------------------------
inline double
pujCGAL::Final::ContourInterpolator::
perimeter( const TContour& c )
{
  const std::size_t n = c.size( );
  if( n < 2 ) return 0.0;
  double L = 0.0;
  for( std::size_t i = 0; i < n; ++i )
  {
    const double dx = c[ ( i + 1 ) % n ].x( ) - c[ i ].x( );
    const double dy = c[ ( i + 1 ) % n ].y( ) - c[ i ].y( );
    L += std::hypot( dx, dy );
  }
  return L;
}

// -------------------------------------------------------------------------
inline void
pujCGAL::Final::ContourInterpolator::
ensure_ccw( TContour& c )
{
  // Shoelace sign is the orientation; reverse a clockwise loop so that
  // correspondence indices and (later) outward normals are consistent.
  if( signed_area( c ) < 0.0 )
    std::reverse( c.begin( ), c.end( ) );
}

// -------------------------------------------------------------------------
inline pujCGAL::Final::ContourInterpolator::TContour
pujCGAL::Final::ContourInterpolator::
resample_uniform( const TContour& c, std::size_t N )
{
  const std::size_t n = c.size( );
  if( n < 2 || N == 0 ) return c;

  std::vector< double > elen( n );
  double L = 0.0;
  for( std::size_t i = 0; i < n; ++i )
  {
    const double dx = c[ ( i + 1 ) % n ].x( ) - c[ i ].x( );
    const double dy = c[ ( i + 1 ) % n ].y( ) - c[ i ].y( );
    elen[ i ] = std::hypot( dx, dy );
    L += elen[ i ];
  }
  if( L <= 0.0 ) return c;

  TContour out;
  out.reserve( N );
  const double step = L / static_cast< double >( N );
  std::size_t e = 0;     // current edge index
  double acc = 0.0;      // arc length at the start of edge e
  for( std::size_t j = 0; j < N; ++j )
  {
    const double target = j * step;
    while( e + 1 < n && acc + elen[ e ] < target ) { acc += elen[ e ]; ++e; }
    const double rem = target - acc;
    const double f = ( elen[ e ] > 0.0 ) ? ( rem / elen[ e ] ) : 0.0;
    const TPoint& p0 = c[ e ];
    const TPoint& p1 = c[ ( e + 1 ) % n ];
    out.emplace_back( p0.x( ) + f * ( p1.x( ) - p0.x( ) ),
                      p0.y( ) + f * ( p1.y( ) - p0.y( ) ) );
  }
  return out;
}

// -------------------------------------------------------------------------
inline pujCGAL::Final::ContourInterpolator::TContour
pujCGAL::Final::ContourInterpolator::
resample_adaptive( const TContour& c, std::size_t N, double lambda )
{
  const std::size_t n = c.size( );
  if( n < 3 || N == 0 ) return resample_uniform( c, N );

  // Turning-angle curvature at each vertex: kappa_i = 2*theta_i / (|u|+|v|).
  std::vector< double > kappa( n, 0.0 );
  for( std::size_t i = 0; i < n; ++i )
  {
    const TPoint& pm = c[ ( i + n - 1 ) % n ];
    const TPoint& p  = c[ i ];
    const TPoint& pp = c[ ( i + 1 ) % n ];
    const double ux = p.x( ) - pm.x( ), uy = p.y( ) - pm.y( );
    const double vx = pp.x( ) - p.x( ), vy = pp.y( ) - p.y( );
    const double lu = std::hypot( ux, uy ), lv = std::hypot( vx, vy );
    if( lu <= 0.0 || lv <= 0.0 ) continue;
    const double cross = ux * vy - uy * vx;
    const double dot   = ux * vx + uy * vy;
    const double theta = std::atan2( cross, dot ); // signed turning angle
    kappa[ i ] = 2.0 * std::fabs( theta ) / ( lu + lv );
  }

  // Warp each edge length by its average density (1 + lambda*|kappa|) and
  // sample uniformly in this warped parameter -> denser near high curvature.
  std::vector< double > w( n );
  double W = 0.0;
  for( std::size_t i = 0; i < n; ++i )
  {
    const double dx = c[ ( i + 1 ) % n ].x( ) - c[ i ].x( );
    const double dy = c[ ( i + 1 ) % n ].y( ) - c[ i ].y( );
    const double elen = std::hypot( dx, dy );
    const double dens = 1.0 + lambda * 0.5 * ( kappa[ i ] + kappa[ ( i + 1 ) % n ] );
    w[ i ] = elen * dens;
    W += w[ i ];
  }
  if( W <= 0.0 ) return resample_uniform( c, N );

  TContour out;
  out.reserve( N );
  const double step = W / static_cast< double >( N );
  std::size_t e = 0;
  double acc = 0.0;
  for( std::size_t j = 0; j < N; ++j )
  {
    const double target = j * step;
    while( e + 1 < n && acc + w[ e ] < target ) { acc += w[ e ]; ++e; }
    const double rem = target - acc;
    const double f = ( w[ e ] > 0.0 ) ? ( rem / w[ e ] ) : 0.0;
    const TPoint& p0 = c[ e ];
    const TPoint& p1 = c[ ( e + 1 ) % n ];
    out.emplace_back( p0.x( ) + f * ( p1.x( ) - p0.x( ) ),
                      p0.y( ) + f * ( p1.y( ) - p0.y( ) ) );
  }
  return out;
}

// -------------------------------------------------------------------------
// Cyclic alignment: minimize sum_i ||A_i - B_{(i+sigma) mod N}||^2 over the
// shift sigma, on centroid-centered copies (so the cost reflects shape, not
// position). Also tries the reversed traversal of B. O(N^2), N small.
inline void
pujCGAL::Final::ContourInterpolator::
cyclic_align( const TContour& A, const TContour& B,
              std::size_t& sigma, bool& reversed )
{
  sigma = 0; reversed = false;
  const std::size_t N = A.size( );
  if( N == 0 || B.size( ) != N ) return;

  const TPoint cA = centroid( A );
  const TPoint cB = centroid( B );
  std::vector< double > ax( N ), ay( N ), bx( N ), by( N );
  for( std::size_t i = 0; i < N; ++i )
  {
    ax[ i ] = A[ i ].x( ) - cA.x( );  ay[ i ] = A[ i ].y( ) - cA.y( );
    bx[ i ] = B[ i ].x( ) - cB.x( );  by[ i ] = B[ i ].y( ) - cB.y( );
  }

  double best = std::numeric_limits< double >::max( );
  for( int rev = 0; rev < 2; ++rev )
    for( std::size_t s = 0; s < N; ++s )
    {
      double cost = 0.0;
      for( std::size_t i = 0; i < N; ++i )
      {
        const std::size_t k  = ( i + s ) % N;
        const std::size_t bi = rev ? ( N - 1 - k ) : k;
        const double dx = ax[ i ] - bx[ bi ];
        const double dy = ay[ i ] - by[ bi ];
        cost += dx * dx + dy * dy;
        if( cost >= best ) break;                 // early prune
      }
      if( cost < best ) { best = cost; sigma = s; reversed = ( rev == 1 ); }
    }
}

namespace pujCGAL { namespace Final { namespace ci_detail
{
  // Reindex B (length N) so that result[i] corresponds to A[i].
  inline ContourInterpolator::TContour
  apply_shift( const ContourInterpolator::TContour& B,
               std::size_t sigma, bool reversed )
  {
    const std::size_t N = B.size( );
    ContourInterpolator::TContour out( N );
    for( std::size_t i = 0; i < N; ++i )
    {
      const std::size_t k  = ( i + sigma ) % N;
      const std::size_t bi = reversed ? ( N - 1 - k ) : k;
      out[ i ] = B[ bi ];
    }
    return out;
  }

  // Natural cubic spline second derivatives (M) via the Thomas algorithm.
  inline std::vector< double >
  spline_M( const std::vector< double >& z, const std::vector< double >& y )
  {
    const int n = (int) z.size( );
    std::vector< double > M( n, 0.0 );
    if( n < 3 ) return M;                          // natural ends => M=0
    std::vector< double > b( n, 0.0 ), c( n, 0.0 ), d( n, 0.0 );
    for( int i = 1; i < n - 1; ++i )
    {
      const double h0 = z[ i ] - z[ i - 1 ];
      const double h1 = z[ i + 1 ] - z[ i ];
      b[ i ] = 2.0 * ( h0 + h1 );
      c[ i ] = h1;
      d[ i ] = 6.0 * ( ( y[ i + 1 ] - y[ i ] ) / h1
                     - ( y[ i ]     - y[ i - 1 ] ) / h0 );
      // sub-diagonal a[i] = h0 used in the elimination below.
    }
    // Forward elimination (rows 1..n-2; M[0]=M[n-1]=0).
    for( int i = 2; i < n - 1; ++i )
    {
      const double a_i = z[ i ] - z[ i - 1 ];      // sub-diagonal
      const double m = a_i / b[ i - 1 ];
      b[ i ] -= m * c[ i - 1 ];
      d[ i ] -= m * d[ i - 1 ];
    }
    if( n - 2 >= 1 ) M[ n - 2 ] = d[ n - 2 ] / b[ n - 2 ];
    for( int i = n - 3; i >= 1; --i )
      M[ i ] = ( d[ i ] - c[ i ] * M[ i + 1 ] ) / b[ i ];
    return M;
  }

  inline double
  spline_eval( const std::vector< double >& z, const std::vector< double >& y,
               const std::vector< double >& M, double zq )
  {
    const int n = (int) z.size( );
    if( n == 1 ) return y[ 0 ];
    int k = 0;
    if( zq <= z.front( ) )      k = 0;
    else if( zq >= z.back( ) )  k = n - 2;
    else { while( k < n - 2 && z[ k + 1 ] < zq ) ++k; }
    const double h = z[ k + 1 ] - z[ k ];
    const double A = ( z[ k + 1 ] - zq ) / h;
    const double B = ( zq - z[ k ] ) / h;
    return A * y[ k ] + B * y[ k + 1 ]
         + ( ( A * A * A - A ) * M[ k ] + ( B * B * B - B ) * M[ k + 1 ] )
           * ( h * h ) / 6.0;
  }

  // Least-squares polynomial fit of degree d via the normal equations.
  inline std::vector< double >
  polyfit( const std::vector< double >& x, const std::vector< double >& y, int d )
  {
    const int n = (int) x.size( );
    if( d > n - 1 ) d = n - 1;
    if( d < 0 ) d = 0;
    const int m = d + 1;
    std::vector< double > S( 2 * d + 1, 0.0 );     // power sums of x
    for( int i = 0; i < n; ++i )
    {
      double p = 1.0;
      for( int k = 0; k <= 2 * d; ++k ) { S[ k ] += p; p *= x[ i ]; }
    }
    std::vector< std::vector< double > > Aug( m, std::vector< double >( m + 1, 0.0 ) );
    for( int r = 0; r < m; ++r )
    {
      for( int c = 0; c < m; ++c ) Aug[ r ][ c ] = S[ r + c ];
      double t = 0.0;
      for( int i = 0; i < n; ++i )
      {
        double p = 1.0;
        for( int k = 0; k < r; ++k ) p *= x[ i ];
        t += y[ i ] * p;
      }
      Aug[ r ][ m ] = t;
    }
    // Gaussian elimination with partial pivoting.
    for( int col = 0; col < m; ++col )
    {
      int piv = col;
      for( int r = col + 1; r < m; ++r )
        if( std::fabs( Aug[ r ][ col ] ) > std::fabs( Aug[ piv ][ col ] ) ) piv = r;
      std::swap( Aug[ col ], Aug[ piv ] );
      const double diag = Aug[ col ][ col ];
      if( std::fabs( diag ) < 1e-300 ) continue;
      for( int r = 0; r < m; ++r )
      {
        if( r == col ) continue;
        const double f = Aug[ r ][ col ] / diag;
        for( int c = col; c <= m; ++c ) Aug[ r ][ c ] -= f * Aug[ col ][ c ];
      }
    }
    std::vector< double > coef( m, 0.0 );
    for( int r = 0; r < m; ++r )
      if( std::fabs( Aug[ r ][ r ] ) > 1e-300 ) coef[ r ] = Aug[ r ][ m ] / Aug[ r ][ r ];
    return coef;                                    // coef[k] is the x^k term
  }

  inline double polyeval( const std::vector< double >& coef, double x )
  {
    double v = 0.0;
    for( int k = (int) coef.size( ) - 1; k >= 0; --k ) v = v * x + coef[ k ];
    return v;
  }

  inline double lin_eval( const std::vector< double >& z,
                          const std::vector< double >& y, double zq )
  {
    const int n = (int) z.size( );
    if( n == 1 ) return y[ 0 ];
    if( zq <= z.front( ) ) return y.front( );
    if( zq >= z.back( ) )  return y.back( );
    int k = 0;
    while( k < n - 2 && z[ k + 1 ] < zq ) ++k;
    const double f = ( zq - z[ k ] ) / ( z[ k + 1 ] - z[ k ] );
    return y[ k ] + f * ( y[ k + 1 ] - y[ k ] );
  }
} } } // namespace pujCGAL::Final::ci_detail

// =============================================================================
// M1 - linear (vertex-correspondence LERP)
// =============================================================================
inline pujCGAL::Final::ContourInterpolator::TContour
pujCGAL::Final::ContourInterpolator::
interpolate( const TContour& A0, const TContour& B0, double t,
             bool curvature_adaptive, double lambda )
{
  if( A0.size( ) < 3 || B0.size( ) < 3 ) return A0;

  // 1. Normalize both contours to CCW so correspondence is consistent.
  TContour A = A0, B = B0;
  ensure_ccw( A );
  ensure_ccw( B );

  // 2. Resample both to a common N (capped at 400 for the O(N^2) alignment).
  std::size_t N = std::max( A.size( ), B.size( ) );
  if( N > 400 ) N = 400;
  if( N < 3 )   N = 3;
  const TContour Ar = curvature_adaptive ? resample_adaptive( A, N, lambda )
                                         : resample_uniform( A, N );
  const TContour Br = curvature_adaptive ? resample_adaptive( B, N, lambda )
                                         : resample_uniform( B, N );

  // 3. Cyclic alignment, then 4. linear blend of aligned equal-length contours.
  std::size_t sigma = 0; bool rev = false;
  cyclic_align( Ar, Br, sigma, rev );
  const TContour Ba = ci_detail::apply_shift( Br, sigma, rev );

  TContour C;
  C.reserve( N );
  for( std::size_t i = 0; i < N; ++i )
    C.emplace_back( ( 1.0 - t ) * Ar[ i ].x( ) + t * Ba[ i ].x( ),
                    ( 1.0 - t ) * Ar[ i ].y( ) + t * Ba[ i ].y( ) );
  return C;
}

// =============================================================================
// Pairwise method switch (Linear | Sdf)
// =============================================================================
inline pujCGAL::Final::ContourInterpolator::TContour
pujCGAL::Final::ContourInterpolator::
interpolate( const TContour& A, const TContour& B, double t,
             InterpKind kind, bool align )
{
  switch( kind )
  {
    case InterpKind::Sdf:
      return DistanceField::interpolate( A, B, t, align );
    case InterpKind::Linear:
      return interpolate( A, B, t );
    default:
      std::cerr << "[warn] Polynomial/Spline are not pairwise methods; "
                   "using Linear.\n";
      return interpolate( A, B, t );
  }
}

// =============================================================================
// M2 - series interpolation along z (polynomial / natural cubic spline)
// =============================================================================
inline std::vector< pujCGAL::Final::ContourInterpolator::TContour >
pujCGAL::Final::ContourInterpolator::
interpolate_series( const std::vector< TContour >& slices,
                    const std::vector< double >&   zs,
                    const std::vector< double >&   query_z,
                    InterpKind kind )
{
  std::vector< TContour > out;
  const std::size_t m = slices.size( );
  if( m < 2 || zs.size( ) != m || query_z.empty( ) ) return out;

  // 1. Common N, resample, and chain the cyclic correspondence so that vertex
  //    index i traces a single trajectory gamma_i(z) across all slices.
  std::size_t N = 0;
  for( const auto& s : slices ) N = std::max( N, s.size( ) );
  if( N > 400 ) N = 400;
  if( N < 3 )   N = 3;

  std::vector< TContour > al( m );
  {
    TContour c0 = slices[ 0 ];
    ensure_ccw( c0 );
    al[ 0 ] = resample_uniform( c0, N );
  }
  for( std::size_t k = 1; k < m; ++k )
  {
    TContour ck = slices[ k ];
    ensure_ccw( ck );
    TContour rk = resample_uniform( ck, N );
    std::size_t sigma = 0; bool rev = false;
    cyclic_align( al[ k - 1 ], rk, sigma, rev );    // align to previous frame
    al[ k ] = ci_detail::apply_shift( rk, sigma, rev );
  }

  // 2-3. Per trajectory and coordinate, fit the requested 1D model in z and
  //      evaluate at every query height.
  out.assign( query_z.size( ), TContour( ) );
  for( auto& q : out ) q.reserve( N );

  const int deg = std::min( 3, (int) m - 1 );        // default cubic, clamped
  std::vector< double > xi( m ), yi( m );
  for( std::size_t i = 0; i < N; ++i )
  {
    for( std::size_t k = 0; k < m; ++k )
    { xi[ k ] = al[ k ][ i ].x( ); yi[ k ] = al[ k ][ i ].y( ); }

    if( kind == InterpKind::Spline )
    {
      const auto Mx = ci_detail::spline_M( zs, xi );
      const auto My = ci_detail::spline_M( zs, yi );
      for( std::size_t qz = 0; qz < query_z.size( ); ++qz )
        out[ qz ].emplace_back( ci_detail::spline_eval( zs, xi, Mx, query_z[ qz ] ),
                                ci_detail::spline_eval( zs, yi, My, query_z[ qz ] ) );
    }
    else if( kind == InterpKind::Polynomial )
    {
      const auto cx = ci_detail::polyfit( zs, xi, deg );
      const auto cy = ci_detail::polyfit( zs, yi, deg );
      for( std::size_t qz = 0; qz < query_z.size( ); ++qz )
        out[ qz ].emplace_back( ci_detail::polyeval( cx, query_z[ qz ] ),
                                ci_detail::polyeval( cy, query_z[ qz ] ) );
    }
    else // Linear (piecewise) along z; Sdf is not a series method.
    {
      for( std::size_t qz = 0; qz < query_z.size( ); ++qz )
        out[ qz ].emplace_back( ci_detail::lin_eval( zs, xi, query_z[ qz ] ),
                                ci_detail::lin_eval( zs, yi, query_z[ qz ] ) );
    }
  }
  return out;
}

#endif // __ContourInterpolator__hxx__

// eof - ContourInterpolator.hxx
