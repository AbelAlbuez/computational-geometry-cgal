#ifndef __pujCGAL__IO__hxx__
#define __pujCGAL__IO__hxx__

#include <fstream>
#include <sstream>
#include <vector>

// -------------------------------------------------------------------------
template< class _TMesh >
bool pujCGAL::IO::
read( _TMesh& M, const std::string& fname )
{
  using TPoint  = typename _TMesh::Point;
  using TKernel = typename CGAL::Kernel_traits< TPoint >::Kernel;
  using TReal   = typename TKernel::RT;
  using TVertex = typename _TMesh::Vertex_index;

  std::ifstream ifs( fname.c_str( ), std::ios::binary | std::ios::ate );
  if( !ifs )
    return( false );

  std::streamsize ifs_size = ifs.tellg( );
  ifs.seekg( 0, std::ios::beg );
  std::vector< char > buffer( ifs_size );
  if( !( ifs.read( buffer.data( ), ifs_size ) ) )
    return( false );

  std::istringstream iss( std::string( buffer.begin( ), buffer.end( ) ) );
  ifs.close( );

  std::vector< TVertex > P;
  std::vector< std::tuple< std::size_t, std::size_t, std::size_t > > F;

  std::string line;
  while( std::getline( iss, line ) )
  {
    if( line[ 0 ] == 'v' )
    {
      TReal x, y, z;
      std::istringstream( line.substr( 1 ) ) >> x >> y >> z;
      P.push_back( M.add_vertex( TPoint( x, y, z ) ) );
    }
    else if( line[ 0 ] == 'f' )
    {
      std::size_t a, b, c;
      std::istringstream( line.substr( 1 ) ) >> a >> b >> c;
      F.push_back( std::make_tuple( a - 1, b - 1, c - 1 ) );
    } // end if
  } // end while

  for( const auto& f: F )
  {
    TVertex a = P[ std::get< 0 >( f ) ];
    TVertex b = P[ std::get< 1 >( f ) ];
    TVertex c = P[ std::get< 2 >( f ) ];
    M.add_face( a, b, c );
  } // end for
  return( true );
}

// -------------------------------------------------------------------------
template< class _TDelaunay >
bool pujCGAL::IO::save( const _TDelaunay& T, const std::string& fname )
{
  std::stringstream out;
  std::size_t i = 0;
  for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
  {
    auto h = v->info( );
    v->info( ) = ++i;
    out
      << "v "
      << v->point( ).x( ) << " " << v->point( ).y( ) << " " << h
      << std::endl;
  } // end for

  for( auto f = T.finite_faces_begin( ); f != T.finite_faces_end( ); ++f )
    out
      << "f "
      << std::size_t( f->vertex( 0 )->info( ) ) << " "
      << std::size_t( f->vertex( 1 )->info( ) ) << " "
      << std::size_t( f->vertex( 2 )->info( ) ) << std::endl;

  // Real write
  std::ofstream rout( fname );
  if( rout )
  {
    rout.write( out.str( ).data( ), out.str( ).size( ) );
    rout.close( );
    return( true );
  }
  else
    return( false );
}

#endif // __pujCGAL__IO__hxx__

// eof - IO.hxx
