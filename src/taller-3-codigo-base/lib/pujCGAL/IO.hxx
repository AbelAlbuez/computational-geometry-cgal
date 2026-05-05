#ifndef __pujCGAL__IO__hxx__
#define __pujCGAL__IO__hxx__

#include <fstream>

// -------------------------------------------------------------------------
template< class _TDelaunay >
bool pujCGAL::IO::save( const _TDelaunay& T, const std::string& fname )
{
  std::ofstream out( fname );
  if( out )
  {
    std::size_t i = 0;
    for(
      auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v
      )
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
    out.close( );
    return( true );
  }
  else
    return( false );
}

#endif // __pujCGAL__IO__hxx__

// eof - IO.hxx
