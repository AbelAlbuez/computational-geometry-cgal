// =========================================================================
// @author Santiago Gil Gallego (santiago_gil@javeriana.edu.co)
// @author Abel albueez (aa-albuezs@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__IO_DualGraph__hxx__
#define __pujCGAL__IO_DualGraph__hxx__

#include <fstream>

// -------------------------------------------------------------------------
template< class TKernel >
bool pujCGAL::IO::
save( const std::string& fname, const pujCGAL::DualGraph< TKernel >& dual )
{
  std::ofstream ofs( fname.c_str( ) );
  if( ofs )
  {
    ofs << "# Dual graph: "
        << dual.num_nodes( ) << " nodes, "
        << dual.num_internal_edges( ) << " internal edges, "
        << dual.num_external_edges( ) << " external edges" << std::endl;

    for(
      auto bIt = dual.barycenters_begin( ); bIt != dual.barycenters_end( ); ++bIt
      )
      ofs << "v "
          << ( *bIt )[ 0 ] << " "
          << ( *bIt )[ 1 ] << " 0" << std::endl;

    ofs << "v "
        << dual.point_infinity( )[ 0 ] << " "
        << dual.point_infinity( )[ 1 ] << " 0" << std::endl;

    ofs << std::endl;

    for(
      auto eIt = dual.internal_edges_begin( ); eIt != dual.internal_edges_end( ); ++eIt
      )
      ofs << "l "
          << ( ( *eIt ).first  + 1 ) << " "
          << ( ( *eIt ).second + 1 ) << std::endl;

    for(
      auto eIt = dual.external_edges_begin( ); eIt != dual.external_edges_end( ); ++eIt
      )
      ofs << "l "
          << ( ( *eIt ).first  + 1 ) << " "
          << ( ( *eIt ).second + 1 ) << std::endl;

    ofs << std::endl;

    ofs.close( );
    return( true );
  }
  else
    return( false );
}

#endif // __pujCGAL__IO_DualGraph__hxx__

// eof - IO_DualGraph.hxx
