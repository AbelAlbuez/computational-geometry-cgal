// =========================================================================
// @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__Triangulation__hxx__
#define __pujCGAL__Triangulation__hxx__

// -------------------------------------------------------------------------
template< class _TKernel >
pujCGAL::Triangulation< _TKernel >::
Triangulation( )
{
  this->clear( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
pujCGAL::Triangulation< _TKernel >::
~Triangulation( )
{
  this->clear( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::Triangulation< _TKernel >::
clear( )
{
  this->m_Points.clear( );
  this->m_Topology.clear( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::Triangulation< _TKernel >::
add_point( const TPoint& p )
{
  this->m_Points.push_back( p );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::Triangulation< _TKernel >::
add_triangle( const TIndex& a, const TIndex& b, const TIndex& c )
{
  this->m_Topology.push_back( { a, b, c } );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Triangulation< _TKernel >::
geometry_begin( )
{
  return( this->m_Points.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Triangulation< _TKernel >::
geometry_begin( ) const
{
  return( this->m_Points.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Triangulation< _TKernel >::
geometry_end( )
{
  return( this->m_Points.end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Triangulation< _TKernel >::
geometry_end( ) const
{
  return( this->m_Points.end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Triangulation< _TKernel >::
topology_begin( )
{
  return( this->m_Topology.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Triangulation< _TKernel >::
topology_begin( ) const
{
  return( this->m_Topology.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Triangulation< _TKernel >::
topology_end( )
{
  return( this->m_Topology.end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Triangulation< _TKernel >::
topology_end( ) const
{
  return( this->m_Topology.end( ) );
}

#endif // __pujCGAL__Triangulation__hxx__

// eof - Triangulation.hxx
