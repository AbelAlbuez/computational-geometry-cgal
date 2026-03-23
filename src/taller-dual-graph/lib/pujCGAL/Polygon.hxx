// =========================================================================
// @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__Polygon__hxx__
#define __pujCGAL__Polygon__hxx__

// -------------------------------------------------------------------------
template< class _TKernel >
pujCGAL::Polygon< _TKernel >::
Polygon( )
{
}

// -------------------------------------------------------------------------
template< class _TKernel >
pujCGAL::Polygon< _TKernel >::
~Polygon( )
{
  if( this->m_Traits != nullptr )
    delete this->m_Traits;
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::Polygon< _TKernel >::
add_point( const TPoint& p )
{
  this->m_Points.push_back( p );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::Polygon< _TKernel >::
build( )
{
  this->m_Traits = new TTraits( CGAL::make_property_map( this->m_Points ) );
  for( std::size_t i = 0; i < this->m_Points.size( ); ++i )
    this->m_Topology.push_back( i );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::Polygon< _TKernel >::
guarantee_CW( )
{
  if( this->area( ) > TReal( 0 ) )
    this->invert( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::Polygon< _TKernel >::
guarantee_CCW( )
{
  if( this->area( ) < TReal( 0 ) )
    this->invert( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::Polygon< _TKernel >::
invert( )
{
  std::reverse(
    this->m_Topology.vertices_begin( ), this->m_Topology.vertices_end( )
    );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Polygon< _TKernel >::
geometry_begin( )
{
  return( this->m_Points.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Polygon< _TKernel >::
geometry_begin( ) const
{
  return( this->m_Points.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Polygon< _TKernel >::
geometry_end( )
{
  return( this->m_Points.end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Polygon< _TKernel >::
geometry_end( ) const
{
  return( this->m_Points.end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Polygon< _TKernel >::
topology_begin( )
{
  return( this->m_Topology.vertices_begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Polygon< _TKernel >::
topology_begin( ) const
{
  return( this->m_Topology.vertices_begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Polygon< _TKernel >::
topology_end( )
{
  return( this->m_Topology.vertices_end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::Polygon< _TKernel >::
topology_end( ) const
{
  return( this->m_Topology.vertices_end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
const typename pujCGAL::Polygon< _TKernel >::
TTraits& pujCGAL::Polygon< _TKernel >::
traits( ) const
{
  return( *( this->m_Traits ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
const typename pujCGAL::Polygon< _TKernel >::
TPoint& pujCGAL::Polygon< _TKernel >::
point( const TIndex& i ) const
{
  return( this->m_Points[ i ] );
}

// -------------------------------------------------------------------------
template< class _TKernel >
typename pujCGAL::Polygon< _TKernel >::
TReal pujCGAL::Polygon< _TKernel >::
area( ) const
{
  TReal a = 0;
  for(
    auto pIt = this->topology_begin( ); pIt != this->topology_end( ); ++pIt
    )
  {
    auto qIt = std::next( pIt );
    if( qIt == this->topology_end( ) )
      qIt = this->topology_begin( );

    auto p = this->m_Points[ *pIt ];
    auto q = this->m_Points[ *qIt ];

    a += p[ 0 ] * q[ 1 ];
    a -= q[ 0 ] * p[ 1 ];
  } // end for
  return( a * TReal( 0.5 ) );
}

#endif // __pujCGAL__Polygon__hxx__

// eof - Polygon.hxx
