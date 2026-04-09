// =========================================================================
// @author Santiago Gil Gallego (santiago_gil@javeriana.edu.co)
// @author Abel albueez (aa-albuezs@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__DualGraph__hxx__
#define __pujCGAL__DualGraph__hxx__

#include <algorithm>
#include <iterator>
#include <map>
#include <vector>

#include <pujCGAL/Triangulation.h>

// -------------------------------------------------------------------------
template< class _TKernel >
pujCGAL::DualGraph< _TKernel >::
DualGraph( )
{
  this->clear( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
pujCGAL::DualGraph< _TKernel >::
~DualGraph( )
{
  this->clear( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::DualGraph< _TKernel >::
clear( )
{
  this->m_Barycenters.clear( );
  this->m_PointInfinity = TPoint( 0, 0 );
  this->m_IsBoundary.clear( );
  this->m_InternalEdges.clear( );
  this->m_ExternalEdges.clear( );
  this->m_AdjacencyMatrix.clear( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
typename pujCGAL::DualGraph< _TKernel >::TIndex
pujCGAL::DualGraph< _TKernel >::
num_nodes( ) const
{
  return( this->m_Barycenters.size( ) + 1 );
}

// -------------------------------------------------------------------------
template< class _TKernel >
typename pujCGAL::DualGraph< _TKernel >::TIndex
pujCGAL::DualGraph< _TKernel >::
num_internal_edges( ) const
{
  return( this->m_InternalEdges.size( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
typename pujCGAL::DualGraph< _TKernel >::TIndex
pujCGAL::DualGraph< _TKernel >::
num_external_edges( ) const
{
  return( this->m_ExternalEdges.size( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
bool pujCGAL::DualGraph< _TKernel >::
is_boundary( const TIndex& i ) const
{
  return( this->m_IsBoundary[ i ] );
}

// -------------------------------------------------------------------------
template< class _TKernel >
const typename pujCGAL::DualGraph< _TKernel >::TPoint&
pujCGAL::DualGraph< _TKernel >::
point_infinity( ) const
{
  return( this->m_PointInfinity );
}

// -------------------------------------------------------------------------
template< class _TKernel >
const typename pujCGAL::DualGraph< _TKernel >::TPoint&
pujCGAL::DualGraph< _TKernel >::
barycenter( const TIndex& i ) const
{
  return( this->m_Barycenters[ i ] );
}

// -------------------------------------------------------------------------
template< class _TKernel >
bool pujCGAL::DualGraph< _TKernel >::
adjacent( const TIndex& i, const TIndex& j ) const
{
  return( this->m_AdjacencyMatrix[ i ][ j ] );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::DualGraph< _TKernel >::
barycenters_begin( ) const
{
  return( this->m_Barycenters.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::DualGraph< _TKernel >::
barycenters_end( ) const
{
  return( this->m_Barycenters.end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::DualGraph< _TKernel >::
internal_edges_begin( ) const
{
  return( this->m_InternalEdges.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::DualGraph< _TKernel >::
internal_edges_end( ) const
{
  return( this->m_InternalEdges.end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::DualGraph< _TKernel >::
external_edges_begin( ) const
{
  return( this->m_ExternalEdges.begin( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
auto pujCGAL::DualGraph< _TKernel >::
external_edges_end( ) const
{
  return( this->m_ExternalEdges.end( ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void pujCGAL::
build_dual_graph(
  pujCGAL::DualGraph< _TKernel >& dual,
  const pujCGAL::Triangulation< _TKernel >& mesh,
  const std::vector< typename _TKernel::Point_2 >& barycenters
  )
{
  using TIndex   = typename pujCGAL::DualGraph< _TKernel >::TIndex;
  using TPoint   = typename _TKernel::Point_2;
  using TReal    = typename _TKernel::RT;
  using TEdge    = typename pujCGAL::DualGraph< _TKernel >::TEdge;
  using TEdgeKey = std::pair< TIndex, TIndex >;
  using TEdgeMap = std::map< TEdgeKey, std::vector< TIndex > >;

  // -- paso 1: instanciamos el objeto
  dual.clear( );

  // -- paso 2: Copiamos los barycenters
  dual.m_Barycenters = barycenters;

  // -- Contamos los triángulos que vienen de la TOPOLOGÍA (NO GEOMETRÍA)
  const TIndex n_tri = static_cast< TIndex >(
    std::distance( mesh.topology_begin( ), mesh.topology_end( ) )
    );
  const TIndex p_inf_idx = n_tri; // P∞ ocula el índice n_tri en la matriz

  // -- paso 3: Inicializamos el m_IsBoundary
  dual.m_IsBoundary.assign( n_tri, false );

  // -- paso 4: construimos la función de mapa (toma: sorted vertex pair -> devuelve: triangle indices)
  TEdgeMap edge_map;
  TIndex tri_idx = 0;
  for(
    auto tIt = mesh.topology_begin( );
    tIt != mesh.topology_end( );
    ++tIt, ++tri_idx
    )
  {
    const TIndex a = ( *tIt )[ 0 ];
    const TIndex b = ( *tIt )[ 1 ];
    const TIndex c = ( *tIt )[ 2 ];

    edge_map[ std::make_pair( std::min( a, b ), std::max( a, b ) ) ].push_back( tri_idx );
    edge_map[ std::make_pair( std::min( b, c ), std::max( b, c ) ) ].push_back( tri_idx );
    edge_map[ std::make_pair( std::min( a, c ), std::max( a, c ) ) ].push_back( tri_idx );
  } // for

  // -- paso 5: Calcula la matriz de adjacencia, clasifica bordes, llena lista de bordes
  dual.m_AdjacencyMatrix.assign(
    n_tri + 1, std::vector< bool >( n_tri + 1, false )
    );

  for( const auto& entry : edge_map )
  {
    const auto& tris = entry.second;
    if( tris.size( ) == 2 )
    {
      TIndex i = tris[ 0 ];
      TIndex j = tris[ 1 ];
      dual.m_AdjacencyMatrix[ i ][ j ] = true;
      dual.m_AdjacencyMatrix[ j ][ i ] = true;
      dual.m_InternalEdges.push_back( { i, j } );
    }
    else if( tris.size( ) == 1 )
    {
      TIndex i = tris[ 0 ];
      dual.m_IsBoundary[ i ] = true;
      dual.m_AdjacencyMatrix[ i ][ p_inf_idx ] = true;
      dual.m_AdjacencyMatrix[ p_inf_idx ][ i ] = true;
      dual.m_ExternalEdges.push_back( { i, p_inf_idx } );
    } // if
  } // for

  // -- paso 6: Computa P∞ de la caja envolvente
  auto gIt = mesh.geometry_begin( );
  TReal xmin = ( *gIt )[ 0 ], xmax = ( *gIt )[ 0 ];
  TReal ymin = ( *gIt )[ 1 ], ymax = ( *gIt )[ 1 ];
  for( ++gIt; gIt != mesh.geometry_end( ); ++gIt )
  {
    TReal x = ( *gIt )[ 0 ];
    TReal y = ( *gIt )[ 1 ];
    if( x < xmin ) xmin = x;
    if( x > xmax ) xmax = x;
    if( y < ymin ) ymin = y;
    if( y > ymax ) ymax = y;
  } //  for

  TReal cx = ( xmin + xmax ) / TReal( 2 );
  TReal cy = ( ymin + ymax ) / TReal( 2 );
  TReal dx = xmax - xmin;
  TReal dy = ymax - ymin;
  dual.m_PointInfinity = TPoint( cx + TReal( 2 ) * dx, cy - TReal( 2 ) * dy );
}

#endif // __pujCGAL__DualGraph__hxx__

// eof - DualGraph.hxx
