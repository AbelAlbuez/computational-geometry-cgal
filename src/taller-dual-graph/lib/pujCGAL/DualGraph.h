// =========================================================================
// @author Santiago Gil Gallego (santiago_gil@javeriana.edu.co)
// @author Abel albueez (aa-albuezs@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__DualGraph__h__
#define __pujCGAL__DualGraph__h__

#include <map>
#include <utility>
#include <vector>

namespace pujCGAL
{
  template< class _TKernel >
  class Triangulation;

  /**
  Aca va el uso del template de triangulacion y dentro va el cálculo
  del grafo dual
   */
  template< class _TKernel >
  class DualGraph
  {
  public:
    using TKernel = _TKernel;
    using TReal   = typename TKernel::RT;
    using TPoint  = typename TKernel::Point_2;
    using TIndex  = std::size_t;
    using TEdge   = std::pair< TIndex, TIndex >;

  public:
    DualGraph( );
    virtual ~DualGraph( );

    void clear( );

    TIndex num_nodes( ) const;
    TIndex num_internal_edges( ) const;
    TIndex num_external_edges( ) const;

    bool          is_boundary( const TIndex& i ) const;
    const TPoint& point_infinity( ) const;
    const TPoint& barycenter( const TIndex& i ) const;
    bool          adjacent( const TIndex& i, const TIndex& j ) const;

    auto barycenters_begin( ) const;
    auto barycenters_end( ) const;

    auto internal_edges_begin( ) const;
    auto internal_edges_end( ) const;

    auto external_edges_begin( ) const;
    auto external_edges_end( ) const;

    template< class _TK >
    friend void build_dual_graph(
      DualGraph< _TK >&,
      const Triangulation< _TK >&,
      const std::vector< typename _TK::Point_2 >&
      );

  protected:
    std::vector< TPoint >              m_Barycenters;
    TPoint                             m_PointInfinity;
    std::vector< bool >                m_IsBoundary;
    std::vector< TEdge >               m_InternalEdges;
    std::vector< TEdge >               m_ExternalEdges;
    std::vector< std::vector< bool > > m_AdjacencyMatrix;
  };

  /**
  acá estamos usando el template de llamar a la función
   */
  template< class _TKernel >
  void build_dual_graph(
    DualGraph< _TKernel >& dual,
    const Triangulation< _TKernel >& mesh,
    const std::vector< typename _TKernel::Point_2 >& barycenters
    );

} // end namespace

#include <pujCGAL/DualGraph.hxx>

#endif // __pujCGAL__DualGraph__h__

// eof - DualGraph.h
