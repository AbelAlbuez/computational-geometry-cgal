// =========================================================================
// @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__Polygon__h__
#define __pujCGAL__Polygon__h__

#include <vector>

#include <CGAL/property_map.h>
#include <CGAL/Partition_traits_2.h>

namespace pujCGAL
{
  /**
   */
  template< class _TKernel >
  class Polygon
  {
  public:
    using TKernel = _TKernel;
    using TReal = typename TKernel::RT;
    using TPoint = typename TKernel::Point_2;
    using TMap = CGAL::Pointer_property_map< TPoint >;
    using TTraits = CGAL::Partition_traits_2< TKernel, typename TMap::type >;
    using TIndex = typename TTraits::Point_2;
    using TTopology = typename TTraits::Polygon_2;

  public:
    Polygon( );
    virtual ~Polygon( );

    void add_point( const TPoint& p );
    void build( );

    void guarantee_CW( );
    void guarantee_CCW( );
    void invert( );

    auto geometry_begin( );
    auto geometry_begin( ) const;
    auto geometry_end( );
    auto geometry_end( ) const;

    auto topology_begin( );
    auto topology_begin( ) const;
    auto topology_end( );
    auto topology_end( ) const;

    const TTraits& traits( ) const;
    const TPoint& point( const TIndex& i ) const;
    TReal area( ) const;

  protected:
    std::vector< TPoint > m_Points;
    TTopology m_Topology;
    TTraits* m_Traits { nullptr };
  };

  /**
   */
  template< class _TTriangulation, class _TPolygon >
  void triangulate( _TTriangulation& mesh, const _TPolygon& polygon );

  /**
   */
  template< class _TTriangulation, class _TPolygon >
  void triangulate(
    _TTriangulation& T,
    const _TPolygon& P,
    const typename _TPolygon::TTopology& M
    );

} // end namespace

#include <pujCGAL/Polygon.hxx>
#include <pujCGAL/Polygon_Functions.hxx>

#endif // __pujCGAL__Polygon__h__

// eof - Polygon.h
