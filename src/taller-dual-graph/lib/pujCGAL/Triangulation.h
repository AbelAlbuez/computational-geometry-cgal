// =========================================================================
// @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__Triangulation__h__
#define __pujCGAL__Triangulation__h__

#include <array>
#include <vector>

namespace pujCGAL
{
  /**
   */
  template< class _TKernel >
  class Triangulation
  {
  public:
    using TKernel = _TKernel;
    using TReal = typename TKernel::RT;
    using TPoint = typename TKernel::Point_2;
    using TIndex = std::size_t;
    using TPoints = std::vector< TPoint >;
    using TTriangle = std::array< TIndex, 3 >;
    using TTopology = std::vector< TTriangle >;

  public:
    Triangulation( );
    virtual ~Triangulation( );

    void clear( );
    void add_point( const TPoint& p );
    void add_triangle( const TIndex& a, const TIndex& b, const TIndex& c );

    auto geometry_begin( );
    auto geometry_begin( ) const;
    auto geometry_end( );
    auto geometry_end( ) const;

    auto topology_begin( );
    auto topology_begin( ) const;
    auto topology_end( );
    auto topology_end( ) const;

  protected:
    TPoints   m_Points;
    TTopology m_Topology;
  };
} // end namespace

#include <pujCGAL/Triangulation.hxx>

#endif // __pujCGAL__Triangulation__h__

// eof - Triangulation.h
