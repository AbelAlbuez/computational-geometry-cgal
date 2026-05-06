#ifndef __pujCGAL__Heightmap__h__
#define __pujCGAL__Heightmap__h__

#include <cstdint>
#include <string>
#include <utility>

namespace pujCGAL
{
  /**
   */
  template< class _TReal = double >
  class Heightmap
  {
  public:
    using TReal = _TReal;
    using TPoint = std::pair< TReal, TReal >;
    using TIndex = std::pair< std::size_t, std::size_t >;

  public:
    Heightmap( );
    virtual ~Heightmap( );

    const std::size_t& width( ) const;
    const std::size_t& height( ) const;

    void set_origin( const TReal& x, const TReal& y );
    void set_dimensions( const TReal& w, const TReal& h );

    const TReal& operator()( const std::size_t& w, const std::size_t& h ) const;
    TPoint point( const std::size_t& w, const std::size_t& h ) const;
    TIndex index( const TReal& x, const TReal& y );

    void read_from_png( const std::string& fname );

  protected:
    TReal*      m_Buffer { nullptr };
    std::size_t m_Width  { 0 };
    std::size_t m_Height { 0 };

    TReal m_Origin[ 2 ]     { 0, 0 };
    TReal m_Dimensions[ 2 ] { 1, 1 };
  };
} // end namespace

#include <pujCGAL/Heightmap.hxx>

#endif // __pujCGAL__Heightmap__h__

// eof - Heightmap.h
