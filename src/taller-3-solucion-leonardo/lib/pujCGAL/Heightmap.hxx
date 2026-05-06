#ifndef __pujCGAL__Heightmap__hxx__
#define __pujCGAL__Heightmap__hxx__

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <png.h>

// -------------------------------------------------------------------------
template< class _TReal >
pujCGAL::Heightmap< _TReal >::
Heightmap( )
{
}

// -------------------------------------------------------------------------
template< class _TReal >
pujCGAL::Heightmap< _TReal >::
~Heightmap( )
{
  if( this->m_Buffer != nullptr )
    std::free( this->m_Buffer );
}

// -------------------------------------------------------------------------
template< class _TReal >
const std::size_t& pujCGAL::Heightmap< _TReal >::
width( ) const
{
  return( this->m_Width );
}

// -------------------------------------------------------------------------
template< class _TReal >
const std::size_t& pujCGAL::Heightmap< _TReal >::
height( ) const
{
  return( this->m_Height );
}

// -------------------------------------------------------------------------
template< class _TReal >
void pujCGAL::Heightmap< _TReal >::
set_origin( const TReal& x, const TReal& y )
{
  this->m_Origin[ 0 ] = x;
  this->m_Origin[ 1 ] = y;
}

// -------------------------------------------------------------------------
template< class _TReal >
void pujCGAL::Heightmap< _TReal >::
set_dimensions( const TReal& w, const TReal& h )
{
  this->m_Dimensions[ 0 ] = w;
  this->m_Dimensions[ 1 ] = h;
}

// -------------------------------------------------------------------------
template< class _TReal >
const typename pujCGAL::Heightmap< _TReal >::
TReal& pujCGAL::Heightmap< _TReal >::
operator()( const std::size_t& w, const std::size_t& h ) const
{
  return( *( this->m_Buffer + ( ( h * this->m_Width ) + w ) ) );
}


// -------------------------------------------------------------------------
template< class _TReal >
typename pujCGAL::Heightmap< _TReal >::
TPoint pujCGAL::Heightmap< _TReal >::
point( const std::size_t& w, const std::size_t& h ) const
{
  return(
    std::make_pair(
      ( this->m_Dimensions[ 0 ] * TReal( w ) / TReal( this->m_Width ) )
      +
      this->m_Origin[ 0 ],
      ( this->m_Dimensions[ 1 ] * TReal( h ) / TReal( this->m_Height ) )
      +
      this->m_Origin[ 1 ]
      )
    );
}

// -------------------------------------------------------------------------
template< class _TReal >
typename pujCGAL::Heightmap< _TReal >::
TIndex pujCGAL::Heightmap< _TReal >::
index( const TReal& x, const TReal& y )
{
  return(
    std::make_pair(
      std::size_t( TReal( this->m_Width ) * ( ( x - this->m_Origin[ 0 ] ) / this->m_Dimensions[ 0 ] ) ),
      std::size_t( TReal( this->m_Height ) * ( ( y - this->m_Origin[ 1 ] ) / this->m_Dimensions[ 1 ] ) )
      )
    );
}

// -------------------------------------------------------------------------
template< class _TReal >
void pujCGAL::Heightmap< _TReal >::
read_from_png( const std::string& fname )
{
  std::FILE* fp = std::fopen( fname.c_str( ), "rb" );
  if( !fp )
    return;

  // 1. Initialize PNG structures
  png_structp png = png_create_read_struct(
    PNG_LIBPNG_VER_STRING, NULL, NULL, NULL
    );
  png_infop info = png_create_info_struct( png );

  // 2. Set error handling (required by libpng)
  if( setjmp( png_jmpbuf( png ) ) )
  {
    png_destroy_read_struct( &png, &info, NULL );
    std::fclose( fp );
    return;
  } // end if

  // 3. Initialize I/O
  png_init_io( png, fp );
  png_read_info( png, info );

  // 4. Get image details
  std::uint8_t bit_depth = png_get_bit_depth( png, info );
  if( bit_depth != 8 && bit_depth != 16 )
    return;
  bit_depth >>= 3;

  this->m_Width = png_get_image_width( png, info );
  this->m_Height = png_get_image_height( png, info );
  std::uint8_t color_type = png_get_color_type( png, info );

  std::size_t components;
  switch( color_type )
  {
  case PNG_COLOR_TYPE_GRAY:       components = 1; break;
  case PNG_COLOR_TYPE_GRAY_ALPHA: components = 2; break;
  case PNG_COLOR_TYPE_RGB:        components = 3; break;
  case PNG_COLOR_TYPE_RGB_ALPHA:  components = 4; break;
  default:                        components = 0; break;
  } // end switch

  // 5. Update info (if necessary)
  png_read_update_info( png, info );

  // 6. Allocate memory for image
  png_bytep* rows
    =
    reinterpret_cast< png_bytep* >(
      std::calloc( this->m_Height, sizeof( png_bytep ) )
      );
  for( std::size_t y = 0; y < this->m_Height; ++y )
    rows[ y ]
      =
      reinterpret_cast< png_byte* >(
        std::malloc( png_get_rowbytes( png, info ) )
        );

  // 7. Read the actual image data
  png_read_image( png, rows );

  // 8. Prepare TReal memory buffer
  std::size_t buffer_size = this->m_Width * this->m_Height;
  this->m_Buffer
    =
    reinterpret_cast< TReal* >(
      std::calloc( buffer_size, sizeof( TReal ) )
      );

  // 9. Copy into a buffer and cleanup
  TReal maxV = std::numeric_limits< TReal >::lowest( );
  TReal* image = this->m_Buffer;
  for( std::size_t y = 0; y < this->m_Height; ++y )
  {
    png_byte* r = rows[ y ];
    for( std::size_t x = 0; x < this->m_Width; ++x )
    {
      *image = 0;
      for( std::size_t c = 0; c < components; ++c )
      {
        TReal v = 0;
        if( bit_depth == 1 )
          v = TReal( *( reinterpret_cast< std::uint8_t* >( r ) ) );
        else if( bit_depth == 2 )
          v = TReal( *( reinterpret_cast< std::uint16_t* >( r ) ) );
        *image += v * v;
        r += bit_depth * sizeof( png_byte );
      } // end for
      *image = std::sqrt( *image );
      maxV = ( maxV < *image )? *image: maxV;
      image++;
    } // end for
    std::free( rows[ y ] );
  } // end for
  std::free( rows );
  std::fclose( fp );
  png_destroy_read_struct( &png, &info, NULL );

  // 10. normalize map
  image = this->m_Buffer;
  for( std::size_t i = 0; i < buffer_size; ++i )
  {
    *image /= maxV;
    image++;
  } // end for
}

#endif // __pujCGAL__Heightmap__hxx__

// eof - Heightmap.hxx
