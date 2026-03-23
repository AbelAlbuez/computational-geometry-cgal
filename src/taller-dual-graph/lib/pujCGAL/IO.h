#ifndef __pujCGAL__IO__h__
#define __pujCGAL__IO__h__

#include <string>

namespace pujCGAL
{
  template< class _TKernel >
  class Polygon;

  template< class _TKernel >
  class Triangulation;

  namespace IO
  {
    /**
     */
    template< class TKernel >
    bool read( const std::string& fname, Polygon< TKernel >& polygon );
 
    /**
     */
    template< class TKernel >
    bool save( const std::string& fname, const Triangulation< TKernel >& mesh );
  } // end namespace
} // end namespace

#include <pujCGAL/IO.hxx>

#endif // __pujCGAL__IO__h__

// eof - IO.h
