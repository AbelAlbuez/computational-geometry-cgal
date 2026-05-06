#ifndef __pujCGAL__IO__h__
#define __pujCGAL__IO__h__

#include <string>

namespace pujCGAL
{
  namespace IO
  {
    /**
     */
    template< class _TMesh >
    bool read( _TMesh& M, const std::string& fname );

    /**
     */
    template< class _TDelaunay >
    bool save( const _TDelaunay& T, const std::string& fname );
  } // end namespace
} // end namespace

#include <pujCGAL/IO.hxx>

#endif // __pujCGAL__IO__h__

// eof - IO.h
