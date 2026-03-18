#ifndef __pujCGAL__IO__h__
#define __pujCGAL__IO__h__

namespace pujCGAL
{
  namespace IO
  {
    /**
     */
    template< class TData >
    bool read( const std::string& fname, TData data );
 
    /**
     */
    template< class TData >
    bool save( const std::string& fname, TData data );
  } // end namespace
} // end namespace

#include <pujCGAL/IO.hxx>

#endif // __pujCGAL__IO__h__

// eof - IO.h
