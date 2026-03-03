#ifndef __pujCGAL__IO__h__
#define __pujCGAL__IO__h__

namespace pujCGAL
{
  namespace IO
  {
    template< class TSegmentsIt >
    bool read( const std::string& fname, TSegmentsIt sIt );

    template< class TPointsIt >
    bool save( const std::string& fname, TPointsIt pB, TPointsIt pE );

  } // end namespace
} // end namespace

#include <pujCGAL/IO.hxx>

#endif // __pujCGAL__IO__h__

// eof - IO.h
