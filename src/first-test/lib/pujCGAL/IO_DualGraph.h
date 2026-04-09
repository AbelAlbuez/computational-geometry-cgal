// =========================================================================
// @author Santiago Gil Gallego (santiago_gil@javeriana.edu.co)
// @author Abel albueez (aa-albuezs@javeriana.edu.co)
// =========================================================================
#ifndef __pujCGAL__IO_DualGraph__h__
#define __pujCGAL__IO_DualGraph__h__

#include <string>

#include <pujCGAL/DualGraph.h>
#include <pujCGAL/IO.h>

namespace pujCGAL
{
  namespace IO
  {
    /**
     */
    template< class TKernel >
    bool save( const std::string& fname, const DualGraph< TKernel >& dual );
  } // end namespace
} // end namespace

#include <pujCGAL/IO_DualGraph.hxx>

#endif // __pujCGAL__IO_DualGraph__h__

// eof - IO_DualGraph.h
