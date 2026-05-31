#ifndef __LinearInterpolator__h__
#define __LinearInterpolator__h__

#include "ContourInterpolator.h"

namespace pujCGAL
{
  namespace Final
  {
    /**
     * Vertex-wise linear interpolation between two contours with the same
     * number of vertices.
     */
    class LinearInterpolator
    {
    public:
      using TKernel  = ContourInterpolator::TKernel;
      using TPoint   = ContourInterpolator::TPoint;
      using TContour = ContourInterpolator::TContour;

    public:
      static TContour interpolate( const TContour& A,
                                   const TContour& B,
                                   double t );
    };
  } // end namespace
} // end namespace

#include "LinearInterpolator.hxx"

#endif // __LinearInterpolator__h__

// eof - LinearInterpolator.h
