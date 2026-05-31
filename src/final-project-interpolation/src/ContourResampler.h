#ifndef __ContourResampler__h__
#define __ContourResampler__h__

#include "ContourInterpolator.h"

namespace pujCGAL
{
  namespace Final
  {
    /**
     * Uniformly redistributes vertices along a closed contour by arc length.
     * Produces a new contour with exactly `n` vertices that lies on the same
     * polyline as the input.
     */
    class ContourResampler
    {
    public:
      using TKernel  = ContourInterpolator::TKernel;
      using TPoint   = ContourInterpolator::TPoint;
      using TContour = ContourInterpolator::TContour;

    public:
      static TContour resample( const TContour& contour, int n );
    };
  } // end namespace
} // end namespace

#include "ContourResampler.hxx"

#endif // __ContourResampler__h__

// eof - ContourResampler.h
