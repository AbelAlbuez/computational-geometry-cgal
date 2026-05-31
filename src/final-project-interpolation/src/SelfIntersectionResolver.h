#ifndef __SelfIntersectionResolver__h__
#define __SelfIntersectionResolver__h__

#include "ContourInterpolator.h"

namespace pujCGAL
{
  namespace Final
  {
    /**
     * Detects self-intersections of a closed polygonal contour using the
     * Bentley-Ottmann sweep from pujCGAL, splits the contour at the crossing
     * points, and returns the resulting sub-loop with the largest area.
     */
    class SelfIntersectionResolver
    {
    public:
      using TKernel  = ContourInterpolator::TKernel;
      using TPoint   = ContourInterpolator::TPoint;
      using TContour = ContourInterpolator::TContour;
      using TSegment = TKernel::Segment_2;

    public:
      static TContour resolve( const TContour& contour );

      /// Convenience: tells whether `contour` has any self-intersection.
      static bool has_self_intersections( const TContour& contour );
    };
  } // end namespace
} // end namespace

#include "SelfIntersectionResolver.hxx"

#endif // __SelfIntersectionResolver__h__

// eof - SelfIntersectionResolver.h
