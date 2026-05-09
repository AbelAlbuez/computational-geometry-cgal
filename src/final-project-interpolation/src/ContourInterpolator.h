#ifndef __ContourInterpolator__h__
#define __ContourInterpolator__h__

#include <string>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace pujCGAL
{
  namespace Final
  {
    /**
     * Geometric interpolator between two 2D tumor contours extracted from
     * consecutive axial MRI slices. The detailed implementation lives in the
     * .hxx file; this header only fixes the public surface and the kernel
     * choice required by the project statement.
     */
    class ContourInterpolator
    {
    public:
      using TKernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
      using TPoint   = TKernel::Point_2;
      using TContour = std::vector< TPoint >;

    public:
      ContourInterpolator( ) = default;

      // -- I/O ------------------------------------------------------------
      static TContour read_obj( const std::string& fname );
      static void     write_obj( const std::string& fname,
                                 const TContour& contour );

      // -- Interpolation (to be implemented) ------------------------------
      // Return a contour halfway between A and B (t in [0, 1]).
      static TContour interpolate( const TContour& A,
                                   const TContour& B,
                                   double t );
    };
  } // end namespace
} // end namespace

#include "ContourInterpolator.hxx"

#endif // __ContourInterpolator__h__

// eof - ContourInterpolator.h
