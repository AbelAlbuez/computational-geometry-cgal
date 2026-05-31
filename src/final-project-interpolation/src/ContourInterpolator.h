#ifndef __ContourInterpolator__h__
#define __ContourInterpolator__h__

#include <string>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "DistanceField.h"   // standalone SDF module (same kernel); no back-dependency

namespace pujCGAL
{
  namespace Final
  {
    /**
     * Geometric interpolator between 2D tumor contours extracted from
     * consecutive axial MRI slices. Implementation lives in the .hxx file.
     *
     *   M1 Linear   : pairwise vertex-correspondence LERP.
     *   M2 Series   : polynomial / natural-cubic-spline along z (>= 2 slices).
     *   M3 Shape    : signed-distance-field blend (delegated to DistanceField).
     */
    class ContourInterpolator
    {
    public:
      using TKernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
      using TPoint   = TKernel::Point_2;
      using TContour = std::vector< TPoint >;

      // Which interpolation model to use. Linear/Sdf are pairwise (A,B);
      // Polynomial/Spline act on a series of slices along z.
      enum class InterpKind { Linear, Polynomial, Spline, Sdf };

    public:
      ContourInterpolator( ) = default;

      // -- I/O ------------------------------------------------------------
      static TContour read_obj( const std::string& fname );
      static void     write_obj( const std::string& fname,
                                 const TContour& contour );

      // -- Geometry helpers (public for reuse / testing) ------------------
      static double  signed_area( const TContour& c );   // shoelace; CCW > 0
      static TPoint  centroid( const TContour& c );       // vertex mean
      static double  perimeter( const TContour& c );      // closed
      static void    ensure_ccw( TContour& c );           // reverse if CW
      static TContour resample_uniform( const TContour& c, std::size_t N );
      static TContour resample_adaptive( const TContour& c, std::size_t N,
                                         double lambda );

      // Best cyclic shift (and traversal direction) of B onto A; A,B must
      // already be resampled to the same length. Used by M1 and the M2 chain.
      static void    cyclic_align( const TContour& A, const TContour& B,
                                   std::size_t& sigma, bool& reversed );

      // -- M1: linear (vertex-correspondence LERP) ------------------------
      static TContour interpolate( const TContour& A,
                                   const TContour& B,
                                   double t,
                                   bool   curvature_adaptive = false,
                                   double lambda             = 2.0 );

      // -- Pairwise method switch (Linear | Sdf) --------------------------
      // Polynomial/Spline are not valid pairwise and fall back to Linear.
      static TContour interpolate( const TContour& A,
                                   const TContour& B,
                                   double t,
                                   InterpKind kind,
                                   bool align = true );

      // -- M2: series interpolation along z -------------------------------
      // Resamples every slice to a common N, propagates the cyclic shift
      // transitively so vertex i traces a trajectory gamma_i(z), fits each
      // coordinate with the requested 1D model, and evaluates at query_z.
      static std::vector< TContour > interpolate_series(
          const std::vector< TContour >& slices,
          const std::vector< double >&   zs,
          const std::vector< double >&   query_z,
          InterpKind kind );
    };
  } // end namespace Final
} // end namespace pujCGAL

#include "ContourInterpolator.hxx"

#endif // __ContourInterpolator__h__

// eof - ContourInterpolator.h
