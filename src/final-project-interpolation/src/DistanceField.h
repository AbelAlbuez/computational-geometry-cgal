#ifndef __DistanceField__h__
#define __DistanceField__h__

#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace pujCGAL
{
  namespace Final
  {
    /**
     * Shape-based ("fancy") interpolation via signed distance fields -- M3.
     *
     * The object we manipulate is the signed distance transform of a contour:
     *   Phi_C(q) = +/- min_{c in dC} || q - c ||,   (+ inside, - outside).
     * The map q -> min_c ||q-c|| is exactly the distance to the contour, and
     * the partition of the plane into "nearest segment" cells is the VORONOI
     * DIAGRAM of the contour's edges; the locus where the nearest point is
     * non-unique (the ridge set of Phi) is the contour's MEDIAL AXIS. That is
     * why this interpolation is the most computational-geometry-flavored one.
     *
     * Pipeline:  build Phi_A, Phi_B over a common grid (closest distance from a
     * CGAL::AABB_tree of the segments, sign from Polygon_2::bounded_side) ->
     * blend Phi_t = (1-t)Phi_A + t Phi_B -> extract the zero level set with
     * marching squares. Optional rigid pre-alignment (centroid + principal
     * axis) so features translate/rotate instead of fading in and out.
     */
    class DistanceField
    {
    public:
      using TKernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
      using TPoint   = TKernel::Point_2;
      using TContour = std::vector< TPoint >;

    public:
      // Pairwise SDF interpolation of two single-loop contours. Returns the
      // largest extracted component of the blended zero level set.
      static TContour interpolate( const TContour& A, const TContour& B,
                                   double t, bool align = true,
                                   double spacing = 1.0 );

      // Multi-loop variant (handles topology changes). Returns *all* extracted
      // loops of the blended zero set, largest first. Used by the topology
      // robustness test.
      static std::vector< TContour > interpolate_loops(
          const std::vector< TContour >& A_loops,
          const std::vector< TContour >& B_loops,
          double t, bool align = true, double spacing = 1.0 );

      // -- helpers (public for testing) -----------------------------------
      static double signed_area( const TContour& c );
      static TPoint centroid( const TContour& c );
    };
  } // end namespace Final
} // end namespace pujCGAL

#include "DistanceField.hxx"

#endif // __DistanceField__h__

// eof - DistanceField.h
