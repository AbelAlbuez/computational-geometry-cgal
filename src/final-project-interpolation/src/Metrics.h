#ifndef __Metrics__h__
#define __Metrics__h__

#include <vector>

#include "Reconstructor.h"

namespace pujCGAL
{
  namespace Final
  {
    /**
     * Stage 5: descriptors of a reconstructed surface and a contour-stack
     * cross-check.
     *
     *   volume_area  : PMP::volume (mm^3, after orient_to_bound_a_volume) and
     *                  PMP::area (mm^2).
     *   stack_volume : trapezoidal integral of the shoelace areas * dz -- an
     *                  independent volume estimate from the raw slices.
     *   symmetry     : reflective (bilateral) symmetry score. The inertia-tensor
     *                  principal axes give three candidate mirror planes; each
     *                  reflection is scored by a volumetric Dice (grid +
     *                  Side_of_triangle_mesh) and a symmetric Hausdorff distance.
     *   mean curvature: discrete cotangent-Laplacian H (Meyer et al. 2003) --
     *                  hand-rolled because CGAL 5.6 lacks the 6.x one-call
     *                  interpolated_corrected_curvatures; summary statistics.
     */
    class Metrics
    {
    public:
      using Mesh     = Reconstructor::Mesh;          // Surface_mesh<Point_3>, EPIC
      using Slice    = Reconstructor::Slice;
      using TContour = ContourInterpolator::TContour;

      struct VolArea  { double volume_mm3 = 0, area_mm2 = 0; };
      struct Symmetry { int axis = -1; double dice = 0, hausdorff = 0; };
      struct CurvStats{ double mean = 0, stddev = 0, min = 0, max = 0,
                               median = 0, p90 = 0; std::size_t n = 0; };

    public:
      static VolArea  volume_area( Mesh& m );        // reorients m in place
      static double   stack_volume( const std::vector< Slice >& stack );

      static Symmetry  symmetry_score( const Mesh& m, int grid_res = 40 );
      static CurvStats mean_curvature_stats( const Mesh& m );
    };
  } // end namespace Final
} // end namespace pujCGAL

#include "Metrics.hxx"

#endif // __Metrics__h__

// eof - Metrics.h
