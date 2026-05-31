#ifndef __Reconstructor__h__
#define __Reconstructor__h__

#include <string>
#include <vector>

#include <CGAL/Surface_mesh.h>

#include "ContourInterpolator.h"

namespace pujCGAL
{
  namespace Final
  {
    /**
     * Stage 3-4: turn a densified stack of 2D contours into a closed 3D surface
     * (per interpolation method) via CGAL Poisson reconstruction.
     *
     *   load_stack   : read originals from a directory (height = index * dz mm).
     *   densify      : insert intermediate contours with M1/M2/M3.
     *   poisson      : assemble an oriented point cloud (in-plane outward normal
     *                  seed + jet_estimate_normals + mst_orient_normals, sign
     *                  repaired by the seed) and run
     *                  poisson_surface_reconstruction_delaunay.
     *   Delaunay does real work inside the Poisson solver (3D Delaunay of the
     *   cloud); the optional Boissonnat band exposes a 2D-Delaunay/tiling
     *   baseline between adjacent original slices.
     */
    class Reconstructor
    {
    public:
      using TKernel   = ContourInterpolator::TKernel;       // EPIC
      using Point2    = TKernel::Point_2;
      using Point3    = TKernel::Point_3;
      using Vector3   = TKernel::Vector_3;
      using Mesh      = CGAL::Surface_mesh< Point3 >;
      using TContour  = ContourInterpolator::TContour;
      using InterpKind = ContourInterpolator::InterpKind;

      struct Slice { TContour contour; double z; };

    public:
      // Originals from a directory of *.obj. Height z = slice_index * dz (mm);
      // if the file name carries no index, the ordinal position is used.
      static std::vector< Slice > load_stack( const std::string& dir, double dz );

      // Insert K-1 intermediate contours per gap using the chosen method.
      static std::vector< Slice > densify( const std::vector< Slice >& orig,
                                           InterpKind kind, int K );

      // Oriented point cloud -> Poisson -> closed Surface_mesh. Returns false
      // if the cloud is too small or the solver produced no faces.
      static bool poisson_reconstruct( const std::vector< Slice >& stack,
                                       Mesh& out );

      // Write .off / .ply (chosen by extension) via CGAL::IO::write_polygon_mesh.
      static bool write_mesh( const std::string& fname, const Mesh& m );

      // closed (is_closed) and orientable (orient_to_bound_a_volume then
      // does_bound_a_volume). Reorients m in place.
      static bool verify_closed_orientable( Mesh& m, bool& closed,
                                            bool& bounds_volume );

      static double volume( const Mesh& m );   // PMP::volume (signed)

      // OPTIONAL baseline: explicit capped contour-stitching band between
      // adjacent original slices (Fuchs-Kedem-Uselton-style correspondence
      // tiling + centroid caps), for comparison against Poisson.
      static bool boissonnat_band( const std::vector< Slice >& orig, Mesh& out );
    };
  } // end namespace Final
} // end namespace pujCGAL

#include "Reconstructor.hxx"

#endif // __Reconstructor__h__

// eof - Reconstructor.h
