// =========================================================================
// Visualizador VTK — parcial 1 concavidades (6 PNG + GIF). No es visualizer.cxx.
// =========================================================================
#ifndef __concavidades_viz_h__
#define __concavidades_viz_h__

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

struct ConcavidadVizInput
{
  std::vector< std::pair< double, double > > cloud;
  std::vector< std::pair< double, double > > ch_ccw;
  std::vector< std::pair< double, double > > adj_ccw;
  std::vector< std::pair< double, double > > mesh_points;
  std::vector< std::array< std::size_t, 3 > >   mesh_tris;
  std::vector< std::pair< double, double > > barycenters;
  std::pair< double, double >                p_inf;
  std::vector< std::pair< std::size_t, std::size_t > > dual_int_edges;
  std::vector< std::pair< std::size_t, std::size_t > > dual_ext_edges;
  std::vector< bool > is_boundary_tri;
  std::vector< bool > is_bolsillo;
  std::vector< int >  pocket_component;
  int                 num_concavities { 0 };
};

void
run_concavidades_viz( const char* argv0, const ConcavidadVizInput& in );

#endif
