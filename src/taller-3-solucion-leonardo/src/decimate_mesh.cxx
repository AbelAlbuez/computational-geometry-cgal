#include <algorithm>
#include <queue>
#include <utility>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <pujCGAL/IO.h>

/**
 * Mean and sample variance is computed using Welford method
 * https://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance
 */
template< class _TMesh >
std::pair<
  typename CGAL::Kernel_traits< typename _TMesh::Point >::Kernel::RT,
  typename CGAL::Kernel_traits< typename _TMesh::Point >::Kernel::RT
  >
compute_neighborhood_stats(
  const _TMesh& mesh,
  const typename _TMesh::Vertex_index& vIdx,
  const std::size_t& order
  );

// -------------------------------------------------------------------------
int main( int argc, char** argv )
{
  using TKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using TReal   = TKernel::RT;
  using TPoint  = TKernel::Point_3;
  using TMesh   = CGAL::Surface_mesh< TPoint >;
  using TVertex = TMesh::Vertex_index;

  using TDelaunayV
    =
    CGAL::Triangulation_vertex_base_with_info_2< TReal, TKernel >;
  using TDelaunayDS = CGAL::Triangulation_data_structure_2< TDelaunayV >;
  using TDelaunay   = CGAL::Delaunay_triangulation_2< TKernel, TDelaunayDS >;
  using TDelaunayP  = TKernel::Point_2;

  if( argc < 5 )
  {
    std::cerr
      << "Usage: " << argv[ 0 ] << " input.obj output.obj order gamma"
      << std::endl;
    return( EXIT_FAILURE );
  } // end if
  std::string input_fname = argv[ 1 ];
  std::string output_fname = argv[ 2 ];
  std::size_t order; std::istringstream( argv[ 3 ] ) >> order;
  TReal gamma; std::istringstream( argv[ 4 ] ) >> gamma;

  // Read input
  TMesh mesh;
  pujCGAL::IO::read( mesh, input_fname );

  // Get points that should "survive"
  std::vector< TDelaunayP > delaunay_points;
  std::vector< TReal > delaunay_heights;
  for(
    auto vIt = mesh.vertices( ).begin( );
    vIt != mesh.vertices( ).end( );
    ++vIt
    )
  {
    auto stats = compute_neighborhood_stats( mesh, *vIt, order );

    TReal d = stats.first - mesh.point( *vIt )[ 2 ];
    TReal s = std::sqrt( stats.second ) * gamma;
    if( !( std::fabs( d ) < s ) )
    {
      auto p = mesh.point( *vIt );
      delaunay_points.push_back( TDelaunayP( p[ 0 ], p[ 1 ] ) );
      delaunay_heights.push_back( p[ 2 ] );
    } // end if
  } // end for

  // Compute new delaunay mesh
  TDelaunay delaunay;
  delaunay.insert( delaunay_points.begin( ), delaunay_points.end( ) );
  auto dvIt = delaunay.finite_vertices_begin( );
  auto hIt = delaunay_heights.begin( );
  for( ; dvIt != delaunay.finite_vertices_end( ); ++dvIt, ++hIt )
    dvIt->info( ) = *hIt;

  // Save it and finish
  pujCGAL::IO::save( delaunay, output_fname );
  return( EXIT_SUCCESS );
}

// -------------------------------------------------------------------------
template< class _TMesh >
std::pair<
  typename CGAL::Kernel_traits< typename _TMesh::Point >::Kernel::RT,
  typename CGAL::Kernel_traits< typename _TMesh::Point >::Kernel::RT
  >
compute_neighborhood_stats(
  const _TMesh& mesh,
  const typename _TMesh::Vertex_index& vIdx,
  const std::size_t& order
  )
{
  using TPoint      = typename _TMesh::Point;
  using TKernel     = typename CGAL::Kernel_traits< TPoint >::Kernel;
  using TReal       = typename TKernel::RT;
  using TVertex     = typename _TMesh::Vertex_index;
  using TCirculator = CGAL::Vertex_around_target_circulator< _TMesh >;

  std::size_t N = mesh.vertices( ).size( ); // Vertex count
  std::vector< bool > marks( N, false );    // Marks structure

  // Variables for the Welford method
  std::size_t K = 0; // Data count
  TReal M = 0;       // Mean
  TReal S = 0;       // Sample variance

  // BFS (breadth-first search) control queue
  std::queue< std::pair< TVertex, std::size_t > > q;
  q.push( std::make_pair( vIdx, std::size_t( 0 ) ) );

  // Main loop
  while( q.size( ) > 0 )
  {
    auto n = q.front( );
    q.pop( );

    // Does the vertex fall inside the desired neighborhood?
    if( n.second < order )
    {
      // Has it been marked as visited?
      if( !( marks[ std::size_t( n.first ) ] ) )
      {
        // Mark it
        marks[ std::size_t( n.first ) ] = true;

        // Update stats according to Welford algorithm
        K++;
        TReal z = mesh.point( n.first )[ 2 ];
        TReal D = z - M;
        M += ( z - M ) / TReal( K );
        S += ( z - M ) * D;

        // Visit neighboring vertices
        TCirculator cIt( mesh.halfedge( n.first ), mesh );
        TCirculator cItEnd( cIt );
        do
        {
          // Enqueue it with an updated depth
          if( !( marks[ std::size_t( *cIt ) ] ) )
            q.push( std::make_pair( *cIt, n.second + 1 ) );
          cIt++;
        } while( cIt != cItEnd );
      } // end if
    } // end if
  } // end while

  return( std::make_pair( M, S / TReal( K - 1 ) ) );
}

// eof - decimate_mesh.cxx
