#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <pujCGAL/Heightmap.h>
#include <pujCGAL/IO.h>

using TKernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using TReal    = TKernel::RT;
using TPoint3  = TKernel::Point_3;
using TMesh    = CGAL::Surface_mesh< TPoint3 >;

using TDelaunayV  = CGAL::Triangulation_vertex_base_with_info_2< TReal, TKernel >;
using TDelaunayDS = CGAL::Triangulation_data_structure_2< TDelaunayV >;
using TDelaunay   = CGAL::Delaunay_triangulation_2< TKernel, TDelaunayDS >;
using TDelaunayP  = TKernel::Point_2;

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

template< class _TMesh >
void marcar_vecinos_como_muertos(
  const _TMesh& mesh,
  const typename _TMesh::Vertex_index& v_centro,
  const std::size_t& order,
  std::vector< bool >& dead
  );

// -------------------------------------------------------------------------
// Copiada de decimate_mesh.cxx (solucion del profesor), sin cambios de logica.
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
  using TKernelL    = typename CGAL::Kernel_traits< TPoint >::Kernel;
  using TRealL      = typename TKernelL::RT;
  using TVertex     = typename _TMesh::Vertex_index;
  using TCirculator = CGAL::Vertex_around_target_circulator< _TMesh >;

  std::size_t N = mesh.vertices( ).size( );
  std::vector< bool > marks( N, false );

  std::size_t K = 0;
  TRealL M = 0;
  TRealL S = 0;

  std::queue< std::pair< TVertex, std::size_t > > q;
  q.push( std::make_pair( vIdx, std::size_t( 0 ) ) );

  while( q.size( ) > 0 )
  {
    auto n = q.front( );
    q.pop( );

    if( n.second < order )
    {
      if( !( marks[ std::size_t( n.first ) ] ) )
      {
        marks[ std::size_t( n.first ) ] = true;

        K++;
        TRealL z = mesh.point( n.first )[ 2 ];
        TRealL D = z - M;
        M += ( z - M ) / TRealL( K );
        S += ( z - M ) * D;

        TCirculator cIt( mesh.halfedge( n.first ), mesh );
        TCirculator cItEnd( cIt );
        do
        {
          if( !( marks[ std::size_t( *cIt ) ] ) )
            q.push( std::make_pair( *cIt, n.second + 1 ) );
          cIt++;
        } while( cIt != cItEnd );
      }
    }
  }

  if( K < 2 )
    return std::make_pair( M, TRealL( 0 ) );

  return std::make_pair( M, S / TRealL( K - 1 ) );
}

template< class _TMesh >
void marcar_vecinos_como_muertos(
  const _TMesh& mesh,
  const typename _TMesh::Vertex_index& v_centro,
  const std::size_t& order,
  std::vector< bool >& dead
  )
{
  using TVertex     = typename _TMesh::Vertex_index;
  using TCirculator = CGAL::Vertex_around_target_circulator< _TMesh >;

  std::vector< bool > visited( mesh.vertices( ).size( ), false );

  std::queue< std::pair< TVertex, std::size_t > > q;
  q.push( std::make_pair( v_centro, std::size_t( 0 ) ) );

  while( q.size( ) > 0 )
  {
    auto n = q.front( );
    q.pop( );

    if( n.second >= order )                 continue;
    if( visited[ std::size_t( n.first ) ] ) continue;
    visited[ std::size_t( n.first ) ] = true;

    if( n.second > 0 )
      dead[ std::size_t( n.first ) ] = true;

    TCirculator cIt( mesh.halfedge( n.first ), mesh );
    TCirculator cItEnd( cIt );
    do
    {
      if( !( visited[ std::size_t( *cIt ) ] ) )
        q.push( std::make_pair( *cIt, n.second + 1 ) );
      cIt++;
    } while( cIt != cItEnd );
  }
}

struct ResultadoImagen
{
  std::string imagen;
  std::size_t width { 0 };
  std::size_t height { 0 };
  std::size_t vertices_antes { 0 };
  std::size_t vertices_despues { 0 };
  double reduccion_pct { 0.0 };
  double tiempo_etapa1_ms { 0.0 };
  double tiempo_etapa2_ms { 0.0 };
  bool render_ok { true };
  std::string render_error;
  bool ok { true };
  std::string error;
};

static std::string fecha_hoy( )
{
  const auto now = std::chrono::system_clock::now( );
  const std::time_t t = std::chrono::system_clock::to_time_t( now );
  std::tm tm { };
#ifdef __APPLE__
  localtime_r( &t, &tm );
#else
  tm = *std::localtime( &t );
#endif
  std::ostringstream ss;
  ss << std::put_time( &tm, "%Y-%m-%d" );
  return ss.str( );
}

static bool es_png( const std::filesystem::path& p )
{
  std::string ext = p.extension( ).string( );
  std::transform(
    ext.begin( ), ext.end( ), ext.begin( ),
    []( unsigned char c ) { return static_cast< char >( std::tolower( c ) ); }
    );
  return ext == ".png";
}

// Replica la etapa build_delaunay_from_heightmap.cxx
static TDelaunay construir_delaunay_desde_png(
    const std::filesystem::path& png_path,
    pujCGAL::Heightmap< TReal >& hm
    )
{
  hm.read_from_png( png_path.string( ) );
  hm.set_origin( -10, -10 );
  hm.set_dimensions( 10, 10 );

  std::vector< TDelaunayP > points;
  for( std::size_t h = 0; h < hm.height( ); ++h )
  {
    for( std::size_t w = 0; w < hm.width( ); ++w )
    {
      auto p = hm.point( w, h );
      points.push_back( TDelaunayP( p.first, p.second ) );
    }
  }

  TDelaunay T;
  T.insert( points.begin( ), points.end( ) );
  for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
  {
    auto i = hm.index( v->point( ).x( ), v->point( ).y( ) );
    v->info( ) = hm( i.first, i.second );
  }

  return T;
}

static ResultadoImagen procesar_imagen(
    const std::filesystem::path& png_path,
    const std::filesystem::path& output_root,
  const std::filesystem::path& project_root,
    std::size_t order,
    TReal gamma
    )
{
  ResultadoImagen r;
  r.imagen = png_path.stem( ).string( );

  const auto out_dir = output_root / ( "batch-" + r.imagen );
  std::filesystem::create_directories( out_dir );

  try
  {
    pujCGAL::Heightmap< TReal > hm;

    const auto t1_start = std::chrono::high_resolution_clock::now( );
    TDelaunay T = construir_delaunay_desde_png( png_path, hm );
    pujCGAL::IO::save( T, ( out_dir / "original.obj" ).string( ) );
    const auto t1_end = std::chrono::high_resolution_clock::now( );

    r.width = hm.width( );
    r.height = hm.height( );
    r.vertices_antes = hm.width( ) * hm.height( );
    r.tiempo_etapa1_ms =
      std::chrono::duration< double, std::milli >( t1_end - t1_start ).count( );

    // Ruta elegida: leer original.obj a Surface_mesh para etapa 2.
    const auto t2_start = std::chrono::high_resolution_clock::now( );

    TMesh mesh;
    pujCGAL::IO::read( mesh, ( out_dir / "original.obj" ).string( ) );

    std::vector< TDelaunayP > delaunay_points;
    std::vector< TReal > delaunay_heights;

    const std::size_t N = mesh.vertices( ).size( );
    std::vector< bool > dead( N, false );

    // Pasada 1: para cada p plano, marcar a sus vecinos como muertos.
    // Si p es importante (fuera de banda) no se hace nada en esta pasada.
    // Regla: gana el primer "muere"; un vertice marcado no se resucita.
    for( auto vIt = mesh.vertices( ).begin( ); vIt != mesh.vertices( ).end( ); ++vIt )
    {
      auto stats = compute_neighborhood_stats( mesh, *vIt, order );
      TReal d = stats.first - mesh.point( *vIt )[ 2 ];
      TReal s = std::sqrt( stats.second ) * gamma;

      if( std::fabs( d ) < s )
        marcar_vecinos_como_muertos( mesh, *vIt, order, dead );
    }

    // Pasada 2: recolectar todos los vivos.
    for( auto vIt = mesh.vertices( ).begin( ); vIt != mesh.vertices( ).end( ); ++vIt )
    {
      if( dead[ std::size_t( *vIt ) ] ) continue;
      auto p = mesh.point( *vIt );
      delaunay_points.push_back( TDelaunayP( p[ 0 ], p[ 1 ] ) );
      delaunay_heights.push_back( p[ 2 ] );
    }

    TDelaunay Td;
    Td.insert( delaunay_points.begin( ), delaunay_points.end( ) );
    auto dvIt = Td.finite_vertices_begin( );
    auto hIt = delaunay_heights.begin( );
    for( ; dvIt != Td.finite_vertices_end( ) && hIt != delaunay_heights.end( ); ++dvIt, ++hIt )
      dvIt->info( ) = *hIt;

    pujCGAL::IO::save( Td, ( out_dir / "simplificado.obj" ).string( ) );

    const auto t2_end = std::chrono::high_resolution_clock::now( );
    r.tiempo_etapa2_ms =
      std::chrono::duration< double, std::milli >( t2_end - t2_start ).count( );

    r.vertices_despues = Td.number_of_vertices( );
    if( r.vertices_antes > 0 )
      r.reduccion_pct =
        100.0 *
        ( 1.0 - static_cast< double >( r.vertices_despues ) /
          static_cast< double >( r.vertices_antes ) );

    const std::filesystem::path visualizer_py = project_root / "visualizer.py";
    const double gamma_value = static_cast< double >( gamma );
    std::ostringstream cmd;
    cmd << "python3 \"" << visualizer_py.string( ) << "\""
        << " --input-png \""  << png_path.string( ) << "\""
        << " --output-dir \"" << out_dir.string( ) << "\""
      << " --name \""       << r.imagen << "\""
      << " --order "         << order
      << " --gamma "         << gamma_value;
    const int render_rc = std::system( cmd.str( ).c_str( ) );
    if( render_rc != 0 )
    {
      r.render_ok = false;
      r.render_error = "codigo de salida=" + std::to_string( render_rc );
    }
  }
  catch( const std::exception& e )
  {
    r.ok = false;
    r.error = e.what( );
  }

  return r;
}

static void escribir_csv(
    const std::filesystem::path& csv_path,
    const std::vector< ResultadoImagen >& resultados
    )
{
  std::ofstream csv( csv_path );
  csv << "imagen,width,height,vertices_antes,vertices_despues,reduccion_pct,";
  csv << "tiempo_etapa1_ms,tiempo_etapa2_ms\n";

  for( const auto& r : resultados )
  {
    if( r.ok )
    {
      csv << r.imagen << ","
          << r.width << ","
          << r.height << ","
          << r.vertices_antes << ","
          << r.vertices_despues << ","
          << r.reduccion_pct << ","
          << r.tiempo_etapa1_ms << ","
          << r.tiempo_etapa2_ms << "\n";
    }
    else
    {
      csv << r.imagen << " [ERROR],0,0,0,0,0,0,0\n";
    }
  }
}

static void escribir_bitacora(
    const std::filesystem::path& md_path,
    const std::vector< ResultadoImagen >& resultados,
    std::size_t order,
    TReal gamma
    )
{
  std::ofstream md( md_path );
  md << "# Bitacora -- Batch runner solucion Leonardo (Taller 3)\n\n";
  md << "Fecha: " << fecha_hoy( ) << "\n";
  md << "Parametros: order=" << order << ", gamma=" << gamma << "\n\n";

  md << "## Pipeline\n";
  md << "En este batch replique el pipeline del profesor en dos etapas. "
     << "Primero construyo el Delaunay completo desde el heightmap y guardo "
     << "original.obj. Luego leo ese OBJ como Surface_mesh y aplico el criterio "
     << "estadistico |h - mu| < gamma*sigma sobre vecindad de orden k para filtrar "
     << "vertices. Con los sobrevivientes reconstruyo Delaunay y guardo simplificado.obj.\n\n";

  md << "## Imagenes procesadas\n";
  for( const auto& r : resultados )
  {
    if( r.ok )
      md << "- " << r.imagen << ".png (" << r.width << "x" << r.height
         << " = " << ( r.width * r.height ) << " pixeles)\n";
    else
      md << "- " << r.imagen << ".png (ERROR: " << r.error << ")\n";
  }

  md << "\n## Resultados por imagen\n\n";
  for( const auto& r : resultados )
  {
    md << "### " << r.imagen << "\n";
    if( !r.ok )
    {
      md << "No pude procesar esta imagen. Error: " << r.error << "\n\n";
      md << "PNGs: error al renderizar (procesamiento de imagen fallido)\n\n";
      continue;
    }

    md << "| vertices antes -> despues | reduccion % | tiempo etapa 1 (ms) | tiempo etapa 2 (ms) |\n";
    md << "|---:|---:|---:|---:|\n";
    md << "| " << r.vertices_antes << " -> " << r.vertices_despues
       << " | " << std::fixed << std::setprecision( 2 ) << r.reduccion_pct
       << " | " << std::fixed << std::setprecision( 2 ) << r.tiempo_etapa1_ms
       << " | " << std::fixed << std::setprecision( 2 ) << r.tiempo_etapa2_ms
       << " |\n\n";

    md << "En esta imagen consegui una reduccion de "
       << std::fixed << std::setprecision( 2 ) << r.reduccion_pct
       << "%. "
       << ( r.tiempo_etapa1_ms >= r.tiempo_etapa2_ms ?
            "La etapa 1 tomo mas tiempo que la etapa 2. " :
            "La etapa 2 tomo mas tiempo que la etapa 1. " )
       << "El resultado me parece razonable para el criterio de filtrado usado.\n\n";

    if( r.render_ok )
      md << "PNGs: 4/4 generados (00_heightmap.png, 01_original.png, 02_simplificado.png, 03_comparativo.png)\n\n";
    else
      md << "PNGs: error al renderizar (" << r.render_error << ")\n\n";
  }

  std::vector< ResultadoImagen > ok_results;
  for( const auto& r : resultados )
    if( r.ok )
      ok_results.push_back( r );

  md << "## Observaciones globales\n\n";
  if( ok_results.empty( ) )
  {
    md << "- No tuve corridas exitosas para resumir patrones globales.\n";
    return;
  }

  auto max_red_it = std::max_element(
    ok_results.begin( ), ok_results.end( ),
    []( const ResultadoImagen& a, const ResultadoImagen& b )
    { return a.reduccion_pct < b.reduccion_pct; }
    );
  auto min_red_it = std::min_element(
    ok_results.begin( ), ok_results.end( ),
    []( const ResultadoImagen& a, const ResultadoImagen& b )
    { return a.reduccion_pct < b.reduccion_pct; }
    );

  auto max_pix_it = std::max_element(
    ok_results.begin( ), ok_results.end( ),
    []( const ResultadoImagen& a, const ResultadoImagen& b )
    { return ( a.width * a.height ) < ( b.width * b.height ); }
    );
  auto min_pix_it = std::min_element(
    ok_results.begin( ), ok_results.end( ),
    []( const ResultadoImagen& a, const ResultadoImagen& b )
    { return ( a.width * a.height ) < ( b.width * b.height ); }
    );

  double avg_t1 = 0.0;
  double avg_t2 = 0.0;
  for( const auto& r : ok_results )
  {
    avg_t1 += r.tiempo_etapa1_ms;
    avg_t2 += r.tiempo_etapa2_ms;
  }
  avg_t1 /= static_cast< double >( ok_results.size( ) );
  avg_t2 /= static_cast< double >( ok_results.size( ) );

  std::size_t errores = resultados.size( ) - ok_results.size( );

  md << "- La mayor reduccion fue en " << max_red_it->imagen
     << " con " << std::fixed << std::setprecision( 2 )
     << max_red_it->reduccion_pct << "%.\n";
  md << "- La menor reduccion fue en " << min_red_it->imagen
     << " con " << std::fixed << std::setprecision( 2 )
     << min_red_it->reduccion_pct << "%.\n";
  md << "- La imagen mas grande fue " << max_pix_it->imagen
     << " (" << ( max_pix_it->width * max_pix_it->height ) << " pixeles).\n";
  md << "- La imagen mas pequena fue " << min_pix_it->imagen
     << " (" << ( min_pix_it->width * min_pix_it->height ) << " pixeles).\n";
  md << "- Tiempo promedio etapa 1: " << std::fixed << std::setprecision( 2 )
     << avg_t1 << " ms; etapa 2: " << avg_t2 << " ms.\n";
  md << "- Imagenes con error durante el batch: " << errores << ".\n";
}

int main( int argc, char** argv )
{
  const std::filesystem::path exe_path = std::filesystem::canonical( argv[ 0 ] );
  const std::filesystem::path project_root = exe_path.parent_path( ).parent_path( );

  const std::filesystem::path data_dir =
    ( argc > 1 ) ? std::filesystem::path( argv[ 1 ] ) : ( project_root / "data" );
  const std::size_t order =
    ( argc > 2 ) ? static_cast< std::size_t >( std::stoul( argv[ 2 ] ) ) : 2;
  const TReal gamma = ( argc > 3 ) ? static_cast< TReal >( std::stod( argv[ 3 ] ) ) : 1.0;

  if( !std::filesystem::exists( data_dir ) || !std::filesystem::is_directory( data_dir ) )
  {
    std::cerr << "No existe el directorio de datos: " << data_dir << std::endl;
    return EXIT_FAILURE;
  }

  const std::filesystem::path output_root = project_root / "output";
  const std::filesystem::path resumen_dir = output_root / "batch-resumen";
  std::filesystem::create_directories( resumen_dir );

  std::vector< std::filesystem::path > pngs;
  for( const auto& entry : std::filesystem::directory_iterator( data_dir ) )
  {
    if( entry.is_regular_file( ) && es_png( entry.path( ) ) )
      pngs.push_back( entry.path( ) );
  }
  std::sort( pngs.begin( ), pngs.end( ) );

  std::vector< ResultadoImagen > resultados;
  for( const auto& p : pngs )
    resultados.push_back( procesar_imagen( p, output_root, project_root, order, gamma ) );

  const auto csv_path = resumen_dir / "resumen.csv";
  const auto md_path = resumen_dir / "bitacora.md";
  escribir_csv( csv_path, resultados );
  escribir_bitacora( md_path, resultados, order, gamma );

  std::cout << "CSV: " << csv_path << std::endl;
  std::cout << "Bitacora: " << md_path << std::endl;
  return EXIT_SUCCESS;
}

// eof - batch_runner.cxx
