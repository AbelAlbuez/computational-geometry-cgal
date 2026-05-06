#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <pujCGAL/Algoritmos.h>
#include <pujCGAL/Heightmap.h>
#include <pujCGAL/IO.h>

namespace
{
  using TKernel          = CGAL::Exact_predicates_inexact_constructions_kernel;
  using TReal            = typename TKernel::RT;
  using TVertices        = CGAL::Triangulation_vertex_base_with_info_2< TReal, TKernel >;
  using TTriangulationDS = CGAL::Triangulation_data_structure_2< TVertices >;
  using TDelaunay        = CGAL::Delaunay_triangulation_2< TKernel, TTriangulationDS >;
  using TPoint           = TKernel::Point_2;

  struct ResultadoPorImagen
  {
    std::string imagen;
    std::size_t width { 0 };
    std::size_t height { 0 };
    pujCGAL::ResultadoSimplificacion resultado;
  };

  static void construir_malla_completa(
      TDelaunay& T,
      const pujCGAL::Heightmap< TReal >& hm
      )
  {
    std::vector< TPoint > points;
    points.reserve( hm.width( ) * hm.height( ) );
    for( std::size_t h = 0; h < hm.height( ); ++h )
      for( std::size_t w = 0; w < hm.width( ); ++w )
      {
        const auto p = hm.point( w, h );
        points.push_back( TPoint( p.first, p.second ) );
      }

    T.insert( points.begin( ), points.end( ) );
    for( auto v = T.finite_vertices_begin( ); v != T.finite_vertices_end( ); ++v )
    {
      const auto idx = hm.index( v->point( ).x( ), v->point( ).y( ) );
      v->info( ) = hm( idx.first, idx.second );
    }
  }

  static std::string hoy_yyyy_mm_dd( )
  {
    const auto now = std::chrono::system_clock::now( );
    const std::time_t t = std::chrono::system_clock::to_time_t( now );
    std::tm tm { };
#ifdef __APPLE__
    localtime_r( &t, &tm );
#else
    tm = *std::localtime( &t );
#endif
    std::ostringstream out;
    out << std::put_time( &tm, "%Y-%m-%d" );
    return out.str( );
  }

  static void escribir_csv(
      const std::filesystem::path& ruta,
      const std::vector< ResultadoPorImagen >& filas
      )
  {
    std::ofstream csv( ruta );
    csv << "imagen,algoritmo,tiempo_ms,vertices_antes,vertices_despues,triangulos,";
    csv << "reduccion_pct,err_L_inf,err_L2_RMSE,err_promedio,ang_min,ang_max\n";
    for( const auto& f : filas )
    {
      const auto& r = f.resultado;
      csv << f.imagen << ","
          << r.nombre_algoritmo << ","
          << r.tiempo_milisegundos << ","
          << r.vertices_antes << ","
          << r.vertices_despues << ","
          << r.triangulos_despues << ","
          << r.reduccion_porcentaje << ","
          << r.error_L_infinito << ","
          << r.error_L2_RMSE << ","
          << r.error_promedio << ","
          << r.angulo_minimo_grados << ","
          << r.angulo_maximo_grados << "\n";
    }
  }

  static std::string comentario_por_imagen(
      const std::vector< ResultadoPorImagen >& filas_img
      )
  {
    if( filas_img.empty( ) )
      return "No tuve resultados para esta imagen.";

    const auto best_reduccion = *std::max_element(
      filas_img.begin( ), filas_img.end( ),
      []( const auto& a, const auto& b )
      { return a.resultado.reduccion_porcentaje < b.resultado.reduccion_porcentaje; }
      );
    const auto best_rmse = *std::min_element(
      filas_img.begin( ), filas_img.end( ),
      []( const auto& a, const auto& b )
      { return a.resultado.error_L2_RMSE < b.resultado.error_L2_RMSE; }
      );

    bool hay_degenerados = false;
    bool hay_reduccion_cero = false;
    for( const auto& f : filas_img )
    {
      if( f.resultado.angulo_minimo_grados < 1.0 )
        hay_degenerados = true;
      if( f.resultado.reduccion_porcentaje <= 0.1 )
        hay_reduccion_cero = true;
    }

    std::ostringstream out;
    out << "En esta imagen, la mayor reduccion la logre con "
        << best_reduccion.resultado.nombre_algoritmo
        << " (" << std::fixed << std::setprecision( 2 )
        << best_reduccion.resultado.reduccion_porcentaje << "%). ";
    out << "El menor error (RMSE) lo obtuve con "
        << best_rmse.resultado.nombre_algoritmo
        << " (" << std::fixed << std::setprecision( 4 )
        << best_rmse.resultado.error_L2_RMSE << "). ";

    if( hay_degenerados )
      out << "Detecte angulos minimos muy bajos en al menos una malla. ";
    if( hay_reduccion_cero )
      out << "Tambien aparecio al menos un caso con reduccion casi nula.";
    if( !hay_degenerados && !hay_reduccion_cero )
      out << "No observe casos raros de reduccion nula ni angulos extremos.";

    return out.str( );
  }

  static void escribir_bitacora(
      const std::filesystem::path& ruta,
      const std::vector< std::string >& imagenes,
      const std::map< std::string, std::pair< std::size_t, std::size_t > >& dimensiones,
      const std::vector< ResultadoPorImagen >& filas,
      double epsilon,
      int orden_k,
      int paso_submuestreo
      )
  {
    std::ofstream md( ruta );
    md << "# Bitacora -- Batch runner Taller 3\n\n";
    md << "Fecha: " << hoy_yyyy_mm_dd( ) << "\n";
    md << "Parametros: epsilon=" << epsilon
       << ", orden_k=" << orden_k
       << ", paso_submuestreo=" << paso_submuestreo << "\n\n";

    md << "## Imagenes procesadas\n";
    for( const auto& nombre : imagenes )
    {
      const auto it = dimensiones.find( nombre );
      const std::size_t w = ( it == dimensiones.end( ) ) ? 0 : it->second.first;
      const std::size_t h = ( it == dimensiones.end( ) ) ? 0 : it->second.second;
      md << "- " << nombre << ".png (" << w << "x" << h
         << " = " << ( w * h ) << " pixeles)\n";
    }

    md << "\n## Resultados por imagen\n\n";
    for( const auto& nombre : imagenes )
    {
      std::vector< ResultadoPorImagen > filas_img;
      for( const auto& f : filas )
        if( f.imagen == nombre )
          filas_img.push_back( f );

      md << "### " << nombre << "\n";
      md << "| algoritmo | tiempo (ms) | vertices antes -> despues | reduccion % | L_inf | RMSE | angulo min | angulo max |\n";
      md << "|---|---:|---:|---:|---:|---:|---:|---:|\n";
      for( const auto& f : filas_img )
      {
        const auto& r = f.resultado;
        md << "| " << r.nombre_algoritmo
           << " | " << std::fixed << std::setprecision( 2 ) << r.tiempo_milisegundos
           << " | " << r.vertices_antes << " -> " << r.vertices_despues
           << " | " << std::fixed << std::setprecision( 2 ) << r.reduccion_porcentaje
           << " | " << std::fixed << std::setprecision( 4 ) << r.error_L_infinito
           << " | " << std::fixed << std::setprecision( 4 ) << r.error_L2_RMSE
           << " | " << std::fixed << std::setprecision( 2 ) << r.angulo_minimo_grados
           << " | " << std::fixed << std::setprecision( 2 ) << r.angulo_maximo_grados
           << " |\n";
      }

      md << "\n" << comentario_por_imagen( filas_img ) << "\n\n";
    }

    // Resumen global por algoritmo.
    struct Acumulado
    {
      int n { 0 };
      double tiempo { 0.0 };
      double reduccion { 0.0 };
      double rmse { 0.0 };
      double ang_min { 0.0 };
    };
    std::map< std::string, Acumulado > acc;
    for( const auto& f : filas )
    {
      auto& a = acc[ f.resultado.nombre_algoritmo ];
      ++a.n;
      a.tiempo += f.resultado.tiempo_milisegundos;
      a.reduccion += f.resultado.reduccion_porcentaje;
      a.rmse += f.resultado.error_L2_RMSE;
      a.ang_min += f.resultado.angulo_minimo_grados;
    }

    auto mejor_reduccion = acc.begin( );
    auto mejor_rmse = acc.begin( );
    auto mas_rapido = acc.begin( );
    for( auto it = acc.begin( ); it != acc.end( ); ++it )
    {
      if( it->second.reduccion / it->second.n > mejor_reduccion->second.reduccion / mejor_reduccion->second.n )
        mejor_reduccion = it;
      if( it->second.rmse / it->second.n < mejor_rmse->second.rmse / mejor_rmse->second.n )
        mejor_rmse = it;
      if( it->second.tiempo / it->second.n < mas_rapido->second.tiempo / mas_rapido->second.n )
        mas_rapido = it;
    }

    int conteo_angulos_bajos = 0;
    for( const auto& f : filas )
      if( f.resultado.angulo_minimo_grados < 1.0 )
        ++conteo_angulos_bajos;

    md << "## Observaciones globales\n\n";
    if( !acc.empty( ) )
    {
      md << "- En promedio, el algoritmo con mayor reduccion fue "
         << mejor_reduccion->first << " ("
         << std::fixed << std::setprecision( 2 )
         << ( mejor_reduccion->second.reduccion / mejor_reduccion->second.n ) << "%).\n";
      md << "- El menor error promedio (RMSE) lo obtuve con "
         << mejor_rmse->first << " ("
         << std::fixed << std::setprecision( 4 )
         << ( mejor_rmse->second.rmse / mejor_rmse->second.n ) << ").\n";
      md << "- El tiempo promedio mas bajo fue el de "
         << mas_rapido->first << " ("
         << std::fixed << std::setprecision( 2 )
         << ( mas_rapido->second.tiempo / mas_rapido->second.n ) << " ms).\n";
      md << "- Detecte " << conteo_angulos_bajos
         << " corridas con angulo minimo menor a 1 grado, un indicador de triangulos muy agudos.\n";
      md << "- En general, la comparacion confirma el compromiso esperado entre reduccion agresiva y error de reconstruccion.\n";
    }
  }

  static bool es_png( const std::filesystem::path& p )
  {
    if( !p.has_extension( ) )
      return false;
    std::string ext = p.extension( ).string( );
    std::transform( ext.begin( ), ext.end( ), ext.begin( ),
      []( unsigned char c ) { return static_cast< char >( std::tolower( c ) ); } );
    return ext == ".png";
  }
}

int main( int argc, char** argv )
{
  try
  {
    const std::filesystem::path exe_path = std::filesystem::canonical( argv[ 0 ] );
    const std::filesystem::path project_root = exe_path.parent_path( ).parent_path( );

    const std::filesystem::path data_dir =
      ( argc > 1 ) ? std::filesystem::path( argv[ 1 ] ) : ( project_root / "data" );
    const double epsilon = ( argc > 2 ) ? std::stod( argv[ 2 ] ) : 10.0;
    const int orden_k = ( argc > 3 ) ? std::stoi( argv[ 3 ] ) : 2;
    const int paso_submuestreo = ( argc > 4 ) ? std::stoi( argv[ 4 ] ) : 8;

    if( !std::filesystem::exists( data_dir ) || !std::filesystem::is_directory( data_dir ) )
    {
      std::cerr << "No existe el directorio de datos: " << data_dir << std::endl;
      return EXIT_FAILURE;
    }

    std::vector< std::filesystem::path > imagenes;
    for( const auto& entry : std::filesystem::directory_iterator( data_dir ) )
    {
      if( entry.is_regular_file( ) && es_png( entry.path( ) ) )
        imagenes.push_back( entry.path( ) );
    }
    std::sort( imagenes.begin( ), imagenes.end( ) );

    if( imagenes.empty( ) )
    {
      std::cerr << "No encontre archivos .png en: " << data_dir << std::endl;
      return EXIT_FAILURE;
    }

    const std::filesystem::path output_root = project_root / "output";
    const std::filesystem::path resumen_dir = output_root / "batch-resumen";
    std::filesystem::create_directories( resumen_dir );

    std::vector< ResultadoPorImagen > filas_globales;
    std::vector< std::string > nombres_imagenes;
    std::map< std::string, std::pair< std::size_t, std::size_t > > dimensiones;

    for( const auto& png_path : imagenes )
    {
      const std::string nombre = png_path.stem( ).string( );
      nombres_imagenes.push_back( nombre );

      pujCGAL::Heightmap< TReal > hm;
      hm.read_from_png( png_path.string( ) );
      hm.set_origin( -10, -10 );
      hm.set_dimensions( 10, 10 );
      dimensiones[ nombre ] = { hm.width( ), hm.height( ) };

      const std::filesystem::path out_dir = output_root / ( "batch-" + nombre );
      std::filesystem::create_directories( out_dir );

      // Original completo.
      {
        TDelaunay T;
        construir_malla_completa( T, hm );
        pujCGAL::IO::save( T, ( out_dir / "original.obj" ).string( ) );
      }

      // taller3_L1 -> simplificado.obj
      {
        TDelaunay T;
        auto r = pujCGAL::simplificar_decimacion_L1( T, hm, epsilon, orden_k );
        r.nombre_algoritmo = "taller3_L1";
        pujCGAL::calcular_metricas( r, T, hm );
        pujCGAL::IO::save( T, ( out_dir / "simplificado.obj" ).string( ) );
        filas_globales.push_back( { nombre, hm.width( ), hm.height( ), r } );
      }

      // submuestreo_uniforme
      {
        TDelaunay T;
        auto r = pujCGAL::simplificar_submuestreo_uniforme( T, hm, paso_submuestreo );
        r.nombre_algoritmo = "submuestreo_uniforme";
        pujCGAL::calcular_metricas( r, T, hm );
        pujCGAL::IO::save( T, ( out_dir / "1_submuestreo.obj" ).string( ) );
        filas_globales.push_back( { nombre, hm.width( ), hm.height( ), r } );
      }

      // decimacion_L1 (comparador)
      {
        TDelaunay T;
        auto r = pujCGAL::simplificar_decimacion_L1( T, hm, epsilon, orden_k );
        r.nombre_algoritmo = "decimacion_L1";
        pujCGAL::calcular_metricas( r, T, hm );
        pujCGAL::IO::save( T, ( out_dir / "2_decimacion_L1.obj" ).string( ) );
        filas_globales.push_back( { nombre, hm.width( ), hm.height( ), r } );
      }

      // insercion_greedy
      {
        TDelaunay T;
        auto r = pujCGAL::simplificar_insercion_greedy( T, hm, epsilon );
        r.nombre_algoritmo = "insercion_greedy";
        pujCGAL::calcular_metricas( r, T, hm );
        pujCGAL::IO::save( T, ( out_dir / "3_greedy.obj" ).string( ) );
        filas_globales.push_back( { nombre, hm.width( ), hm.height( ), r } );
      }

      std::cout << "Procesada: " << png_path.filename( ).string( ) << std::endl;
    }

    const auto csv_path = resumen_dir / "resumen.csv";
    const auto md_path = resumen_dir / "bitacora.md";
    escribir_csv( csv_path, filas_globales );
    escribir_bitacora(
      md_path,
      nombres_imagenes,
      dimensiones,
      filas_globales,
      epsilon,
      orden_k,
      paso_submuestreo
      );

    std::cout << "Resumen CSV: " << csv_path << std::endl;
    std::cout << "Bitacora: " << md_path << std::endl;
  }
  catch( const std::exception& e )
  {
    std::cerr << "Error: " << e.what( ) << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

// eof - batch_runner.cxx
