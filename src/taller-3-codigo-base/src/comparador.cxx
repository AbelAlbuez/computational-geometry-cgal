// =========================================================================
// comparador.cxx
//
// Compara distintos algoritmos de simplificacion de heightmap sobre la
// misma imagen de entrada. Para cada algoritmo mide:
//   - tiempo de ejecucion
//   - reduccion de vertices
//   - error vertical L_infinito, L2 (RMSE) y promedio
//   - angulos minimo y maximo (calidad de la malla)
//
// Salida:
//   - tabla legible por consola
//   - CSV con todos los resultados
//   - una malla .obj por algoritmo, para inspeccion visual
//
// Uso:
//   ./comparador input.png [epsilon] [orden_k] [paso_submuestreo]
// =========================================================================

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <pujCGAL/Heightmap.h>
#include <pujCGAL/IO.h>
#include <pujCGAL/Algoritmos.h>

// -------------------------------------------------------------------------
// Imprime una fila del resultado en formato legible.
// -------------------------------------------------------------------------
static void imprimir_resultado(
    const pujCGAL::ResultadoSimplificacion& r
    )
{
  std::cout << "\n[" << r.nombre_algoritmo << "]\n";
  std::cout << "  Tiempo:                " << std::fixed << std::setprecision( 1 )
            << r.tiempo_milisegundos        << " ms\n";
  std::cout << "  Vertices:              " << r.vertices_antes
            << " -> "                       << r.vertices_despues
            << "  (reduccion "              << std::setprecision( 2 )
            << r.reduccion_porcentaje       << "%)\n";
  std::cout << "  Triangulos:            " << r.triangulos_despues   << "\n";
  std::cout << "  Error L_inf vertical:  " << std::setprecision( 4 )
            << r.error_L_infinito           << "\n";
  std::cout << "  Error L2 RMSE:         " << r.error_L2_RMSE        << "\n";
  std::cout << "  Error promedio:        " << r.error_promedio       << "\n";
  std::cout << "  Angulo minimo:         " << std::setprecision( 2 )
            << r.angulo_minimo_grados       << " deg\n";
  std::cout << "  Angulo maximo:         " << r.angulo_maximo_grados << " deg\n";
}

// -------------------------------------------------------------------------
// Escribe el CSV consolidado con los resultados de todos los algoritmos.
// -------------------------------------------------------------------------
static void escribir_csv(
    const std::filesystem::path& ruta,
    const std::vector< pujCGAL::ResultadoSimplificacion >& resultados
    )
{
  std::ofstream csv( ruta );
  csv << "algoritmo,tiempo_ms,vertices_antes,vertices_despues,triangulos,"
      << "reduccion_pct,err_L_inf,err_L2_RMSE,err_promedio,"
      << "ang_min,ang_max\n";
  for( const auto& r : resultados )
    csv << r.nombre_algoritmo       << ","
        << r.tiempo_milisegundos    << ","
        << r.vertices_antes         << ","
        << r.vertices_despues       << ","
        << r.triangulos_despues     << ","
        << r.reduccion_porcentaje   << ","
        << r.error_L_infinito       << ","
        << r.error_L2_RMSE          << ","
        << r.error_promedio         << ","
        << r.angulo_minimo_grados   << ","
        << r.angulo_maximo_grados   << "\n";
}

// -------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // Tipos: igual configuracion que heightmap.cxx para que las mallas
  // de comparador y taller3 sean intercambiables.
  using TKernel          = CGAL::Exact_predicates_inexact_constructions_kernel;
  using TReal            = typename TKernel::RT;
  using TVertices        = CGAL::Triangulation_vertex_base_with_info_2< TReal, TKernel >;
  using TTriangulationDS = CGAL::Triangulation_data_structure_2< TVertices >;
  using TDelaunay        = CGAL::Delaunay_triangulation_2< TKernel, TTriangulationDS >;

  if( argc < 2 )
  {
    std::cerr << "Uso: " << argv[ 0 ]
              << " input.png [epsilon] [orden_k] [paso_submuestreo]"
              << std::endl;
    return( EXIT_FAILURE );
  } // end if

  const double epsilon = ( argc > 2 ) ? std::stod( argv[ 2 ] ) : 10.0;
  const int    orden_k = ( argc > 3 ) ? std::stoi( argv[ 3 ] ) : 2;
  const int    paso    = ( argc > 4 ) ? std::stoi( argv[ 4 ] ) : 8;

  // Lectura del heightmap (compartida entre los algoritmos).
  pujCGAL::Heightmap< TReal > hm;
  hm.read_from_png( argv[ 1 ] );
  hm.set_origin    ( -10, -10 );
  hm.set_dimensions(  10,  10 );

  std::cout << "================================================================\n";
  std::cout << "Comparador de algoritmos de simplificacion de heightmap\n";
  std::cout << "================================================================\n";
  std::cout << "Imagen:      " << argv[ 1 ]                  << "\n";
  std::cout << "Heightmap:   " << hm.width( ) << " x " << hm.height( )
            << " = " << hm.width( ) * hm.height( ) << " pixeles\n";
  std::cout << "Parametros:  epsilon=" << epsilon
            << ", orden_k="            << orden_k
            << ", paso_submuestreo="   << paso << "\n";

  // Carpeta de salida.
  const std::filesystem::path input_path( argv[ 1 ] );
  const std::string           test_name = input_path.stem( ).string( );
  const std::filesystem::path exe_path  = std::filesystem::canonical( argv[ 0 ] );
  const std::filesystem::path output_dir
    = exe_path.parent_path( ).parent_path( )
    / "output"
    / ( "comparacion-" + test_name );
  std::filesystem::create_directories( output_dir );

  std::vector< pujCGAL::ResultadoSimplificacion > resultados;

  // ---------------------------------------------------------------------
  // Algoritmo 1: submuestreo uniforme. Linea base trivial.
  // ---------------------------------------------------------------------
  {
    std::cout << "\n>> Ejecutando submuestreo uniforme..." << std::endl;
    TDelaunay T;
    auto r = pujCGAL::simplificar_submuestreo_uniforme( T, hm, paso );
    pujCGAL::calcular_metricas( r, T, hm );
    pujCGAL::IO::save( T, ( output_dir / "1_submuestreo.obj" ).string( ) );
    resultados.push_back( r );
  }

  // ---------------------------------------------------------------------
  // Algoritmo 2: decimacion L1-promedio (algoritmo del taller).
  // ---------------------------------------------------------------------
  {
    std::cout << ">> Ejecutando decimacion L1-promedio..." << std::endl;
    TDelaunay T;
    auto r = pujCGAL::simplificar_decimacion_L1( T, hm, epsilon, orden_k );
    pujCGAL::calcular_metricas( r, T, hm );
    pujCGAL::IO::save( T, ( output_dir / "2_decimacion_L1.obj" ).string( ) );
    resultados.push_back( r );
  }

  // ---------------------------------------------------------------------
  // Algoritmo 3: insercion golosa Garland-Heckbert.
  // ---------------------------------------------------------------------
  {
    std::cout << ">> Ejecutando insercion golosa Garland-Heckbert..." << std::endl;
    TDelaunay T;
    // El umbral se interpreta en el mismo espacio de altura que el
    // algoritmo del taller, asi que reutilizamos `epsilon` como cota L_inf.
    auto r = pujCGAL::simplificar_insercion_greedy( T, hm, epsilon );
    pujCGAL::calcular_metricas( r, T, hm );
    pujCGAL::IO::save( T, ( output_dir / "3_greedy.obj" ).string( ) );
    resultados.push_back( r );
  }

  // ---- Reportes ---------------------------------------------------------
  std::cout << "\n================================================================\n";
  std::cout << "Resultados\n";
  std::cout << "================================================================";
  for( const auto& r : resultados ) imprimir_resultado( r );

  const auto csv_path = output_dir / "comparacion.csv";
  escribir_csv( csv_path, resultados );
  std::cout << "\n----------------------------------------------------------------\n";
  std::cout << "CSV:    " << csv_path.string( )    << "\n";
  std::cout << "Mallas: " << output_dir.string( )  << "\n";

  return( EXIT_SUCCESS );
}

// eof - comparador.cxx
