// =============================================================================
// Final project: Geometric Interpolation of Tumor Contours
// =============================================================================

#include <iostream>
#include <string>
#include <cstdlib> // Para std::system

#include "ContourInterpolator.h"

int main( int argc, char** argv )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ]
              << " slice_A.obj slice_B.obj" << std::endl;
    return 1;
  }

  using TInterp = pujCGAL::Final::ContourInterpolator;

  const std::string fa = argv[ 1 ];
  const std::string fb = argv[ 2 ];

  // 1. Cargar contornos
  auto A = TInterp::read_obj( fa );
  auto B = TInterp::read_obj( fb );

  std::cout << "Slice A: " << fa << "  (" << A.size( ) << " vertices)\n";
  std::cout << "Slice B: " << fb << "  (" << B.size( ) << " vertices)\n";

  if( A.empty( ) || B.empty( ) )
  {
    std::cerr << "Error: uno de los contornos esta vacio." << std::endl;
    return 2;
  }

  // =======================================================================
  // A. Interpolación Geométrica (C++ / CGAL)
  // =======================================================================
  std::cout << "\n[*] Ejecutando Interpolacion Geometrica (CGAL)..." << std::endl;
  
  // Interpolación al 50% (t = 0.5)
  auto mid_contour = TInterp::interpolate(A, B, 0.5);
  
  // Generar nombre de salida para el .obj interpolado
  std::string obj_out = fa.substr(0, fa.find_last_of('/')) + "/geom_interp_mid.obj";
  TInterp::write_obj(obj_out, mid_contour);
  std::cout << "[*] Contorno geometrico guardado en: " << obj_out << std::endl;

  // =======================================================================
  // B. Interpolación IA Baseline (Python)
  // =======================================================================
  std::cout << "\n[*] Llamando al orquestador Python para Baseline IA..." << std::endl;
  
  // Construir el comando. Asume que se corre desde la carpeta build/ o raíz
  // y que el entorno virtual de Python está activado.
  std::string python_cmd = "python ../scripts/inferencia_ia_contours.py " + fa + " " + fb;
  
  int py_result = std::system(python_cmd.c_str());
  
  if (py_result != 0) {
      std::cerr << "[!] Error al ejecutar el script de Python. Codigo: " << py_result << std::endl;
  } else {
      std::cout << "[*] Inferencia IA completada con exito." << std::endl;
  }

  return 0;
}