// =============================================================================
// Orquestador de Inteligencia Artificial para Interpolación de Contornos
// Autor: Jesús Alberto Puenayan Quiceno
// =============================================================================

#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Uso: " << argv[0] << " <slice_A.obj> <slice_B.obj> <output.obj>\n";
        return 1;
    }

    std::string sliceA = argv[1];
    std::string sliceB = argv[2];
    std::string outputObj = argv[3];

    // Llamamos directamente al entorno de Python para ejecutar la red neuronal
    std::string command = "python3 scripts/inferencia_ia_contours.py " + 
                          sliceA + " " + sliceB + " " + outputObj;
    
    std::cout << "[*] Iniciando Pipeline de IA (VFI)...\n";
    int ret = std::system(command.c_str());
    
    if (ret == 0) {
        std::cout << "[*] Inferencia y exportacion a OBJ completada con exito.\n";
    } else {
        std::cerr << "[!] Error en la inferencia IA. Codigo de salida: " << ret << "\n";
    }

    return ret;
}