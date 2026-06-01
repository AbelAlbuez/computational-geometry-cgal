# scripts/run_batch_ia.py
import subprocess
from pathlib import Path

# Configuración de rutas
ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "build" / "interpolador_ia"
CONTOURS_DIR = ROOT / "data" / "contours"

# Las 10 carpetas de tu dataset de pruebas
CASOS = [
    "BraTS-GLI-00008-100",
    "BraTS-GLI-00008-101",
    "BraTS-GLI-00008-103",
    "BraTS-GLI-00009-100",
    "BraTS-GLI-00009-101",
    "BraTS-GLI-00020-100",
    "BraTS-GLI-00020-101",
    "BraTS-GLI-00027-100",
    "BraTS-GLI-00027-101",
    "BraTS-GLI-00046-100"
]

def procesar_lote():
    if not BINARY.exists():
        print(f"[!] Error: No se encontró el ejecutable {BINARY}.")
        print("    Asegúrate de compilar tu proyecto con 'make' en la carpeta build/ primero.")
        return

    print("=" * 80)
    print("🚀 INICIANDO PRUEBAS EN LOTE - INTERPOLACIÓN IA (VFI)")
    print("=" * 80)

    casos_exitosos = 0

    for caso in CASOS:
        case_dir = CONTOURS_DIR / caso
        if not case_dir.exists():
            print(f"[-] {caso}: Directorio no encontrado. Omitiendo.")
            continue

        # Buscar todos los archivos de slices reales (ignorando los intermedios ya generados)
        slices = sorted([f for f in case_dir.glob("slice_*.obj")])

        if len(slices) < 2:
            print(f"[-] {caso}: No hay suficientes slices para interpolar. Omitiendo.")
            continue

        # Estrategia: Tomar los dos slices del medio del volumen
        # Esto garantiza que el tumor tenga un área significativa para evaluar bien a la IA
        mid_idx = len(slices) // 2
        slice_a = slices[mid_idx - 1]
        slice_b = slices[mid_idx]
        out_obj = case_dir / "ia_interp_mid.obj"

        print(f"[*] Procesando {caso} | Usando cortes: {slice_a.name} y {slice_b.name}")

        # Construir y ejecutar el comando de C++
        cmd = [str(BINARY), str(slice_a), str(slice_b), str(out_obj)]
        
        # Ejecutamos el binario y capturamos la salida
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"    ✅ Éxito: {out_obj.name} generado correctamente.")
            casos_exitosos += 1
        else:
            print(f"    ❌ Error en la inferencia:")
            print(f"       {result.stderr.strip()}")

    print("=" * 80)
    print(f"📊 RESULTADO FINAL: {casos_exitosos}/{len(CASOS)} casos procesados con éxito.")
    print("=" * 80)

if __name__ == "__main__":
    procesar_lote()