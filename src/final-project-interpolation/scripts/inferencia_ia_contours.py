import os
import time
import torch
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import imageio.v3 as iio
import cv2
from pathlib import Path
from skimage.draw import polygon

# ============================================================
# 1. WARP Y ARQUITECTURAS (Tomadas de experimento_mitad_resolucion.py)
# ============================================================
def warp_backward_diferenciable(img, flow):
    img, flow = img.float(), flow.float()
    B, C, H, W = img.size()
    _, _, Hf, Wf = flow.size()
    if H != Hf or W != Wf:
        flow = F.interpolate(flow, size=(H, W), mode="bilinear", align_corners=False)
    xs = torch.linspace(-1, 1, W, device=img.device).view(1, 1, 1, W).expand(B, -1, H, -1)
    ys = torch.linspace(-1, 1, H, device=img.device).view(1, 1, H, 1).expand(B, -1, -1, W)
    grid = torch.cat([xs, ys], dim=1)
    flow_norm = torch.zeros_like(flow)
    flow_norm[:, 0, :, :] = flow[:, 0, :, :] / ((W - 1) / 2.0)
    flow_norm[:, 1, :, :] = flow[:, 1, :, :] / ((H - 1) / 2.0)
    grid = (grid + flow_norm).permute(0, 2, 3, 1)
    return F.grid_sample(img, grid, mode='bilinear', padding_mode='border', align_corners=True)

class BaselineIFBlock(nn.Module):
    def __init__(self, c_in, c_base=64):
        super().__init__()
        self.conv0 = nn.Sequential(nn.Conv2d(c_in, c_base//2, 3, 2, 1), nn.PReLU())
        self.deep = nn.Sequential(nn.Conv2d(c_base//2, c_base, 3, 2, 1), nn.PReLU(), nn.Conv2d(c_base, c_base, 3, 1, 1), nn.PReLU(), nn.Conv2d(c_base, c_base, 3, 1, 1), nn.PReLU())
        self.proy = nn.ConvTranspose2d(c_base, 5, 4, 2, 1)
    def forward(self, x, prev_flow=None, scale=1.0):
        x = x.float()
        if scale != 1.0: x = F.interpolate(x, scale_factor=1.0/scale, mode="bilinear", align_corners=False)
        if prev_flow is not None:
            prev_flow = F.interpolate(prev_flow.float(), scale_factor=1.0/scale, mode="bilinear", align_corners=False) * (1.0/scale)
            x = torch.cat([x, prev_flow], dim=1)
        f = self.conv0(x)
        out = self.proy(self.deep(f))
        if scale != 1.0: out = F.interpolate(out, scale_factor=scale, mode="bilinear", align_corners=False) * scale
        return out

class BaselineIFNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.b0 = BaselineIFBlock(6, 240)
        self.b1 = BaselineIFBlock(17, 150)
        self.b2 = BaselineIFBlock(17, 90)
        self.b3 = BaselineIFBlock(17, 60)
    def forward(self, m0, m1):
        m0, m1 = m0.float()[:, :3], m1.float()[:, :3]
        B, C, H, W = m0.shape
        x = torch.cat([m0, m1], dim=1)
        r0 = F.interpolate(self.b0(x, None, 4.0), size=(H, W), mode="bilinear", align_corners=False)
        flow, mask = r0[:, :4], torch.sigmoid(r0[:, 4:5])
        w0, w1 = warp_backward_diferenciable(m0, flow[:, 0:2]), warp_backward_diferenciable(m1, flow[:, 2:4])
        r1 = F.interpolate(self.b1(torch.cat([m0, m1, w0, w1, mask], dim=1), flow, 2.0), size=(H, W), mode="bilinear", align_corners=False)
        flow, mask = flow + r1[:, :4], torch.sigmoid(r1[:, 4:5])
        w0, w1 = warp_backward_diferenciable(m0, flow[:, 0:2]), warp_backward_diferenciable(m1, flow[:, 2:4])
        r2 = F.interpolate(self.b2(torch.cat([m0, m1, w0, w1, mask], dim=1), flow, 1.0), size=(H, W), mode="bilinear", align_corners=False)
        flow, mask = flow + r2[:, :4], torch.sigmoid(r2[:, 4:5])
        w0, w1 = warp_backward_diferenciable(m0, flow[:, 0:2]), warp_backward_diferenciable(m1, flow[:, 2:4])
        r3 = F.interpolate(self.b3(torch.cat([m0, m1, w0, w1, mask], dim=1), flow, 0.5), size=(H, W), mode="bilinear", align_corners=False)
        flow, mask = flow + r3[:, :4], torch.sigmoid(r3[:, 4:5])
        f0, f1 = warp_backward_diferenciable(m0, flow[:, 0:2]), warp_backward_diferenciable(m1, flow[:, 2:4])
        return f0 * mask + f1 * (1.0 - mask), flow

# ============================================================
# 2. FUNCIONES DE APOYO PARA CONTORNOS Y PADDING
# ============================================================
def obj_to_mask(obj_path: Path, shape=(240, 240)):
    """Convierte un contorno .obj en una máscara de 3 canales (RGB-like) para la red."""
    coords = []
    with open(obj_path, 'r') as f:
        for line in f:
            if line.startswith('v '):
                _, x, y = line.split()
                coords.append([float(y), float(x)]) # skimage usa (row, col)
    
    coords = np.array(coords)
    mask = np.zeros(shape, dtype=np.uint8)
    
    if len(coords) > 2:
        rr, cc = polygon(coords[:, 0], coords[:, 1], shape)
        mask[rr, cc] = 255
        
    # Convertir a 3 canales para que sea compatible con BaselineIFNet (espera 6 canales concat)
    mask_3c = np.stack([mask]*3, axis=-1)
    return mask_3c

def pad_image(img_tensor, multiplier=32):
    _, _, h, w = img_tensor.size()
    ph = ((h - 1) // multiplier + 1) * multiplier
    pw = ((w - 1) // multiplier + 1) * multiplier
    padding = (0, pw - w, 0, ph - h)
    return F.pad(img_tensor, padding), (h, w)

# ============================================================
# 3. MOTOR DE INFERENCIA PARA PARES DE SLICES (.OBJ)
# ============================================================
def generar_interpolacion_ia(ruta_obj_a, ruta_obj_b, ruta_png_salida, ruta_modelo):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n[*] INICIANDO INFERENCIA IA (ESCALA 1.0)")
    print(f"[*] Dispositivo: {device}")
    
    # Cargar el modelo Baseline
    model = BaselineIFNet().to(device)
    try:
        model.load_state_dict(torch.load(ruta_modelo, map_location=device))
    except Exception as e:
        print(f"[!] Error al cargar el modelo: {e}")
        return
        
    model.eval()

    # 1. Convertir geometría a píxeles
    mask_a = obj_to_mask(Path(ruta_obj_a))
    mask_b = obj_to_mask(Path(ruta_obj_b))
    
    # 2. Preparar tensores (B, C, H, W) normalizados a [0, 1]
    t0_orig = torch.from_numpy(mask_a).permute(2, 0, 1).unsqueeze(0).float().to(device) / 255.0
    t1_orig = torch.from_numpy(mask_b).permute(2, 0, 1).unsqueeze(0).float().to(device) / 255.0
    
    # Aplicar padding para que las dimensiones sean divisibles por 32
    t0_pad, (h_orig, w_orig) = pad_image(t0_orig)
    t1_pad, _ = pad_image(t1_orig)

    print(f"[*] Procesando máscaras de tamaño: {w_orig}x{h_orig}")
    tiempo_inicio = time.time()

    with torch.no_grad():
        with torch.autocast(device_type="cuda" if "cuda" in str(device) else "cpu", dtype=torch.float16):
            # 3. Inferencia de la Red
            pred_pad, flow_pad = model(t0_pad, t1_pad)
            
            # Quitar el padding
            pred = pred_pad[:, :, :h_orig, :w_orig]

        # 4. Formatear la salida para guardar
        pred_np = pred.squeeze(0).permute(1, 2, 0).cpu().numpy()
        pred_np = (np.clip(pred_np, 0, 1) * 255.0).astype(np.uint8)
        
        # Como es una máscara, podemos tomar solo el primer canal para guardarlo en escala de grises
        pred_gray = pred_np[:, :, 0]
        
        # Guardar resultado
        iio.imwrite(ruta_png_salida, pred_gray)
    
    tiempo_total = time.time() - tiempo_inicio
    print(f"[*] ¡CONVERSIÓN FINALIZADA en {tiempo_total:.3f} segundos!")
    print(f"[*] Imagen generada en: {ruta_png_salida}\n")

if __name__ == "__main__":
    # Ajusta estas rutas a tu estructura del proyecto
    ROOT_DIR = Path(__file__).resolve().parents[1]
    CONTOURS_DIR = ROOT_DIR / "data" / "contours" / "BraTS-GLI-00000-000"
    
    OBJ_A = CONTOURS_DIR / "slice_0070.obj"
    OBJ_B = CONTOURS_DIR / "slice_0071.obj"
    PNG_OUT = CONTOURS_DIR / "baseline_mid_0070_0071.png"
    
    # Asegúrate de colocar la ruta correcta a tu .pth
    MODELO_PATH = ROOT_DIR / "modelos" / "modelo_baseline_100_definitivo.pth" 
    
    if OBJ_A.exists() and OBJ_B.exists():
        generar_interpolacion_ia(str(OBJ_A), str(OBJ_B), str(PNG_OUT), str(MODELO_PATH))
    else:
        print(f"[!] No se encontraron los archivos .obj en: {CONTOURS_DIR}")