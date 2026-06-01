import os
import time
import torch
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import imageio.v3 as iio
from pathlib import Path
from skimage.draw import polygon
from skimage import measure

# ============================================================
# 1. WARP Y ARQUITECTURAS
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
    coords = []
    with open(obj_path, 'r') as f:
        for line in f:
            if line.startswith('v '):
                # Dividimos la línea en partes
                parts = line.split()
                # Tomamos solo X (índice 1) e Y (índice 2), ignorando la Z y el prefijo 'v'
                x = parts[1]
                y = parts[2]
                coords.append([float(y), float(x)]) 
    
    coords = np.array(coords)
    mask = np.zeros(shape, dtype=np.uint8)
    
    if len(coords) > 2:
        rr, cc = polygon(coords[:, 0], coords[:, 1], shape)
        mask[rr, cc] = 255
        
    mask_3c = np.stack([mask]*3, axis=-1)
    return mask_3c

def pad_image(img_tensor, multiplier=32):
    _, _, h, w = img_tensor.size()
    ph = ((h - 1) // multiplier + 1) * multiplier
    pw = ((w - 1) // multiplier + 1) * multiplier
    padding = (0, pw - w, 0, ph - h)
    return F.pad(img_tensor, padding), (h, w)

# ============================================================
# 3. MOTOR DE INFERENCIA Y EXTRACCIÓN A .OBJ
# ============================================================
def mask_to_obj(mask_gray, output_path, z_coord=0.0):
    # Encontrar contornos en la máscara predicha
    contours = measure.find_contours(mask_gray, 127.5)
    if not contours:
        print("[!] No se encontró ningún contorno en la predicción de la IA.")
        return
    
    # Tomar el contorno más largo (para ignorar ruido pequeño)
    contour = max(contours, key=len)
    
    with open(output_path, "w") as f:
        f.write(f"# Contorno IA Interpolado - {len(contour)} vertices\n")
        for pt in contour:
            # find_contours devuelve (fila, columna) que equivale a (y, x)
            y, x = pt
            f.write(f"v {x:.6f} {y:.6f} {z_coord:.6f}\n")
        
        n = len(contour)
        for i in range(1, n):
            f.write(f"l {i} {i+1}\n")
        if n >= 2:
            f.write(f"l {n} 1\n")

def generar_interpolacion_ia(ruta_obj_a, ruta_obj_b, ruta_obj_salida, ruta_modelo):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    model = BaselineIFNet().to(device)
    try:
        model.load_state_dict(torch.load(ruta_modelo, map_location=device))
    except Exception as e:
        print(f"[!] Error al cargar el modelo: {e}")
        return
        
    model.eval()

    mask_a = obj_to_mask(Path(ruta_obj_a))
    mask_b = obj_to_mask(Path(ruta_obj_b))
    
    t0_orig = torch.from_numpy(mask_a).permute(2, 0, 1).unsqueeze(0).float().to(device) / 255.0
    t1_orig = torch.from_numpy(mask_b).permute(2, 0, 1).unsqueeze(0).float().to(device) / 255.0
    
    t0_pad, (h_orig, w_orig) = pad_image(t0_orig)
    t1_pad, _ = pad_image(t1_orig)

    with torch.no_grad():
        with torch.autocast(device_type="cuda" if "cuda" in str(device) else "cpu", dtype=torch.float16):
            pred_pad, flow_pad = model(t0_pad, t1_pad)
            pred = pred_pad[:, :, :h_orig, :w_orig]

        pred_np = pred.squeeze(0).permute(1, 2, 0).cpu().numpy()
        pred_gray = (np.clip(pred_np, 0, 1) * 255.0).astype(np.uint8)[:, :, 0]
        
    # Extraer el contorno y guardarlo como .obj (Z=0 para la evaluación)
    mask_to_obj(pred_gray, ruta_obj_salida, z_coord=0.0)
    print(f"[*] OBJ de IA generado en: {ruta_obj_salida}")

if __name__ == "__main__":
    import argparse
    import sys
    
    ROOT_DIR = Path(__file__).resolve().parents[1]
    DEFAULT_MODEL_PATH = ROOT_DIR / "modelos" / "modelo_baseline_100_definitivo.pth"
    
    parser = argparse.ArgumentParser(description="Generar interpolacion IA a OBJ.")
    parser.add_argument("slice_a", help="Ruta al contorno .obj del slice A")
    parser.add_argument("slice_b", help="Ruta al contorno .obj del slice B")
    parser.add_argument("output_obj", help="Ruta de salida para el .obj interpolado")
    parser.add_argument("--modelo", default=str(DEFAULT_MODEL_PATH), help="Ruta al archivo .pth del modelo")
    
    args = parser.parse_args()
    generar_interpolacion_ia(args.slice_a, args.slice_b, args.output_obj, args.modelo)