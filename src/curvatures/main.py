import math
import numpy as np
import matplotlib.pyplot as plt
from curvatures import I1, I2

def det2(M):
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]

def K_gauss(u, v):
    I1_mat = I1(u, v)
    I2_mat = I2(u, v)
    det_I1 = det2(I1_mat)
    if abs(det_I1) < 1e-10:
        return 0.0
    det_I2 = det2(I2_mat)
    return det_I2 / det_I1

def H_mean(u, v):
    I1_mat = I1(u, v)
    I2_mat = I2(u, v)
    E, F, G = I1_mat[0][0], I1_mat[0][1], I1_mat[1][1]
    L, M, N = I2_mat[0][0], I2_mat[0][1], I2_mat[1][1]
    det_I1 = det2(I1_mat)
    if abs(det_I1) < 1e-10:
        return 0.0
    num = -(2 * F * M - E * N - G * L)
    return num / (2 * det_I1)

# Parámetros de la grilla
u_vals = np.linspace(0.01, math.pi - 0.01, 60)
v_vals = np.linspace(0, 2 * math.pi, 120)
U, V = np.meshgrid(u_vals, v_vals, indexing='ij')

X = np.zeros_like(U)
Y = np.zeros_like(U)
Z = np.zeros_like(U)
K = np.zeros_like(U)
H = np.zeros_like(U)

for i in range(U.shape[0]):
    for j in range(U.shape[1]):
        u = U[i, j]
        v = V[i, j]
        X[i, j] = math.sin(u) * math.cos(v)
        Y[i, j] = math.sin(u) * math.sin(v)
        Z[i, j] = math.cos(u)
        K[i, j] = K_gauss(u, v)
        H[i, j] = H_mean(u, v)

# Resumen numérico
print("=== Curvatura Gaussiana K ===")
print(f"  min : {K.min():.6f}")
print(f"  max : {K.max():.6f}")
print(f"  mean: {K.mean():.6f}")
print("  (Esfera unitaria esperado: K = 1.0 en todos los puntos)")
print()
print("=== Curvatura Media H ===")
print(f"  min : {H.min():.6f}")
print(f"  max : {H.max():.6f}")
print(f"  mean: {H.mean():.6f}")
print("  (Esfera unitaria esperado: H = 1.0 en todos los puntos)")

# Visualización
fig = plt.figure(figsize=(14, 6))

# Subplot K
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
k_plot = ax1.plot_surface(X, Y, Z, facecolors=plt.cm.coolwarm((K - K.min()) / (K.max() - K.min())), rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
ax1.set_title("Curvatura Gaussiana K")
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.set_box_aspect([1, 1, 1])
m1 = plt.cm.ScalarMappable(cmap='coolwarm')
m1.set_array(K)
cbar1 = fig.colorbar(m1, ax=ax1, shrink=0.6, pad=0.1)
cbar1.set_label("K")

# Subplot H
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
h_plot = ax2.plot_surface(X, Y, Z, facecolors=plt.cm.RdYlGn((H - H.min()) / (H.max() - H.min())), rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
ax2.set_title("Curvatura Media H")
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.set_box_aspect([1, 1, 1])
m2 = plt.cm.ScalarMappable(cmap='RdYlGn')
m2.set_array(H)
cbar2 = fig.colorbar(m2, ax=ax2, shrink=0.6, pad=0.1)
cbar2.set_label("H")

plt.tight_layout()
plt.savefig("output_curvatures.png", dpi=150)
plt.close()