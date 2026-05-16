import numpy as np
import matplotlib.pyplot as plt

from curvatures import Sphere, Ellipsoid, Torus, PlanarCircle


def print_summary(name, K, H):
    print(f"\n=== {name} ===")
    print(f"  K - min: {K.min():.4f}  max: {K.max():.4f}  mean: {K.mean():.4f}")
    print(f"  H - min: {H.min():.4f}  max: {H.max():.4f}  mean: {H.mean():.4f}")


def surface_panel(fig, position, X, Y, Z, values, cmap, title, label):
    ax = fig.add_subplot(2, 3, position, projection="3d")
    vmin = float(values.min())
    vmax = float(values.max())
    if abs(vmax - vmin) < 1e-12:
        norm_values = np.full_like(values, 0.5)
    else:
        norm_values = (values - vmin) / (vmax - vmin)

    ax.plot_surface(
        X,
        Y,
        Z,
        facecolors=cmap(norm_values),
        alpha=0.9,
        linewidth=0,
        antialiased=True,
        shade=False,
    )

    mappable = plt.cm.ScalarMappable(cmap=cmap)
    mappable.set_clim(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1e-12)
    mappable.set_array(values)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.65, pad=0.06)
    cbar.set_label(label)

    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_box_aspect([1, 1, 1])


# Grillas de evaluacion
u60 = np.linspace(0.05, np.pi - 0.05, 60)
v120 = np.linspace(0.0, 2 * np.pi, 120)
u_torus = np.linspace(0.0, 2 * np.pi, 80)
v_torus = np.linspace(0.0, 2 * np.pi, 80)
t_circle = np.linspace(0.0, 2 * np.pi, 200)

# Superficies
sphere = Sphere(R=1.0)
ellipsoid = Ellipsoid(a=2.0, b=1.0, c=0.5)
torus = Torus(R=2.0, r=0.8)
circle = PlanarCircle(R=1.0)

# Evaluacion en grillas
X_s, Y_s, Z_s, K_s, H_s = sphere.compute_grid(u60, v120)
X_e, Y_e, Z_e, K_e, H_e = ellipsoid.compute_grid(u60, v120)
X_t, Y_t, Z_t, K_t, H_t = torus.compute_grid(u_torus, v_torus)
x_c, y_c, kappa_c = circle.compute_grid(t_circle)

# Resumen por consola
print_summary("Esfera (R=1)", K_s, H_s)
print_summary("Elipsoide (a=2, b=1, c=0.5)", K_e, H_e)
print_summary("Toro (R=2, r=0.8)", K_t, H_t)
print("\n=== Circulo (R=1) ===")
print(
    f"  kappa - min: {kappa_c.min():.4f}  max: {kappa_c.max():.4f}  mean: {kappa_c.mean():.4f}"
)

# Figura 3D: 2 filas x 3 columnas
fig = plt.figure(figsize=(18, 10))

surface_panel(fig, 1, X_s, Y_s, Z_s, K_s, plt.cm.coolwarm, "Esfera K", "K")
surface_panel(fig, 2, X_e, Y_e, Z_e, K_e, plt.cm.coolwarm, "Elipsoide K", "K")
surface_panel(fig, 3, X_t, Y_t, Z_t, K_t, plt.cm.coolwarm, "Toro K", "K")

surface_panel(fig, 4, X_s, Y_s, Z_s, H_s, plt.cm.RdYlGn, "Esfera H", "H")
surface_panel(fig, 5, X_e, Y_e, Z_e, H_e, plt.cm.RdYlGn, "Elipsoide H", "H")
surface_panel(fig, 6, X_t, Y_t, Z_t, H_t, plt.cm.RdYlGn, "Toro H", "H")

plt.tight_layout()
plt.savefig("output_3d.png", dpi=150)
plt.close(fig)

# Figura separada: circulo coloreado por curvatura
fig2 = plt.figure(figsize=(6, 6))
ax2 = fig2.add_subplot(1, 1, 1)
scatter = ax2.scatter(x_c, y_c, c=kappa_c, cmap="coolwarm", s=18)
ax2.plot(x_c, y_c, color="black", linewidth=0.8, alpha=0.6)
ax2.set_title("Circulo 2D con curvatura kappa")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_aspect("equal", adjustable="box")
plt.colorbar(scatter, ax=ax2, label="kappa")
plt.tight_layout()
plt.savefig("output_circle.png", dpi=150)
plt.close(fig2)
