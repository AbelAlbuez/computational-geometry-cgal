# Curvaturas de superficies y curvas

Cálculo y visualización de la curvatura gaussiana $K$ y la curvatura media $H$ sobre superficies paramétricas en $\mathbb{R}^3$, y de la curvatura $\kappa$ sobre una curva plana, usando las formas fundamentales clásicas de la geometría diferencial.

## Contenido

- [curvatures.py](curvatures.py): clases base y superficies/curvas implementadas.
  - `Surface`: clase abstracta con primera forma fundamental $I$, segunda forma fundamental $II$, curvatura gaussiana $K = \det(II)/\det(I)$ y curvatura media $H = (2FM - EN - GL)/(2\det(I))$.
  - `Sphere(R)`, `Ellipsoid(a,b,c)`, `Torus(R,r)`: parametrizaciones con derivadas analíticas hasta segundo orden.
  - `PlanarCircle(R)`: curva 2D con $\kappa = |x'y'' - y'x''| / (x'^2 + y'^2)^{3/2}$.
- [main.py](main.py): script de evaluación; muestrea las superficies en una grilla $(u, v)$, calcula $K$ y $H$ y genera las figuras.
- [requirements.txt](requirements.txt): dependencias (`numpy`, `matplotlib`).

## Fundamento

Para una superficie $r(u,v)$:

$$
I = \begin{bmatrix} E & F \\ F & G \end{bmatrix},\quad
II = \begin{bmatrix} L & M \\ M & N \end{bmatrix}
$$

con $E = r_u\cdot r_u$, $F = r_u\cdot r_v$, $G = r_v\cdot r_v$, $L = n\cdot r_{uu}$, $M = -n\cdot r_{uv}$, $N = n\cdot r_{vv}$ y $n = (r_u\times r_v)/\|r_u\times r_v\|$.

$$
K = \frac{\det II}{\det I},\qquad H = \frac{2FM - EN - GL}{2\det I}
$$

## Uso

```bash
cd src/curvatures
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python main.py
```

Salida:

- Resumen por consola con min/max/media de $K$, $H$ y $\kappa$.
- `output_3d.png`: panel 2x3 con $K$ y $H$ sobre esfera, elipsoide y toro.
- `output_circle.png`: círculo plano coloreado por $\kappa$.

## Resultados esperados

- Esfera de radio $R$: $K = 1/R^2$, $H = 1/R$ (constantes).
- Toro: $K$ positiva en la cara externa, negativa en la interna, cero en los círculos superior e inferior.
- Círculo de radio $R$: $\kappa = 1/R$ constante.
