import math


def dot(a, b):
    """Producto punto de vectores 3D."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    """Producto cruzado de vectores 3D. Retorna lista [x, y, z]."""
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def norm(a):
    """Norma euclidiana de vector 3D."""
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def norm_cross(a, b):
    """Producto cruzado normalizado (normal unitaria)."""
    c = cross(a, b)
    n = norm(c)
    if n < 1e-12:
        return [0.0, 0.0, 0.0]
    return [c[0] / n, c[1] / n, c[2] / n]


def det2(M):
    """Determinante de matriz 2x2 representada como [[a,b],[c,d]]."""
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]


class Surface:
    def r(self, u, v):
        """Posicion r(u,v). Retorna [x, y, z]."""
        raise NotImplementedError

    def ru(self, u, v):
        """Derivada parcial dr/du. Retorna [x, y, z]."""
        raise NotImplementedError

    def rv(self, u, v):
        """Derivada parcial dr/dv. Retorna [x, y, z]."""
        raise NotImplementedError

    def ruu(self, u, v):
        """Segunda derivada d2r/du2. Retorna [x, y, z]."""
        raise NotImplementedError

    def rvv(self, u, v):
        """Segunda derivada d2r/dv2. Retorna [x, y, z]."""
        raise NotImplementedError

    def ruv(self, u, v):
        """Derivada mixta d2r/dudv. Retorna [x, y, z]."""
        raise NotImplementedError

    def normal(self, u, v):
        """Normal unitaria n = (ru x rv) / |ru x rv|."""
        return norm_cross(self.ru(u, v), self.rv(u, v))

    def I1(self, u, v):
        """Primera forma fundamental [[E,F],[F,G]]."""
        du = self.ru(u, v)
        dv = self.rv(u, v)
        E = dot(du, du)
        F = dot(du, dv)
        G = dot(dv, dv)
        return [[E, F], [F, G]]

    def I2(self, u, v):
        """Segunda forma fundamental [[L,M],[M,N]]."""
        n = self.normal(u, v)
        duu = self.ruu(u, v)
        dvv = self.rvv(u, v)
        duv = self.ruv(u, v)
        L = dot(n, duu)
        M = -dot(n, duv)
        N = dot(n, dvv)
        return [[L, M], [M, N]]

    def gaussian_curvature(self, u, v):
        """K = det(I2) / det(I1)."""
        d1 = det2(self.I1(u, v))
        if abs(d1) < 1e-10:
            return 0.0
        d2 = det2(self.I2(u, v))
        return d2 / d1

    def mean_curvature(self, u, v):
        """H = -(2FM - EN - GL) / (2*det(I1))."""
        i1 = self.I1(u, v)
        i2 = self.I2(u, v)
        E = i1[0][0]
        F = i1[0][1]
        G = i1[1][1]
        L = i2[0][0]
        M = i2[0][1]
        N = i2[1][1]
        d1 = det2(i1)
        if abs(d1) < 1e-10:
            return 0.0
        return (2 * F * M - E * N - G * L) / (2 * d1)

    def compute_grid(self, u_vals, v_vals):
        """
        Evalua X, Y, Z, K, H en la grilla (u_vals x v_vals).
        Retorna X, Y, Z, K, H como arrays numpy 2D.
        """
        import numpy as np

        nu = len(u_vals)
        nv = len(v_vals)
        X = np.zeros((nu, nv))
        Y = np.zeros((nu, nv))
        Z = np.zeros((nu, nv))
        K = np.zeros((nu, nv))
        H = np.zeros((nu, nv))

        for i, u in enumerate(u_vals):
            for j, v in enumerate(v_vals):
                p = self.r(u, v)
                X[i, j] = p[0]
                Y[i, j] = p[1]
                Z[i, j] = p[2]
                K[i, j] = self.gaussian_curvature(u, v)
                H[i, j] = self.mean_curvature(u, v)

        return X, Y, Z, K, H


class Sphere(Surface):
    def __init__(self, R=1.0):
        self.R = R

    def r(self, u, v):
        R = self.R
        return [R * math.sin(u) * math.cos(v), R * math.sin(u) * math.sin(v), R * math.cos(u)]

    def ru(self, u, v):
        R = self.R
        return [R * math.cos(u) * math.cos(v), R * math.cos(u) * math.sin(v), -R * math.sin(u)]

    def rv(self, u, v):
        R = self.R
        return [-R * math.sin(u) * math.sin(v), R * math.sin(u) * math.cos(v), 0.0]

    def ruu(self, u, v):
        R = self.R
        return [-R * math.sin(u) * math.cos(v), -R * math.sin(u) * math.sin(v), -R * math.cos(u)]

    def rvv(self, u, v):
        R = self.R
        return [-R * math.sin(u) * math.cos(v), -R * math.sin(u) * math.sin(v), 0.0]

    def ruv(self, u, v):
        R = self.R
        return [-R * math.cos(u) * math.sin(v), R * math.cos(u) * math.cos(v), 0.0]


class Ellipsoid(Surface):
    def __init__(self, a=2.0, b=1.0, c=0.5):
        self.a = a
        self.b = b
        self.c = c

    def r(self, u, v):
        a = self.a
        b = self.b
        c = self.c
        return [a * math.sin(u) * math.cos(v), b * math.sin(u) * math.sin(v), c * math.cos(u)]

    def ru(self, u, v):
        a = self.a
        b = self.b
        c = self.c
        return [a * math.cos(u) * math.cos(v), b * math.cos(u) * math.sin(v), -c * math.sin(u)]

    def rv(self, u, v):
        a = self.a
        b = self.b
        return [-a * math.sin(u) * math.sin(v), b * math.sin(u) * math.cos(v), 0.0]

    def ruu(self, u, v):
        a = self.a
        b = self.b
        c = self.c
        return [-a * math.sin(u) * math.cos(v), -b * math.sin(u) * math.sin(v), -c * math.cos(u)]

    def rvv(self, u, v):
        a = self.a
        b = self.b
        return [-a * math.sin(u) * math.cos(v), -b * math.sin(u) * math.sin(v), 0.0]

    def ruv(self, u, v):
        a = self.a
        b = self.b
        return [-a * math.cos(u) * math.sin(v), b * math.cos(u) * math.cos(v), 0.0]


class Torus(Surface):
    def __init__(self, R=2.0, r=0.8):
        self.R = R
        self.r_minor = r

    def r(self, u, v):
        R = self.R
        r = self.r_minor
        return [
            (R + r * math.cos(u)) * math.cos(v),
            (R + r * math.cos(u)) * math.sin(v),
            r * math.sin(u),
        ]

    def ru(self, u, v):
        r = self.r_minor
        return [
            -r * math.sin(u) * math.cos(v),
            -r * math.sin(u) * math.sin(v),
            r * math.cos(u),
        ]

    def rv(self, u, v):
        R = self.R
        r = self.r_minor
        return [
            -(R + r * math.cos(u)) * math.sin(v),
            (R + r * math.cos(u)) * math.cos(v),
            0.0,
        ]

    def ruu(self, u, v):
        r = self.r_minor
        return [
            -r * math.cos(u) * math.cos(v),
            -r * math.cos(u) * math.sin(v),
            -r * math.sin(u),
        ]

    def rvv(self, u, v):
        R = self.R
        r = self.r_minor
        return [
            -(R + r * math.cos(u)) * math.cos(v),
            -(R + r * math.cos(u)) * math.sin(v),
            0.0,
        ]

    def ruv(self, u, v):
        r = self.r_minor
        return [
            r * math.sin(u) * math.sin(v),
            -r * math.sin(u) * math.cos(v),
            0.0,
        ]


class PlanarCircle:
    def __init__(self, R=1.0):
        self.R = R

    def gamma(self, t):
        R = self.R
        return [R * math.cos(t), R * math.sin(t)]

    def gamma_p(self, t):
        R = self.R
        return [-R * math.sin(t), R * math.cos(t)]

    def gamma_pp(self, t):
        R = self.R
        return [-R * math.cos(t), -R * math.sin(t)]

    def curvature(self, t):
        dx, dy = self.gamma_p(t)
        ddx, ddy = self.gamma_pp(t)
        num = abs(dx * ddy - dy * ddx)
        den = (dx * dx + dy * dy) ** 1.5
        if den < 1e-10:
            return 0.0
        return num / den

    def compute_grid(self, t_vals):
        import numpy as np

        x_vals = np.zeros(len(t_vals))
        y_vals = np.zeros(len(t_vals))
        kappa_vals = np.zeros(len(t_vals))

        for i, t in enumerate(t_vals):
            p = self.gamma(t)
            x_vals[i] = p[0]
            y_vals[i] = p[1]
            kappa_vals[i] = self.curvature(t)

        return x_vals, y_vals, kappa_vals
