"""Magnetic field models.

Closed forms (PDF section "Ideal and Finite Thick Solenoid Approximations"):
  E1  ideal solenoid            B = mu0 n I
  E2  finite thin solenoid      Bz(z) = mu0 n I/2 [cos th1 + cos th2]
  E3  thick uniform solenoid    Bz(z) = mu0 j/2 sum_ends zeta ln[(R2+sqrt(R2^2+zeta^2))/(R1+sqrt(R1^2+zeta^2))]
  E4  Bitter (J = C/r)          Bz(z) = mu0 C/2 sum_ends [asinh(zeta/R1) - asinh(zeta/R2)]
  E5  single loop on axis       Bz = mu0 I a^2 / (2 (a^2+z^2)^{3/2})
Off-axis fields: the winding is replaced by a Gauss-Legendre set of coaxial
current loops (exact elliptic-integral loop fields), see ``CoilLoops``.
"""
import numpy as np
from .constants import MU_0
from .elliptic import field_from_loops


# ---------------------------------------------------------------- closed forms
def B_ideal(n, I):
    """E1: infinitely long solenoid."""
    return MU_0 * n * I


def B_loop_axis(a, I, z):
    """E5: on-axis field of a single loop."""
    z = np.asarray(z, dtype=float)
    return MU_0 * I * a ** 2 / (2.0 * (a ** 2 + z ** 2) ** 1.5)


def B_thin_axis(R, L, n, I, z):
    """E2: finite thin solenoid (centre at z=0, length L) on axis."""
    z = np.asarray(z, dtype=float)
    z2 = L / 2.0 - z
    z1 = L / 2.0 + z
    return 0.5 * MU_0 * n * I * (z2 / np.sqrt(R ** 2 + z2 ** 2) + z1 / np.sqrt(R ** 2 + z1 ** 2))


def _thick_term(R1, R2, zeta):
    return zeta * np.log((R2 + np.sqrt(R2 ** 2 + zeta ** 2)) / (R1 + np.sqrt(R1 ** 2 + zeta ** 2)))


def B_thick_axis(R1, R2, L, j, z):
    """E3: uniform current density j [A/m^2] thick solenoid on axis."""
    z = np.asarray(z, dtype=float)
    return 0.5 * MU_0 * j * (_thick_term(R1, R2, L / 2.0 - z) + _thick_term(R1, R2, L / 2.0 + z))


def B_thick_center(R1, R2, L, j):
    """PDF closed form at the centre: mu0 j L/2 ln[(R2+sqrt(R2^2+L^2/4))/(R1+sqrt(R1^2+L^2/4))]."""
    h = L / 2.0
    return 0.5 * MU_0 * j * L * np.log((R2 + np.sqrt(R2 ** 2 + h ** 2)) / (R1 + np.sqrt(R1 ** 2 + h ** 2)))


def _bitter_term(R1, R2, zeta):
    return np.arcsinh(zeta / R1) - np.arcsinh(zeta / R2)


def B_bitter_axis(R1, R2, L, C, z):
    """E4: Bitter coil with smeared azimuthal J(r) = C / r [C in A/m] on axis."""
    z = np.asarray(z, dtype=float)
    return 0.5 * MU_0 * C * (_bitter_term(R1, R2, L / 2.0 - z) + _bitter_term(R1, R2, L / 2.0 + z))


def B_bitter_center(R1, R2, L, C):
    """E4 at z=0: mu0 C [asinh(L/2R1) - asinh(L/2R2)]."""
    return MU_0 * C * (np.arcsinh(L / (2.0 * R1)) - np.arcsinh(L / (2.0 * R2)))


# ---------------------------------------------------------------- loop models
def _gauss(n, lo, hi, panels=1):
    x, w = np.polynomial.legendre.leggauss(n)
    edges = np.linspace(lo, hi, panels + 1)
    xs, ws = [], []
    for k in range(panels):
        a, b = edges[k], edges[k + 1]
        xs.append(0.5 * (b - a) * x + 0.5 * (a + b))
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(xs), np.concatenate(ws)


class CoilLoops(object):
    """Axisymmetric winding represented by weighted coaxial loops.

    profile: 'uniform' (J const) or 'bitter' (J = C/r) or 'thin' (R1 == R2).
    ``amp`` is j [A/m^2] for 'uniform', C [A/m] for 'bitter', n*I [A/m] for 'thin'.
    """

    def __init__(self, R1, R2, L, amp, profile="bitter", nr=12, nz=16, z_panels=4):
        self.R1, self.R2, self.L, self.amp, self.profile = R1, R2, L, amp, profile
        z, wz = _gauss(nz, -L / 2.0, L / 2.0, z_panels)
        if profile == "thin" or R2 <= R1:
            r = np.array([R1])
            wr = np.array([1.0])
            Jr = np.array([amp])  # sheet current density K = nI
        else:
            r, wr = _gauss(nr, R1, R2)
            Jr = amp / r if profile == "bitter" else amp * np.ones_like(r)
        RR, ZZ = np.meshgrid(r, z, indexing="ij")
        I = np.outer(Jr * wr, wz)
        self.a = np.ascontiguousarray(RR.ravel())
        self.z = np.ascontiguousarray(ZZ.ravel())
        self.I = np.ascontiguousarray(I.ravel())

    @property
    def ampere_turns(self):
        return float(self.I.sum())

    def field(self, rho, z):
        rho = np.ascontiguousarray(np.atleast_1d(rho).astype(float).ravel())
        z = np.ascontiguousarray(np.atleast_1d(z).astype(float).ravel())
        rho, z = np.broadcast_arrays(rho, z)
        return field_from_loops(np.ascontiguousarray(rho), np.ascontiguousarray(z), self.a, self.z, self.I)

    def field_xyz(self, pts):
        pts = np.asarray(pts, dtype=float)
        rho = np.hypot(pts[:, 0], pts[:, 1])
        Br, Bz = self.field(rho, pts[:, 2])
        with np.errstate(invalid="ignore", divide="ignore"):
            c = np.where(rho > 0, pts[:, 0] / np.where(rho > 0, rho, 1), 0.0)
            s = np.where(rho > 0, pts[:, 1] / np.where(rho > 0, rho, 1), 0.0)
        return np.column_stack([Br * c, Br * s, Bz])


def sphere_points(radius, n=400):
    """Fibonacci points on a sphere surface plus the axis poles and centre."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    th = np.pi * (1.0 + 5 ** 0.5) * i
    pts = radius * np.column_stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)])
    extra = np.array([[0, 0, radius], [0, 0, -radius], [radius, 0, 0], [0, 0, 0.0]])
    return np.vstack([pts, extra])


def homogeneity_ppm(coil, dsv_radius, n=400):
    """Peak-to-peak |B| variation over the DSV sphere surface, in ppm of |B(0)|.

    |B| is (to leading order) a harmonic-like function, so extremes over the
    ball lie on its surface; the centre is included as the reference.
    """
    pts = sphere_points(dsv_radius, n)
    B = coil.field_xyz(pts)
    Bm = np.sqrt((B ** 2).sum(axis=1))
    B0 = Bm[-1]
    return float((Bm.max() - Bm.min()) / B0 * 1e6), float(B0)
