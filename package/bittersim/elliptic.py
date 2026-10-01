"""Complete elliptic integrals and exact circular current-loop fields.

Two equivalent implementations:
* ``loop_field_scipy``: vectorised NumPy using scipy.special.ellipk/ellipe
  (parameter m = k^2 convention, available in SciPy 1.1).
* ``loop_field_numba`` / ``field_from_loops``: Numba nopython kernels using the
  arithmetic-geometric mean (AGM), parallelised over observation points with
  prange.  The AGM result agrees with SciPy to ~1e-15 (see tests).
"""
import math
import numpy as np
from scipy.special import ellipk, ellipe
from numba import njit, prange

MU0 = 4.0e-7 * math.pi


@njit(cache=False)
def ellipke_agm(m):
    """Return (K(m), E(m)) via AGM, m = k^2 in [0, 1)."""
    a = 1.0
    b = math.sqrt(1.0 - m)
    s = 0.5 * m
    w = 0.5
    for _ in range(60):
        c = 0.5 * (a - b)
        an = 0.5 * (a + b)
        b = math.sqrt(a * b)
        a = an
        w *= 2.0
        s += w * c * c
        if abs(c) < 1e-17:
            break
    K = math.pi / (2.0 * a)
    return K, K * (1.0 - s)


@njit(cache=False)
def loop_field_numba(a, I, rho, z):
    """(B_rho, B_z) of a loop radius a, current I, at cylindrical (rho, z) relative to loop centre."""
    if rho < 1e-12 * a:
        return 0.0, MU0 * I * a * a / (2.0 * (a * a + z * z) ** 1.5)
    ap2 = (a + rho) ** 2 + z * z
    am2 = (a - rho) ** 2 + z * z
    m = 4.0 * a * rho / ap2
    if m >= 1.0:
        return 0.0, 0.0  # on the filament itself (singular) - skipped
    K, E = ellipke_agm(m)
    pre = MU0 * I / (2.0 * math.pi * math.sqrt(ap2))
    Bz = pre * (K + (a * a - rho * rho - z * z) / am2 * E)
    Br = pre * z / rho * (-K + (a * a + rho * rho + z * z) / am2 * E)
    return Br, Bz


@njit(parallel=True, cache=False)
def field_from_loops(obs_rho, obs_z, loop_a, loop_z, loop_I):
    """Superpose exact loop fields. Memory O(M): only outputs are stored."""
    M = obs_rho.shape[0]
    N = loop_a.shape[0]
    Br = np.zeros(M)
    Bz = np.zeros(M)
    for i in prange(M):
        sr = 0.0
        sz = 0.0
        for j in range(N):
            br, bz = loop_field_numba(loop_a[j], loop_I[j], obs_rho[i], obs_z[i] - loop_z[j])
            sr += br
            sz += bz
        Br[i] = sr
        Bz[i] = sz
    return Br, Bz


def loop_field_scipy(a, I, rho, z):
    """Vectorised reference implementation with scipy.special.ellipk/ellipe."""
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    z = np.atleast_1d(np.asarray(z, dtype=float))
    rho, z = np.broadcast_arrays(rho, z)
    Br = np.zeros(rho.shape)
    Bz = np.zeros(rho.shape)
    axis = rho < 1e-12 * a
    Bz[axis] = MU0 * I * a * a / (2.0 * (a * a + z[axis] ** 2) ** 1.5)
    r = rho[~axis]
    zz = z[~axis]
    ap2 = (a + r) ** 2 + zz ** 2
    am2 = (a - r) ** 2 + zz ** 2
    m = 4.0 * a * r / ap2
    K = ellipk(m)
    E = ellipe(m)
    pre = MU0 * I / (2.0 * np.pi * np.sqrt(ap2))
    Bz[~axis] = pre * (K + (a * a - r ** 2 - zz ** 2) / am2 * E)
    Br[~axis] = pre * zz / r * (-K + (a * a + r ** 2 + zz ** 2) / am2 * E)
    return Br, Bz
