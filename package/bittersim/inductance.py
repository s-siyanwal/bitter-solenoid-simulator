"""Inductance and stored energy of axisymmetric windings.

E8  Maxwell mutual inductance of coaxial loops:
    M = mu0 sqrt(a b) [(2/k - k) K(k) - (2/k) E(k)],  k^2 = 4ab/((a+b)^2+u^2)
E9  Self-inductance of a winding with N turns and normalised current
    density profile: L = (N/NI)^2 * iint J(r)J(r') 2 int_0^len (len-u) M(r,r',u) du dr dr'
E10 Nagaoka/Lorenz current-sheet inductance (exact):
    L = mu0 pi R^2 N^2/len * K_N,
    K_N = 4/(3 pi k') [ (k'^2/k^2)(K(k)-E(k)) + E(k) - k ],  k^2 = 4R^2/(4R^2+len^2)
Energy W = L I^2 / 2.
"""
import numpy as np
from scipy.special import ellipk, ellipe
from scipy.integrate import quad
from .constants import MU_0


def mutual_loops(a, b, u):
    m = np.minimum(4.0 * a * b / ((a + b) ** 2 + u ** 2), 1.0 - 1e-15)
    k = np.sqrt(m)
    return MU_0 * np.sqrt(a * b) * ((2.0 / k - k) * ellipk(m) - 2.0 / k * ellipe(m))


def nagaoka_inductance(R, length, N):
    m = 4.0 * R ** 2 / (4.0 * R ** 2 + length ** 2)
    k = np.sqrt(m)
    kp = np.sqrt(1.0 - m)
    KN = 4.0 / (3.0 * np.pi * kp) * ((kp ** 2 / m) * (ellipk(m) - ellipe(m)) + ellipe(m) - k)
    return MU_0 * np.pi * R ** 2 * N ** 2 / length * KN


def _pair(a, b, length):
    f = lambda u: (length - u) * mutual_loops(a, b, u)
    # M ~ 1e-7 H: absolute tolerance must be disabled (epsabs=0)
    split = min(length, 4.0 * max(a, b))
    val1, _ = quad(f, 0.0, split, limit=400, epsabs=0.0, epsrel=1e-10)
    val2 = 0.0
    if split < length:
        val2, _ = quad(f, split, length, limit=400, epsabs=0.0, epsrel=1e-10)
    val = val1 + val2
    return 2.0 * val


def sheet_inductance(R, length, N):
    """Thin current sheet via E9 with r = r' = R (log singularity handled by quad)."""
    return (N / length) ** 2 * _pair(R, R, length)


def coil_inductance(R1, R2, length, N, profile="bitter", nr=10):
    """E9 for a thick winding with J ~ 1/r ('bitter') or J = const ('uniform')."""
    x, w = np.polynomial.legendre.leggauss(nr)
    r = 0.5 * (R2 - R1) * x + 0.5 * (R1 + R2)
    wr = 0.5 * (R2 - R1) * w
    J = 1.0 / r if profile == "bitter" else np.ones_like(r)
    NI_per_amp_density = length * np.sum(J * wr)  # NI for unit amplitude
    tot = 0.0
    for i in range(nr):
        for j in range(i, nr):
            v = J[i] * J[j] * wr[i] * wr[j] * _pair(r[i], r[j], length)
            tot += v if i == j else 2.0 * v
    return (N / NI_per_amp_density) ** 2 * tot
