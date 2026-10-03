"""Parallel particle-level emulation that re-derives the continuum results independently.

This is still not molecular dynamics of the magnet (about 1e28 atoms). Three particle
estimators each solve the same physics as one continuum formula, so the two can be compared:

P1  Current elements (magnetostatics). The winding is sampled with N point current elements
    J dV = C dr dphi dz phi_hat. That is uniform in (r, phi, z) because J = C/r and dV = r dr dphi dz.
    Each element adds Biot-Savart dB = mu0/(4 pi) w (phi_hat x d)/|d|^3. Every sample i takes r
    with density ~ a/(r sqrt(r^2 + a^2)), a = L/2 (Halton base 2, inverse CDF via asinh; this is the
    E4 integrand) and z = r s / sqrt(1 - s^2) with s = sin(alpha) uniform in [0, s_max]
    (Halton base 3). Together these importance-sample the isocentre integrand exactly, so B0 has zero
    variance and the error elsewhere in the DSV scales with (rho/r)^2. Each sample then places a ring of K = 16 point elements at angles
    phi0 + 2 pi m/K (random offset phi0 from base 5), plus their z-mirror: 32 elements per sample.
    A plain Halton angle left 10-30 ppm of spurious tesseral error; the ring removes it, since
    angular aliasing scales as (rho/R1)^K. The same samples give NI = int J dA and P = int rho J_cu^2 dV.
    Continuum counterparts: E4 (B0), the loop model (axis and DSV), NI and E31 (P).
P2  Conduction carriers (resistivity). N classical carriers run an exact Ornstein-Uhlenbeck
    velocity process with relaxation time tau(T). The Einstein form of Green-Kubo then gives
    sigma = n e^2 D/(k T), with D taken from the mean-square displacement. tau(T) comes from a Bloch-Gruneisen lattice-scattering
    model: rho = rho_0 + A (T/Theta)^5 J5(Theta/T), with Theta_R = 343 K and RRR = 100 (ASSUMPTIONS),
    calibrated so that rho(20 C) = 1.68e-8 ohm m. Continuum counterpart: linear rho(T) with
    alpha = 0.00393 1/K.
P3  Heat random walks (conduction). Walkers follow the 2-D Bessel process
    dR = dt/(2R) + dW, starting at the cell edge b. They reflect at b and are absorbed at the hole wall a,
    with a Brownian-bridge crossing test. Feynman-Kac gives dT(b) = q_v/(2k) E[exit time].
    Continuum counterpart: E18.

Parallelism: Numba prange. Every particle draws from its own counter-based random stream
(splitmix64 keyed by seed and particle index), and every reduction runs in fixed-size chunks
followed by a serial sum. So the serial build (parallel=False) and the parallel build give
bit-identical results.
"""
import math
import time
import numpy as np
from numba import njit, prange
from .constants import MU_0, RHO_CU_20, ALPHA_CU, K_CU

E_CHARGE = 1.60217662e-19
M_E = 9.1093837e-31
K_B = 1.38064852e-23
N_CU = 8.47e28            # conduction-electron density of copper [1/m^3]
THETA_R = 343.0           # Bloch-Gruneisen temperature for Cu [K]          ASSUMPTION
RRR = 100.0               # residual-resistivity ratio rho(293 K)/rho_0       ASSUMPTION
CHUNK = 4096

_M64 = np.uint64(0xFFFFFFFFFFFFFFFF)
_G = np.uint64(0x9E3779B97F4A7C15)
_C1 = np.uint64(0xBF58476D1CE4E5B9)
_C2 = np.uint64(0x94D049BB133111EB)
_S30 = np.uint64(30)
_S27 = np.uint64(27)
_S31 = np.uint64(31)
_S11 = np.uint64(11)


@njit(cache=False)
def _mix(z):
    z = (z ^ (z >> _S30)) * _C1
    z = (z ^ (z >> _S27)) * _C2
    return z ^ (z >> _S31)


@njit(cache=False)
def _uniform(state):
    """Advance a splitmix64 state; return (new_state, u in (0,1))."""
    state = state + _G
    z = _mix(state)
    u = (np.float64(z >> _S11) + 0.5) * (1.0 / 9007199254740992.0)
    return state, u


@njit(cache=False)
def _normal(state):
    state, u1 = _uniform(state)
    state, u2 = _uniform(state)
    return state, math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


@njit(cache=False)
def _stream(seed, i):
    return _mix(np.uint64(seed) * _G + np.uint64(i) * _C1 + np.uint64(1))


@njit(cache=False)
def _halton(i, base):
    f = 1.0
    r = 0.0
    k = i
    while k > 0:
        f = f / base
        r = r + f * (k % base)
        k = k // base
    return r


# ------------------------------------------------------------------ P1 current elements
def _field_kernel(px, py, pz, R1, R2, L, C, N, out):
    K = 16                                                     # point elements per ring sample
    a = 0.5 * L
    A1 = math.asinh(a / R1)
    Z = A1 - math.asinh(a / R2)                                # r drawn with density ~ a/(r sqrt(r^2+a^2))
    base = C * (2.0 * math.pi / K) / N * 1e-7                  # mu0/(4 pi) * C dphi / N
    for j in prange(px.shape[0]):
        bx = 0.0
        by = 0.0
        bz = 0.0
        for i in range(1, N + 1):
            r = a / math.sinh(A1 - Z * _halton(i, 2))
            hr = math.sqrt(r * r + a * a)
            smax = a / hr
            sv = smax * _halton(i, 3)                          # s = sin(alpha), uniform: importance sampling
            q = math.sqrt(1.0 - sv * sv)
            z0 = r * sv / q
            wz = (Z * r * hr / a) * smax * r / (q * q * q)     # radial weight * s-range * dz/ds
            ph0 = 2.0 * math.pi * _halton(i, 5) / K
            for m in range(2 * K):
                phm = ph0 + 2.0 * math.pi * (m % K) / K
                zm = z0 if m < K else -z0
                c = math.cos(phm)
                s = math.sin(phm)
                dx = px[j] - r * c
                dy = py[j] - r * s
                dz = pz[j] - zm
                d2 = dx * dx + dy * dy + dz * dz
                inv3 = wz / (d2 * math.sqrt(d2))
                # phi_hat = (-s, c, 0); phi_hat x d = (c dz, s dz, -s dy - c dx)
                bx += c * dz * inv3
                by += s * dz * inv3
                bz += (-s * dy - c * dx) * inv3
        out[j, 0] = base * bx
        out[j, 1] = base * by
        out[j, 2] = base * bz


def _moment_kernel(R1, R2, N, out):
    """Chunked sums of 1/r and 1/r^2 over the Halton radii (for NI and P)."""
    nchunk = out.shape[0]
    for c in prange(nchunk):
        s1 = 0.0
        s2 = 0.0
        lo = c * CHUNK + 1
        hi = min(N, (c + 1) * CHUNK)
        for i in range(lo, hi + 1):
            r = R1 + (R2 - R1) * _halton(i, 2)
            s1 += 1.0 / r
            s2 += 1.0 / (r * r)
        out[c, 0] = s1
        out[c, 1] = s2


# ------------------------------------------------------------------ P2 carriers (Green-Kubo)
def _gk_kernel(tau, vth, dt, steps, seed, out):
    """out[i] = x(T)^2 for one carrier: displacement integrated (trapezoid) over an exact OU velocity path."""
    a = math.exp(-dt / tau)
    b = vth * math.sqrt(1.0 - a * a)
    for i in prange(out.shape[0]):
        st = _stream(seed, i)
        st, g = _normal(st)
        v = vth * g
        x = 0.0
        for k in range(steps):
            st, g = _normal(st)
            vn = a * v + b * g
            x += 0.5 * (v + vn) * dt
            v = vn
        out[i] = x * x


# ------------------------------------------------------------------ P3 heat random walks
def _walk_kernel(a, b, dt, seed, max_steps, out):
    sq = math.sqrt(dt)
    for i in prange(out.shape[0]):
        st = _stream(seed, i)
        r = b
        t = 0.0
        for k in range(max_steps):
            st, g = _normal(st)
            rn = r + dt / (2.0 * r) + sq * g
            if rn > b:
                rn = 2.0 * b - rn
            t += dt
            if rn <= a:
                break
            st, u = _uniform(st)
            if u < math.exp(-2.0 * (r - a) * (rn - a) / dt):    # Brownian-bridge crossing
                break
            r = rn
        out[i] = t


_ser = {k: njit(cache=False)(f) for k, f in
        (("field", _field_kernel), ("mom", _moment_kernel), ("gk", _gk_kernel), ("walk", _walk_kernel))}
_par = {k: njit(parallel=True, cache=False)(f) for k, f in
        (("field", _field_kernel), ("mom", _moment_kernel), ("gk", _gk_kernel), ("walk", _walk_kernel))}


def _k(name, parallel):
    return (_par if parallel else _ser)[name]


# ------------------------------------------------------------------ public API
def element_field(pts, R1, R2, L, C, N, parallel=True):
    """P1: B at points (M x 3) from N Halton ring samples (32 point current elements each)."""
    pts = np.ascontiguousarray(np.asarray(pts, dtype=float))
    out = np.zeros((pts.shape[0], 3))
    _k("field", parallel)(pts[:, 0].copy(), pts[:, 1].copy(), pts[:, 2].copy(),
                          float(R1), float(R2), float(L), float(C), int(N), out)
    return out


def element_moments(R1, R2, L, C, N, rho, lam, parallel=True):
    """P1: NI = C L (R2-R1) <1/r>,  P = rho C^2/lam * 2 pi L (R2-R1) <1/r> (copper J = C/(lam r))."""
    nchunk = (int(N) + CHUNK - 1) // CHUNK
    out = np.zeros((nchunk, 2))
    _k("mom", parallel)(float(R1), float(R2), int(N), out)
    m1 = float(np.sum(out[:, 0])) / N
    # J_avg = C/r and dV = r dr dphi dz, so both NI and P integrate 1/r over a uniform radius
    return {"NI": C * L * (R2 - R1) * m1, "P": rho * C ** 2 / lam * 2 * math.pi * L * (R2 - R1) * m1}


def bloch_gruneisen_rho(T_C, theta=THETA_R, rrr=RRR, rho20=RHO_CU_20):
    """P2 lattice model, calibrated at 20 C. Integral by Gauss-Legendre (deterministic)."""
    def J5(y):
        x, w = np.polynomial.legendre.leggauss(64)
        xs = 0.5 * y * (x + 1.0)
        f = xs ** 5 / ((np.exp(xs) - 1.0) * (1.0 - np.exp(-xs)))
        return 0.5 * y * float(np.sum(w * f))
    T20 = 293.15
    g20 = (T20 / theta) ** 5 * J5(theta / T20)
    rho0 = rho20 / rrr
    A = (rho20 - rho0) / g20
    T = float(T_C) + 273.15
    return rho0 + A * (T / theta) ** 5 * J5(theta / T)


def green_kubo_rho(T_C, n_carriers, seed=1, parallel=True, window=200.0, steps_per_tau=10, rho_model=None):
    """P2: carrier ensemble -> resistivity at T_C via the Einstein relation (the integrated
    Green-Kubo form): sigma = n e^2 D / (k T), with <x(T)^2> = 2 D [T - tau (1 - e^{-T/tau})] for a
    stationary OU velocity (evaluated exactly for the discrete trapezoid estimator). Returns (rho_emulated, rho_target, stderr_rel)."""
    rho_t = bloch_gruneisen_rho(T_C) if rho_model is None else rho_model(T_C)
    tau = M_E / (N_CU * E_CHARGE ** 2 * rho_t)
    T = float(T_C) + 273.15
    vth = math.sqrt(K_B * T / M_E)
    dt = tau / steps_per_tau
    steps = int(window * steps_per_tau)
    out = np.zeros(int(n_carriers))
    _k("gk", parallel)(tau, vth, dt, steps, int(seed), out)
    msd = float(np.sum(out)) / out.size
    se = float(np.std(out)) / math.sqrt(out.size)
    # exact expectation of the discrete (trapezoid) estimator for the target tau, so that no
    # time-step bias is left: E[x^2] = vth^2 dt^2 sum_ij w_i w_j a^|i-j|
    a = math.exp(-dt / tau)
    w = np.ones(steps + 1); w[0] = w[-1] = 0.5
    c = np.correlate(w, w, "full")[steps:]               # lag 0..steps
    S = c[0] + 2.0 * float(np.sum(c[1:] * a ** np.arange(1, steps + 1)))
    msd_expected = vth ** 2 * dt ** 2 * S
    # D scales msd linearly at fixed dt/tau, so rho_emulated = rho_target * E[msd] / msd
    return rho_t * msd_expected / msd, rho_t, se / msd


def walk_dT(q_v, a, b, k, n_walkers, seed=2, parallel=True, dt_rel=2e-4):
    """P3: Feynman-Kac conduction rise at the cell edge. dt = dt_rel * (b - a)^2."""
    dt = dt_rel * (b - a) ** 2
    max_steps = int(50.0 * b * b / dt) + 10
    out = np.zeros(int(n_walkers))
    _k("walk", parallel)(float(a), float(b), dt, int(seed), max_steps, out)
    Et = float(np.sum(out)) / out.size
    se = float(np.std(out)) / math.sqrt(out.size)
    return q_v / (2.0 * k) * Et, q_v / (2.0 * k) * se


def axis_and_dsv_points(L, dsv, n_axis=21, n_dsv=400):
    from .fields import sphere_points
    z = np.linspace(-0.5 * L, 0.5 * L, n_axis)
    ax = np.column_stack([0 * z, 0 * z, z])
    return ax, sphere_points(dsv / 2.0, n_dsv)
