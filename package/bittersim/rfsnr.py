"""Quasi-static RF receive model: surface loop, Swiss-roll slab and conducting tissue.

Replaces the E28 heuristic for SNR comparisons. Geometry (all ASSUMPTIONS, see DEFAULTS):
coil = circular loop of radius a in the plane z = 0, its axis normal to B0; a laterally infinite
slab of Swiss rolls (rolls along the coil axis) fills g < z < g + t; tissue (conductivity sigma)
fills z > g + t (with no slab the same space is air); the voxel is on the axis at depth d in tissue.

Q1  Magnetostatic scalar potential, Hankel components psi_k(z) J0(k rho). Free space: e^{+-kz}.
    Uniaxial slab (mu_z = mu, mu_t = 1): psi ~ e^{+-q z}, q = k / sqrt(mu); boundary conditions:
    psi continuous, mu_z d(psi)/dz continuous. `slab_transfer` solves the 4x4 system per k.
Q2  Loop field below the coil: Hz(rho, z) = (I a / 2) int k J1(k a) J0(k rho) tau_k e^{-k (z - z2)} dk.
Q3  Tissue loss (quasi-static, tissue reaction and permittivity neglected): E_phi = -j w A_phi,
    A_phi = mu0 (I a/2) int J1(ka) J1(k rho) tau_k e^{-k(z - z2)} dk, and Hankel-Parseval gives
    R_tissue = sigma w^2 mu0^2 pi a^2 / 4 * int |J1(ka) tau_k|^2 e^{-2 k s} / k^2 dk  (s = tissue offset).
    For a loop lying on a half-space this tends to sigma w^2 mu0^2 a^3 / 3 (int J1^2/x^2 = 4/(3 pi)).
Q4  Slab loss: P = (w mu0 mu''/2) int |Hz|^2 dV (only the axial permeability is lossy), so
    R_slab = w mu0 mu'' * 2 pi (a/2)^2 int dk/k |J1(ka)|^2 int_slab |h_k(z)|^2 dz.
Q5  Coil loss: round wire with skin effect, no proximity effect: R_coil = 2 a rho_cu / (d_w delta).
Q6  Reciprocity: SNR ~ |B1 per amp at voxel| / sqrt(T_c R_coil + T_s R_tissue + T_r R_slab).
    Results are ratios to the same coil with the slab replaced by air (same distances), and to the
    coil moved down onto the tissue (the no-metamaterial best case).
"""
import math
import numpy as np
from scipy.special import j1
from .constants import MU_0, RHO_CU_20, GAMMA_PROTON
from .swissroll import SwissRoll, larmor_hz

DEFAULTS = dict(
    a=0.03,            # coil radius [m]                                 ASSUMPTION
    d_wire=2e-3,       # coil wire diameter [m]                          ASSUMPTION
    gap=2e-3,          # coil-to-slab gap [m]                            ASSUMPTION
    t_slab=0.05,       # Swiss-roll slab thickness [m]                   ASSUMPTION
    depth=0.02,        # voxel depth below the tissue surface [m]        ASSUMPTION
    sigma=0.5,         # tissue conductivity at ~21 MHz [S/m]            ASSUMPTION (order of magnitude)
    T_coil=293.0, T_tissue=310.0, T_slab=293.0,                        # ASSUMPTION
)


def k_grid(a, n=6000):
    """Wavenumber grid with a log-spaced tail (same grid in docs/bittersim.js)."""
    u = np.linspace(0.0, 1.0, n)
    return (60.0 / a) * u ** 2 + 1e-9


def slab_transfer(k, mu, t, z1):
    """Q1: region-3 amplitude tau_k (incident psi = e^{-kz}, slab from z1 to z1+t) and the slab
    coefficients (A, Bp) with psi_slab = A e^{-q(z-z1)} + Bp e^{q(z - z1 - t)}."""
    k = np.asarray(k, dtype=complex)
    mu = complex(mu)
    q = k / np.sqrt(mu)
    E = np.exp(-q * t)
    inc = np.exp(-k * z1)
    eta = mu * q / k                                   # = sqrt(mu)
    # unknowns: r (reflected, as r e^{k(z-z1)}), A, Bp, tau
    # z1:  inc + r = A + Bp E ;  k(-inc + r) = eta k (-A + Bp E)
    # z2:  A E + Bp = tau     ;  eta k (-A E + Bp) = -k tau
    # eliminate r: 2 inc = A (1 + eta) + Bp E (1 - eta)
    # z2: Bp (1 + eta) = A E (eta - 1) + ... -> from tau eqs: Bp(eta + 1) = A E (eta - 1)
    Bp_over_A = E * (eta - 1.0) / (eta + 1.0)
    A = 2.0 * inc / ((1.0 + eta) + Bp_over_A * E * (1.0 - eta))
    Bp = Bp_over_A * A
    tau = A * E + Bp
    return tau, A, Bp, q


def _trap(y, x):
    """Trapezoid rule; keeps the imaginary part for complex integrands (B1 is a phasor)."""
    v = np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return complex(v) if np.iscomplexobj(y) else float(v)


def evaluate(mu=None, p=None, f=None, B0=0.5):
    """Q1-Q6 for one configuration. mu=None means no slab (air)."""
    q_ = dict(DEFAULTS); q_.update(p or {})
    a, g, t, d = q_["a"], q_["gap"], q_["t_slab"], q_["depth"]
    f = larmor_hz(B0) if f is None else f
    w = 2 * math.pi * f
    k = k_grid(a)
    J = j1(k * a)
    z2 = g + t
    if mu is None:
        tau = np.exp(-k * z2) + 0j
        R_slab = 0.0
    else:
        tau, A, Bp, q = slab_transfer(k, mu, t, g)
        # h_k(z) = -d psi/dz inside slab = q (A e^{-q s} - Bp e^{q (s - t)}), s in [0, t]
        s = np.linspace(0.0, t, 41)[:, None]
        h = q[None, :] * (A[None, :] * np.exp(-q[None, :] * s) - Bp[None, :] * np.exp(q[None, :] * (s - t)))
        hz2 = np.trapezoid(np.abs(h) ** 2, s[:, 0], axis=0) if hasattr(np, "trapezoid") else np.trapz(np.abs(h) ** 2, s[:, 0], axis=0)
        integrand = J ** 2 * hz2 / k
        R_slab = w * MU_0 * abs(complex(mu).imag) * 2 * math.pi * (a / 2.0) ** 2 * _trap(integrand, k)
    B1 = MU_0 * (a / 2.0) * _trap(k * J * tau * np.exp(-k * d), k)            # Q2 on axis, per amp
    R_tissue = q_["sigma"] * w ** 2 * MU_0 ** 2 * math.pi * a ** 2 / 4.0 * _trap(np.abs(J * tau) ** 2 / k ** 2, k)
    delta = math.sqrt(2 * RHO_CU_20 / (w * MU_0))
    R_coil = 2 * a * RHO_CU_20 / (q_["d_wire"] * delta)
    noise = math.sqrt(q_["T_coil"] * R_coil + q_["T_tissue"] * R_tissue + q_["T_slab"] * R_slab)
    return {"B1_T_per_A": abs(B1), "R_coil": R_coil, "R_tissue": R_tissue, "R_slab": R_slab,
            "snr_metric": abs(B1) / noise, "f_Hz": f}


def compare(B0=0.5, detune=None, loss_multiplier=50.0, p=None):
    """SNR of coil + Swiss-roll slab relative to (i) the same coil over an air gap of the same
    thickness and (ii) the coil placed directly on the tissue (gap and slab removed)."""
    q_ = dict(DEFAULTS); q_.update(p or {})
    fL = larmor_hz(B0)
    if detune is None:
        detune = best_detune(B0, loss_multiplier, q_)[0]
    sr = SwissRoll(fL * detune, loss_multiplier=loss_multiplier)
    mu = complex(sr.mu(fL))
    slab = evaluate(mu, q_, fL)
    air = evaluate(None, q_, fL)
    contact = evaluate(None, dict(q_, gap=0.0, t_slab=0.0), fL)
    return {"detune": float(detune), "mu": [mu.real, mu.imag], "slab": slab, "air": air, "contact": contact,
            "gain_vs_air": slab["snr_metric"] / air["snr_metric"],
            "gain_vs_contact": slab["snr_metric"] / contact["snr_metric"]}


def best_detune(B0=0.5, loss_multiplier=50.0, p=None, ratios=None):
    fL = larmor_hz(B0)
    ratios = np.linspace(1.0005, 1.10, 200) if ratios is None else ratios
    air = evaluate(None, p, fL)["snr_metric"]
    g = np.array([evaluate(complex(SwissRoll(fL * r, loss_multiplier=loss_multiplier).mu(fL)), p, fL)["snr_metric"] / air
                  for r in ratios])
    i = int(np.argmax(g))
    return float(ratios[i]), float(g[i]), ratios, g
