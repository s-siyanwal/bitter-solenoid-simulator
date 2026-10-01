"""Swiss-roll RF metamaterial (PDF section "Metamaterials in fMRI").

Acts ONLY at RF.  At omega = 0 the model gives mu_eff = 1 exactly, so the
Swiss rolls are modelled as having no effect on the static B0 field or its
homogeneity.  (Not modelled: the weak diamagnetism of copper,
chi ~ -1e-5, which in a real build could add ppm-level B0 distortion.)

E23 Larmor frequency: f0 = gamma B0 / (2 pi)
E24 Faraday EMF: E_ext = -j omega mu0 N S H_ext
E25 Lorentzian effective permeability (PDF):
    mu_eff(w) = 1 - F w^2 / (w^2 - w0^2 + j w Gamma)
E26 Pendry et al. (1999) geometry mapping (roll radius r, N turns, gap d,
    dielectric eps_r, sheet resistance sigma_s, lattice constant a):
    F = pi r^2 / a^2,  w0^2 = d c0^2 / (2 pi^2 eps_r r^3 (N-1)),
    Gamma = 2 sigma_s / (mu0 r (N-1))      [x loss_multiplier, ASSUMPTION]
E27 Sheet resistance with skin effect: sigma_s = rho / (delta (1 - exp(-t/delta))),
    delta = sqrt(2 rho / (omega mu0))
E28 SNR heuristic (ASSUMPTION - the PDF gives no quantitative SNR model):
    reciprocity coupling of a surface coil (radius a_c) to a source at depth l:
    eta0 = (a_c^2/(a_c^2 + l^2))^{3/2}; with a guide of axial permeability mu':
    eta_mm = 1 / (1 + (1/eta0 - 1)/mu')  (magnetic-circuit: guide reluctance
    scales as 1/mu'; eta_mm = 0 for mu' <= 0); added noise xi = kappa mu'';
    G_SNR = (eta_mm / eta0) / sqrt(1 + xi).
"""
import math
import numpy as np
from .constants import MU_0, C_LIGHT, GAMMA_PROTON, RHO_CU_20


def larmor_hz(B0):
    return GAMMA_PROTON * B0 / (2.0 * math.pi)


def faraday_emf(omega, N, S, H):
    return -1j * omega * MU_0 * N * S * H


def mu_eff(omega, F, omega0, Gamma):
    omega = np.asarray(omega, dtype=float)
    return 1.0 - F * omega ** 2 / (omega ** 2 - omega0 ** 2 + 1j * omega * Gamma)


def sheet_resistance(omega, t_foil, rho=RHO_CU_20):
    delta = math.sqrt(2.0 * rho / (omega * MU_0))
    return rho / (delta * (1.0 - math.exp(-t_foil / delta))), delta


class SwissRoll(object):
    """Default geometry (ASSUMPTIONS): r = 5 mm, a = 12 mm lattice, N = 28 turns,
    10 um Cu foil, eps_r = 3 dielectric; the gap d is solved so that w0 equals
    the requested tuning frequency."""

    def __init__(self, f_tune, r=5e-3, a=12e-3, N=28, eps_r=3.0, t_foil=10e-6,
                 loss_multiplier=50.0):
        self.r, self.a, self.N, self.eps_r, self.t_foil = r, a, N, eps_r, t_foil
        self.F = math.pi * r ** 2 / a ** 2
        self.omega0 = 2.0 * math.pi * f_tune
        self.gap = self.omega0 ** 2 * 2.0 * math.pi ** 2 * eps_r * r ** 3 * (N - 1) / C_LIGHT ** 2
        self.sigma_s, self.skin_depth = sheet_resistance(self.omega0, t_foil)
        self.loss_multiplier = loss_multiplier
        self.Gamma = loss_multiplier * 2.0 * self.sigma_s / (MU_0 * r * (N - 1))
        self.Q = self.omega0 / self.Gamma

    def mu(self, f):
        return mu_eff(2.0 * math.pi * np.asarray(f, dtype=float), self.F, self.omega0, self.Gamma)


def snr_gain(mu, depth=0.05, coil_radius=0.03, kappa=0.05):
    mu = np.asarray(mu, dtype=complex)
    eta0 = (coil_radius ** 2 / (coil_radius ** 2 + depth ** 2)) ** 1.5
    mr = mu.real
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = np.where(mr > 0, 1.0 / (1.0 + (1.0 / eta0 - 1.0) / np.where(mr > 0, mr, 1.0)), 0.0)
    xi = kappa * np.abs(mu.imag)
    return (eta / eta0) / np.sqrt(1.0 + xi)


def best_detuning(B0=0.5, **kw):
    """Scan the roll tuning frequency around the Larmor frequency and return
    the tuning that maximises the heuristic SNR gain at f_Larmor."""
    fL = larmor_hz(B0)
    ratios = np.linspace(0.95, 1.10, 3001)
    gains = np.array([float(snr_gain(SwissRoll(fL * q, **kw).mu(fL))) for q in ratios])
    k = int(np.argmax(gains))
    return ratios, gains, ratios[k], gains[k]
