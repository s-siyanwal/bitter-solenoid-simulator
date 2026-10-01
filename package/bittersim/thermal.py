"""Hydraulics and convective heat transfer (PDF section "Thermodynamics ...").

E11 Re = rho_f v D_h / mu_f ;  Pr = cp mu / k ;  h = Nu k / D_h
E12 Dittus-Boelter (PDF, heating): Nu = 0.023 Re^0.8 Pr^0.4
E13 Gnielinski (cross-check):      Nu = (f/8)(Re-1000)Pr / (1 + 12.7 sqrt(f/8)(Pr^(2/3)-1))
E14 Darcy friction: laminar 64/Re; turbulent Petukhov f = (0.790 ln Re - 1.64)^-2
E15 Pressure drop: dp = (f L/D + K_minor) rho v^2 / 2 ; pump power = dp Q / eta_pump
E16 Water energy balance per channel: m_dot cp dT_f/dz = q'  ->  T_f(L) = T_in + q' L/(m_dot cp)
E17 Newton cooling: T_s = T_f + q''/h
E18 Conduction in the copper cell around a hole (annulus a..b, insulated at b,
    uniform source q_v):  dT = q_v/(4k) [2 b^2 ln(b/a) - (b^2 - a^2)]
E19 Lumped transient: C_th dT/dt = P(T) - (T - T_w)/R_th,  P(T) = P20 (1 + alpha (T-20))
"""
import math
import numpy as np
from .materials import water_props


def friction_factor(Re):
    if Re < 2300.0:
        return 64.0 / Re
    return (0.790 * math.log(Re) - 1.64) ** -2


def nusselt_dittus_boelter(Re, Pr):
    return 0.023 * Re ** 0.8 * Pr ** 0.4


def nusselt_gnielinski(Re, Pr):
    f = friction_factor(max(Re, 2300.0))
    return (f / 8.0) * (Re - 1000.0) * Pr / (1.0 + 12.7 * math.sqrt(f / 8.0) * (Pr ** (2.0 / 3.0) - 1.0))


def channel_flow(v, D, length, T_bulk, correlation="dittus-boelter", K_minor=1.5, eta_pump=0.7):
    """Single circular channel: Re, Pr, h, f, dp, per-channel mass flow."""
    w = water_props(T_bulk)
    Re = w["rho"] * v * D / w["mu"]
    Pr = w["Pr"]
    if Re < 2300.0:
        Nu = 4.36                      # laminar, uniform heat flux
        Nu_g = 4.36
    else:
        Nu_g = nusselt_gnielinski(Re, Pr)
        Nu = nusselt_dittus_boelter(Re, Pr) if correlation == "dittus-boelter" else Nu_g
    h = Nu * w["k"] / D
    f = friction_factor(Re)
    dp = (f * length / D + K_minor) * 0.5 * w["rho"] * v ** 2
    area = math.pi * D ** 2 / 4.0
    return {"Re": Re, "Pr": Pr, "Nu": Nu, "Nu_gnielinski": Nu_g, "h": h,
            "h_gnielinski": Nu_g * w["k"] / D, "f": f, "dp": dp,
            "m_dot": w["rho"] * v * area, "Q": v * area, "cp": w["cp"], "rho": w["rho"],
            "eta_pump": eta_pump}


def annulus_conduction_dT(q_v, a, b, k):
    """E18."""
    if b <= a:
        return 0.0
    return q_v / (4.0 * k) * (2.0 * b ** 2 * math.log(b / a) - (b ** 2 - a ** 2))


def axial_profile(T_in, q_line, m_dot, cp, h, perim_cu, length, dT_cond=0.0, n=101):
    """E16/E17 along one channel: returns z, T_fluid(z), T_surface(z)."""
    z = np.linspace(0.0, length, n)
    Tf = T_in + q_line * z / (m_dot * cp)
    Ts = Tf + q_line / (h * perim_cu) + dT_cond
    return z, Tf, Ts


def transient_lumped(P20, alpha, C_th, R_th, T_w, T0=20.0, t_end=600.0, n=601, cooling=True):
    """E19 integrated with scipy.integrate.odeint (available in SciPy 1.1)."""
    from scipy.integrate import odeint
    def rhs(T, t):
        P = P20 * (1.0 + alpha * (T[0] - 20.0))
        q_out = (T[0] - T_w) / R_th if cooling else 0.0
        return [(P - q_out) / C_th]
    t = np.linspace(0.0, t_end, n)
    T = odeint(rhs, [T0], t)[:, 0]
    return t, T
