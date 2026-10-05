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


def channel_flow(v, D, length, T_bulk, correlation="dittus-boelter", K_minor=1.5, eta_pump=0.7,
                fluid=None, friction_multiplier=1.0):
    """Single circular channel: Re, Pr, h, f, dp, per-channel mass flow.

    fluid=None keeps the water correlations. A callable fluid(T) is only used
    by catalog coolants. The published default path does not pass it.

    friction_multiplier scales Darcy f for pressure drop only (MON stacked-channel
    rule: 10–20× smooth friction, conventional h). Nu / h always use the smooth
    correlation. Default multiplier 1.0 preserves the published continuum path.
    """
    w = water_props(T_bulk) if fluid is None else fluid(T_bulk)
    Re = w["rho"] * v * D / w["mu"]
    Pr = w["Pr"]
    if Re < 2300.0:
        Nu = 4.36                      # laminar, uniform heat flux
        Nu_g = 4.36
    else:
        Nu_g = nusselt_gnielinski(Re, Pr)
        Nu = nusselt_dittus_boelter(Re, Pr) if correlation == "dittus-boelter" else Nu_g
    h = Nu * w["k"] / D
    f_smooth = friction_factor(Re)
    mult = float(friction_multiplier)
    if mult <= 0.0:
        raise ValueError("friction_multiplier must be > 0")
    f = f_smooth * mult
    dp = (f * length / D + K_minor) * 0.5 * w["rho"] * v ** 2
    area = math.pi * D ** 2 / 4.0
    return {"Re": Re, "Pr": Pr, "Nu": Nu, "Nu_gnielinski": Nu_g, "h": h,
            "h_gnielinski": Nu_g * w["k"] / D, "f": f, "f_smooth": f_smooth,
            "friction_multiplier": mult, "dp": dp,
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


def lumped_time_constant(P20, alpha, C_th, R_th):
    """E19 is linear in T: dT/dt = -(T - T_ss)/tau, tau = C_th / (1/R_th - alpha P20).

    Returns inf when alpha*P20 >= 1/R_th (thermal runaway, no steady state)."""
    k = 1.0 / R_th - alpha * P20
    return float("inf") if k <= 0.0 else C_th / k


def lumped_steady_state(P20, alpha, R_th, T_w):
    """Steady state of E19 with cooling: P20 (1 + alpha (T-20)) = (T - T_w)/R_th."""
    k = 1.0 / R_th - alpha * P20
    if k <= 0.0:
        return float("inf")
    return (P20 * (1.0 - 20.0 * alpha) + T_w / R_th) / k


def adiabatic_time_to(T_limit, T0, rho20, alpha, J, dens, cp):
    """Local adiabatic heating of copper carrying current density J (no cooling, no conduction):
    dens cp dT/dt = rho20 (1 + alpha (T-20)) J^2  ->  closed-form time from T0 to T_limit.

    Applied at the inner radius (largest J for J ~ 1/r) this is the conservative
    pump-failure bound; the lumped E19 value averages J^2 over the whole stack."""
    if T_limit <= T0:
        return 0.0
    k = rho20 * alpha * J ** 2 / (dens * cp)
    return math.log((T_limit - 20.0 + 1.0 / alpha) / (T0 - 20.0 + 1.0 / alpha)) / k


def transient_summary(res, alpha, T_limit=85.0):
    """Post-process E19 for a design result dict (used by examples/run_all.py)."""
    from .constants import RHO_CU_20, DENS_CU, CP_CU
    P20, C_th, R_th = res["P20_W"], res["C_th_J_K"], res["R_th_K_W"]
    T_w = res["T_in"] + 0.5 * res["dT_water_mixed_K"]
    tau = lumped_time_constant(P20, alpha, C_th, R_th)
    Tss = lumped_steady_state(P20, alpha, R_th, T_w)
    # integrate well past the transient so the 63 % point is measured against the true steady state
    tt, Tc = transient_lumped(P20, alpha, C_th, R_th, T_w, t_end=10.0 * tau, n=20001)
    tau63 = float(np.interp(20.0 + (1.0 - math.exp(-1.0)) * (Tss - 20.0), Tc, tt))
    tf, Tf = transient_lumped(P20, alpha, C_th, R_th, T_w, t_end=6 * 3600.0, n=21601, cooling=False)
    t85 = float(tf[np.argmax(Tf >= T_limit)]) if np.any(Tf >= T_limit) else None
    t85_hot = adiabatic_time_to(T_limit, res["T_hot_C"], RHO_CU_20, alpha,
                                res["J_cu_inner_A_per_mm2"] * 1e6, DENS_CU, CP_CU)
    return {"transient_tau63_s": tau63, "transient_tau_analytic_s": tau, "transient_T_steady_C": Tss,
            "transient_T_end_60s_C": float(np.interp(60.0, tt, Tc)),
            "pump_failure_time_to_85C_s": t85, "pump_failure_hotspot_adiabatic_s": t85_hot}
