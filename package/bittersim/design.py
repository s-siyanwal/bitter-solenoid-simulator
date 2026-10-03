"""Coupled electro-thermal-hydraulic evaluation of a Bitter magnet design.

Electrical model (PDF section "Bitter Magnet Radial Current Distribution"):
  E29 J(r) = V0 / (2 pi rho r)   (V0 = voltage per plate/turn)
  E30 I_plate = V0 d ln(R2/R1) / (2 pi rho)
  E31 P = V0^2 L ln(R2/R1) / (2 pi rho)   (fully dense stack)
With insulation (thickness d_ins) and cooling holes (area fraction f_h) the
copper fill factor is lambda = d/(d+d_ins) * (1 - f_h) (ASSUMPTION: holes are
graded so they do not distort the 1/r profile).  Smeared density
J_avg = C/r, C = lambda V0/(2 pi rho);  NI = C L ln(R2/R1);
P = 2 pi rho C^2 L ln(R2/R1)/lambda;  R_turn = 2 pi rho/(d (1-f_h) ln(R2/R1)).
The current is set so that Bz(0) = B0_target exactly via E4.
"""
import math
import numpy as np
from .constants import MU_0, ALPHA_CU, K_CU, DENS_CU, CP_CU, RHO_CU_20
from .materials import rho_cu
from .fields import CoilLoops, homogeneity_ppm
from . import thermal

DEFAULTS = dict(
    R1=0.05,          # inner (bore) radius [m]   PDF: R1 >= 0.05 m
    R2=0.15,          # outer radius [m]          PDF code initial guess
    L=0.8,            # stack length [m]          PDF code initial guess
    d_plate=2.0e-3,   # copper plate thickness [m]            ASSUMPTION
    d_ins=0.25e-3,    # insulator thickness [m]               ASSUMPTION
    D_hole=5.0e-3,    # cooling hole diameter [m]             PDF code: 5 mm
    v_flow=2.5,       # water velocity [m/s]                  PDF code initial guess
    pitch_factor=2.0, # hole pitch / hole diameter            ASSUMPTION
    T_in=20.0,        # inlet water temperature [C]           PDF code
    B0=0.5,           # target isocentre field [T]            PDF
    dsv=0.03,         # DSV diameter for homogeneity [m]      ASSUMPTION
    K_minor=1.5,      # inlet+outlet minor-loss coefficient   ASSUMPTION
    eta_pump=0.7,     # pump efficiency                       ASSUMPTION
    correlation="dittus-boelter",
)


class BitterDesign(object):
    def __init__(self, **kw):
        p = dict(DEFAULTS)
        unknown = set(kw) - set(p)
        if unknown:
            raise ValueError("unknown parameters: %s" % sorted(unknown))
        p.update(kw)
        self.p = p
        for k, v in p.items():
            setattr(self, k, v)

    def hole_layout(self):
        pitch = self.pitch_factor * self.D_hole
        n_rows = max(1, int(math.floor((self.R2 - self.R1) / pitch + 1e-9)))
        dr = (self.R2 - self.R1) / n_rows
        r_rows = self.R1 + (np.arange(n_rows) + 0.5) * dr
        n_per = np.maximum(1, np.floor(2 * np.pi * r_rows / pitch + 1e-9)).astype(int)
        n_holes = int(n_per.sum())
        f_h = n_holes * math.pi * self.D_hole ** 2 / 4.0 / (math.pi * (self.R2 ** 2 - self.R1 ** 2))
        return {"pitch": pitch, "n_rows": n_rows, "dr_row": dr, "r_rows": r_rows,
                "n_per_row": n_per, "n_holes": n_holes, "hole_fraction": f_h}


def _audit_warnings(res, radial_rho):
    """Audit notes. They do not change the numeric continuum solution."""
    w = []
    if not radial_rho:
        w.append("Uniform rho(T) is used. There is no radial rho(T) feedback into J = C/r. "
                 "The optional correction is off so published numbers stay bit-stable. "
                 "At 0.5 T and about 26 C this error is small.")
    else:
        w.append("Radial rho(T) correction is reported only. P, V and the locked B0 are unchanged.")
    w.append("Cooling holes are modelled as graded round holes that do not disturb the 1/r profile. "
             "This is not a Florida-Bitter plate. Real plates use elongated, staggered holes "
             "because round holes concentrate current and hoop stress.")
    w.append("Helical slit and plate-to-plate contact resistance are omitted from this resistance "
             "(%.4g ohm). Contact and busbar resistance can dominate a coil near 1.76 mohm." % res["R_total_ohm"])
    w.append("Hoop stress bound %.4g MPa is the thin-ring upper estimate sigma = C*B/lambda. "
             "It is far below hard-copper yield. This 0.5 T design is not stress-limited." % res["sigma_hoop_max_MPa"])
    w.append("No critical heat flux, onset of nucleate boiling, or deionized-water chemistry "
             "in the continuum evaluation. A pump-failure time to 85 C is an adiabatic lumped "
             "estimate and is invalid once the wall boils.")
    w.append("Swiss-roll SNR gain is a heuristic. It is not an MRI SNR measurement. mu_eff(0) = 1.")
    near = []
    ppm = res.get("homogeneity_ppm")
    if ppm is not None and ppm >= 98.0:
        near.append("homogeneity")
    if res["R2"] >= 0.295:
        near.append("R2 upper bound")
    if res["R1"] <= 0.0505:
        near.append("R1 lower bound")
    if res["Re"] <= 1.05e4:
        near.append("Re >= 1e4")
    if res["d_plate"] >= 5.8e-3:
        near.append("plate thickness")
    if res["pitch_factor"] >= 7.7:
        near.append("hole pitch")
    if near:
        w.append("Constraints at or near their bounds: %s." % ", ".join(near))
    else:
        w.append("No tracked bound is active at this point. Interactive mode still lists each limit.")
    w.append("No FEniCS model is included. Off-axis checks stay on the loop field and the segment engine.")
    return w


def evaluate_design(design=None, homogeneity=True, n_iter=30, rho_scale=1.0,
                    radial_rho=False, fluid=None, k_solid=None, alpha=None, rho20=None,
                    dens=None, cp_solid=None, **kw):
    """Return a dict of all derived quantities for a design.

    Extra arguments default to the published path. rho_scale=1, radial_rho off,
    fluid None, and OFHC constants leave every continuum number unchanged.
    """
    d = design if design is not None else BitterDesign(**kw)

    def rho_at(T):
        if rho20 is None and alpha is None:
            val = float(rho_cu(T))
        else:
            r0 = RHO_CU_20 if rho20 is None else float(rho20)
            al = ALPHA_CU if alpha is None else float(alpha)
            val = r0 * (1.0 + al * (float(T) - 20.0))
        if rho_scale != 1.0:
            val = val * float(rho_scale)
        return val

    R1, R2, L = d.R1, d.R2, d.L
    lnr = math.log(R2 / R1)
    lay = d.hole_layout()
    f_h = lay["hole_fraction"]
    lam_ax = d.d_plate / (d.d_plate + d.d_ins)
    lam = lam_ax * (1.0 - f_h)
    n_turns = L / (d.d_plate + d.d_ins)

    # current amplitude for the target field (E4)
    C = d.B0 / (MU_0 * (math.asinh(L / (2 * R1)) - math.asinh(L / (2 * R2))))
    NI = C * L * lnr
    I = NI / n_turns

    # hydraulics (per hole; all holes in parallel at the same velocity)
    n_holes = lay["n_holes"]
    A_wet_cu_per_len = math.pi * d.D_hole * lam_ax          # copper part of hole wall
    T_cu = d.T_in + 10.0
    T_bulk = d.T_in + 5.0
    for _ in range(n_iter):
        rho = rho_at(T_cu)
        P = 2 * math.pi * rho * C ** 2 * L * lnr / lam
        fl = thermal.channel_flow(d.v_flow, d.D_hole, L, T_bulk, d.correlation, d.K_minor, d.eta_pump, fluid=fluid)
        m_tot = fl["m_dot"] * n_holes
        dT_mix = P / (m_tot * fl["cp"])
        dT_film_mean = P / (fl["h"] * A_wet_cu_per_len * n_holes * L)
        T_bulk_new = d.T_in + 0.5 * dT_mix
        T_cu_new = T_bulk_new + dT_film_mean
        if abs(T_cu_new - T_cu) < 1e-10 and abs(T_bulk_new - T_bulk) < 1e-10:
            T_cu, T_bulk = T_cu_new, T_bulk_new
            break
        T_cu, T_bulk = T_cu_new, T_bulk_new
    rho = rho_at(T_cu)
    P = 2 * math.pi * rho * C ** 2 * L * lnr / lam
    R_turn = 2 * math.pi * rho / (d.d_plate * (1 - f_h) * lnr)
    R_tot = R_turn * n_turns
    V = I * R_tot

    # hot-spot: innermost hole row, outlet end (z = L), local resistivity iterated
    n0 = lay["n_per_row"][0]
    r_out0 = R1 + lay["dr_row"]
    cell_area = lay["dr_row"] * 2 * math.pi * lay["r_rows"][0] / n0
    b = math.sqrt(cell_area / math.pi)
    a = d.D_hole / 2.0
    k_cond = K_CU if k_solid is None else k_solid
    T_hot = T_cu
    for _ in range(n_iter):
        rho_h = rho_at(T_hot)
        q_line = (2 * math.pi / n0) * rho_h * C ** 2 / lam * math.log(r_out0 / R1)   # W/m per hole
        dT_w = q_line * L / (fl["m_dot"] * fl["cp"])
        q_flux = q_line / A_wet_cu_per_len
        dT_film = q_flux / fl["h"]
        q_v_cu = rho_h * (C / (lam * R1)) ** 2
        dT_cond = thermal.annulus_conduction_dT(q_v_cu, a, b, k_cond)
        T_new = d.T_in + dT_w + dT_film + dT_cond
        if abs(T_new - T_hot) < 1e-10:
            T_hot = T_new
            break
        T_hot = T_new
    T_wall_hot = d.T_in + dT_w + dT_film

    # energy balance check: sum of cell heats over all rows == P (at uniform rho)
    edges = R1 + np.arange(lay["n_rows"] + 1) * lay["dr_row"]
    cell_heat = (2 * math.pi * rho * C ** 2 / lam * np.log(edges[1:] / edges[:-1])) * L
    energy_residual = (cell_heat.sum() - P) / P

    dp = fl["dp"]
    Q_tot = fl["Q"] * n_holes
    P_pump_hyd = dp * Q_tot
    P_pump = P_pump_hyd / d.eta_pump

    # mechanics (closed-form part); profile in mechanics.py
    sigma_hoop_est = C * d.B0 / lam   # E21 with Bz(R1) <= B0 (upper bound)

    # thermal transient parameters (E19). Default density and cp stay OFHC.
    # Catalog alloys pass their own values; omitting them leaves the published path unchanged.
    dens_s = DENS_CU if dens is None else float(dens)
    cp_s = CP_CU if cp_solid is None else float(cp_solid)
    V_cu = math.pi * (R2 ** 2 - R1 ** 2) * L * lam
    C_th = V_cu * dens_s * cp_s
    R_th = 1.0 / (fl["h"] * A_wet_cu_per_len * n_holes * L)
    P20 = 2 * math.pi * rho_at(20.0) * C ** 2 * L * lnr / lam
    adiabatic_rate = P / C_th

    res = dict(d.p)
    res.update({
        "C_A_per_m": C, "NI_At": NI, "n_turns": n_turns, "I_A": I, "V_total_V": V,
        "V_per_turn_V": V / n_turns, "R_total_ohm": R_tot, "R_turn_ohm": R_turn,
        "P_elec_W": P, "rho_cu_mean": rho, "T_cu_mean_C": T_cu,
        "J_cu_inner_A_per_mm2": C / (lam * R1) / 1e6, "J_cu_outer_A_per_mm2": C / (lam * R2) / 1e6,
        "fill_axial": lam_ax, "fill_total": lam, "hole_fraction": f_h,
        "n_holes": n_holes, "n_rows": lay["n_rows"], "holes_inner_row": int(n0),
        "Re": fl["Re"], "Pr": fl["Pr"], "Nu": fl["Nu"], "Nu_gnielinski": fl["Nu_gnielinski"],
        "h_W_m2K": fl["h"], "h_gnielinski_W_m2K": fl["h_gnielinski"], "f_darcy": fl["f"],
        "dp_Pa": dp, "flow_m3_s": Q_tot, "flow_L_min": Q_tot * 60000.0, "m_dot_kg_s": m_tot,
        "dT_water_mixed_K": dT_mix, "T_out_mixed_C": d.T_in + dT_mix,
        "dT_water_inner_K": dT_w, "dT_film_inner_K": dT_film, "dT_cond_inner_K": dT_cond,
        "T_wall_hot_C": T_wall_hot, "T_hot_C": T_hot, "q_flux_inner_W_m2": q_flux,
        "P_pump_hyd_W": P_pump_hyd, "P_pump_W": P_pump, "P_total_W": P + P_pump,
        "energy_residual": energy_residual, "sigma_hoop_max_MPa": sigma_hoop_est / 1e6,
        "C_th_J_K": C_th, "R_th_K_W": R_th, "P20_W": P20, "adiabatic_rate_K_s": adiabatic_rate,
        "f_larmor_MHz": 2.6752218744e8 * d.B0 / (2 * math.pi) / 1e6,
        "mass_cu_kg": V_cu * dens_s,
    })
    if homogeneity:
        coil = CoilLoops(R1, R2, L, C, "bitter", nr=12, nz=16, z_panels=4)
        ppm, B0num = homogeneity_ppm(coil, d.dsv / 2.0)
        res["homogeneity_ppm"] = ppm
        res["B0_numeric_T"] = B0num
    res["warnings"] = _audit_warnings(res, radial_rho)
    if radial_rho:
        contrast = ALPHA_CU * (res["T_hot_C"] - res["T_cu_mean_C"])
        res["radial_rho_dB_over_B"] = -0.5 * contrast
        res["radial_rho_model_grade"] = "approximate"
        res["radial_rho_note"] = (
            "Approximate fractional isocentre shift from a hotter inner radius. "
            "Not applied to P, V, or the locked B0.")
    return res
