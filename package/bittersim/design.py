"""Coupled electro-thermal-hydraulic evaluation of a Bitter magnet design.

Electrical model (PDF section "Bitter Magnet Radial Current Distribution"):
  E29 J(r) = V0 / (2 pi rho r)   (V0 = voltage per plate/turn)
  E30 I_plate = V0 d ln(R2/R1) / (2 pi rho)
  E31 P = V0^2 L ln(R2/R1) / (2 pi rho)   (fully dense stack)
With insulation (thickness d_ins) and cooling holes (area fraction f_h) the
copper fill factor is lambda = d/(d+d_ins) * (1 - f_h).  Hole layouts:
  - "uniform": equal radial pitch, n_holes/row ~ r (constant areal density).
  - "montgomery": density ~ (R1/r)^2 — equal holes per ring, ring spacing ~ r
    (MON p.75–76).
  - "vinokur": Vinokur tanh stretching packing rings toward the bore (BETA Eq 4).
Smeared density J_avg = C/r, C = lambda V0/(2 pi rho);  NI = C L ln(R2/R1);
P_copper = 2 pi rho C^2 L ln(R2/R1)/lambda;  R_turn = 2 pi rho/(d (1-f_h) ln(R2/R1)).
Plate-to-plate contact resistance R_c (default 0) adds to R, V and P and deposits
local heat in the overlap sector.  The current is set so that Bz(0) = B0_target
exactly via E4.
"""
import math
import numpy as np
from .constants import (
    MU_0, ALPHA_CU, K_CU, DENS_CU, CP_CU, RHO_CU_20,
    R_C_OHM_DEFAULT, OVERLAP_DEG_DEFAULT, V_SUPPLY_MAX_DEFAULT,
    FRICTION_MULTIPLIER_SMOOTH,
)
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
    # --- P0 literature-audit parameters (defaults preserve published path) ---
    friction_multiplier=FRICTION_MULTIPLIER_SMOOTH,  # MON stack: use 10–20
    R_c_ohm=R_C_OHM_DEFAULT,                         # per interface [Ohm]
    overlap_deg=OVERLAP_DEG_DEFAULT,                 # overlap sector [deg]
    hole_layout_mode="uniform",                      # uniform|montgomery|vinokur
    vinokur_b=1.5,                                   # Vinokur stretch parameter
    elongated_aspect=1.0,                            # 1=round; >1 elongated ellipse
    V_supply_max=V_SUPPLY_MAX_DEFAULT,               # for R_c headroom report
)


def _ellipse_hole(D, aspect):
    """Area-preserving ellipse from round diameter D and aspect a/b >= 1.

    Returns area, hydraulic diameter, wet perimeter. aspect=1 is a circle.
    """
    aspect = float(aspect)
    if aspect < 1.0:
        raise ValueError("elongated_aspect must be >= 1")
    a = 0.5 * D * math.sqrt(aspect)
    b = 0.5 * D / math.sqrt(aspect)
    area = math.pi * a * b
    # Ramanujan approx for ellipse perimeter
    perim = math.pi * (3.0 * (a + b) - math.sqrt((3.0 * a + b) * (a + 3.0 * b)))
    Dh = 4.0 * area / perim
    return area, Dh, perim


def _vinokur_edges(R1, R2, n_rows, b):
    """Edges packed toward R1 via Vinokur (BETA Eq 4 form: 1 + tanh(b(ζ-1))/tanh b)."""
    b = float(b)
    if abs(b) < 1e-12:
        return R1 + (R2 - R1) * np.linspace(0.0, 1.0, n_rows + 1)
    tanh_b = math.tanh(b)
    zeta = np.linspace(0.0, 1.0, n_rows + 1)
    s = 1.0 + np.tanh(b * (zeta - 1.0)) / tanh_b   # 0 at ζ=0, 1 at ζ=1; packs at inner
    return R1 + (R2 - R1) * s


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
        """Return cooling-hole layout.

        Modes
        -----
        uniform : equal radial pitch; n_per_row ∝ r (constant areal density).
                  This maximises temperature non-uniformity relative to graded
                  layouts (AUDIT D1). Kept as default for published-path stability.
        montgomery : equal holes per ring, geometric ring spacing (∝ r), giving
                     areal density ∝ (R1/r)^2 (MON p.75–76).
        vinokur : Vinokur tanh stretching packing rings toward the bore
                  (BETA Eq 4); equal holes per ring.
        """
        pitch = self.pitch_factor * self.D_hole
        mode = str(self.hole_layout_mode).lower()
        area_h, Dh, perim = _ellipse_hole(self.D_hole, self.elongated_aspect)
        annulus = math.pi * (self.R2 ** 2 - self.R1 ** 2)

        if mode == "uniform":
            n_rows = max(1, int(math.floor((self.R2 - self.R1) / pitch + 1e-9)))
            edges = self.R1 + np.arange(n_rows + 1) * ((self.R2 - self.R1) / n_rows)
            r_rows = 0.5 * (edges[:-1] + edges[1:])
            n_per = np.maximum(1, np.floor(2 * math.pi * r_rows / pitch + 1e-9)).astype(int)
        elif mode == "montgomery":
            # Equal Δln r ⇒ spacing ∝ r; equal hole count per ring (MON p.75–76).
            # Choose n_rows so the innermost spacing is about `pitch`.
            n_rows = max(1, int(math.floor((self.R2 - self.R1) / pitch + 1e-9)))
            edges = self.R1 * (self.R2 / self.R1) ** (np.arange(n_rows + 1) / float(n_rows))
            r_rows = 0.5 * (edges[:-1] + edges[1:])
            n_each = max(1, int(math.floor(2 * math.pi * self.R1 / pitch + 1e-9)))
            n_per = np.full(n_rows, n_each, dtype=int)
        elif mode == "vinokur":
            n_rows = max(1, int(math.floor((self.R2 - self.R1) / pitch + 1e-9)))
            edges = _vinokur_edges(self.R1, self.R2, n_rows, self.vinokur_b)
            r_rows = 0.5 * (edges[:-1] + edges[1:])
            n_each = max(1, int(math.floor(2 * math.pi * self.R1 / pitch + 1e-9)))
            n_per = np.full(n_rows, n_each, dtype=int)
        else:
            raise ValueError("hole_layout_mode must be uniform|montgomery|vinokur, got %r" % mode)

        n_holes = int(n_per.sum())
        f_h = n_holes * area_h / annulus
        dr_rows = np.diff(edges)
        return {
            "pitch": pitch, "n_rows": int(n_rows), "dr_row": float(dr_rows[0]),
            "dr_rows": dr_rows, "edges": edges, "r_rows": r_rows,
            "n_per_row": n_per, "n_holes": n_holes, "hole_fraction": f_h,
            "mode": mode, "area_hole": area_h, "D_h": Dh, "perim_hole": perim,
            "elongated_aspect": float(self.elongated_aspect),
        }


def _audit_warnings(res, radial_rho):
    """Audit notes. They do not change the numeric continuum solution."""
    w = []
    if not radial_rho:
        w.append("Uniform rho(T) is used. There is no radial rho(T) feedback into J = C/r. "
                 "The optional correction is off so published numbers stay bit-stable. "
                 "At 0.5 T and about 26 C this error is small.")
    else:
        w.append("Radial rho(T) correction is reported only. P, V and the locked B0 are unchanged.")
    mode = res.get("hole_layout_mode", "uniform")
    if mode == "uniform":
        w.append("Cooling holes use a uniform-density round layout (equal pitch, n_per_row ~ r). "
                 "This is NOT a graded (MON/BETA) layout. Docs previously said 'graded'; that "
                 "mismatch is fixed. Use hole_layout_mode='montgomery' or 'vinokur' for graded rings. "
                 "Real high-power plates often use elongated, staggered Florida-Bitter holes.")
    else:
        w.append("Cooling holes use graded layout mode=%s (MON/BETA). "
                 "elongated_aspect=%.3g. Not a full Florida-Bitter plate model."
                 % (mode, res.get("elongated_aspect", 1.0)))
    if res.get("R_c_ohm", 0.0) <= 0.0:
        w.append("Plate-to-plate contact resistance R_c = 0 in this evaluation "
                 "(copper R = %.4g ohm). Set R_c_ohm > 0 to include joint resistance in R, V, P "
                 "and overlap-sector heat. V headroom allows max R_c ≈ %.4g ohm/interface."
                 % (res["R_copper_ohm"], res.get("R_c_max_per_interface_ohm", float("nan"))))
    else:
        w.append("Contact resistance R_c = %.4g ohm/interface over %d interfaces "
                 "(R_contact = %.4g ohm, %.1f %% of total R)."
                 % (res["R_c_ohm"], res["n_interfaces"], res["R_contact_ohm"],
                    100.0 * res["R_contact_ohm"] / res["R_total_ohm"]))
    if res.get("friction_multiplier", 1.0) <= 1.0 + 1e-12:
        w.append("Smooth-pipe friction (Petukhov). MON stacked channels run at 10–20×; "
                 "set friction_multiplier=15 (FRICTION_MULTIPLIER_STACK_MON) for that rule. "
                 "h is always the conventional smooth correlation.")
    else:
        w.append("Stacked-channel friction_multiplier = %.3g applied to Darcy f for Δp only; "
                 "h uses the conventional smooth Nu (MON p.104)." % res["friction_multiplier"])
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
    re_floor = res.get("Re_min_correlation", 5500.0)
    if res["Re"] <= 1.05 * re_floor:
        near.append("Re near soft floor (%.0f)" % re_floor)
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
    fluid None, friction_multiplier=1, R_c_ohm=0, hole_layout_mode='uniform',
    and OFHC constants leave every continuum number unchanged.
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
    n_interfaces = max(0, int(round(n_turns)) - 1)
    R_c = float(d.R_c_ohm)
    R_contact = n_interfaces * R_c
    Dh = lay["D_h"]
    perim_hole = lay["perim_hole"]
    # copper wet perimeter fraction of full hole perimeter (axial copper fill)
    A_wet_cu_per_len = perim_hole * lam_ax

    # current amplitude for the target field (E4)
    C = d.B0 / (MU_0 * (math.asinh(L / (2 * R1)) - math.asinh(L / (2 * R2))))
    NI = C * L * lnr
    I = NI / n_turns

    # hydraulics (per hole; all holes in parallel at the same velocity)
    n_holes = lay["n_holes"]
    T_cu = d.T_in + 10.0
    T_bulk = d.T_in + 5.0
    for _ in range(n_iter):
        rho = rho_at(T_cu)
        P_cu = 2 * math.pi * rho * C ** 2 * L * lnr / lam
        fl = thermal.channel_flow(
            d.v_flow, Dh, L, T_bulk, d.correlation, d.K_minor, d.eta_pump,
            fluid=fluid, friction_multiplier=d.friction_multiplier)
        m_tot = fl["m_dot"] * n_holes
        # Contact heat is deposited in overlap sectors, not uniformly in channels.
        # Channel mixed rise uses copper Joule heat only; contact added to hot-spot later.
        dT_mix = P_cu / (m_tot * fl["cp"])
        dT_film_mean = P_cu / (fl["h"] * A_wet_cu_per_len * n_holes * L)
        T_bulk_new = d.T_in + 0.5 * dT_mix
        T_cu_new = T_bulk_new + dT_film_mean
        if abs(T_cu_new - T_cu) < 1e-10 and abs(T_bulk_new - T_bulk) < 1e-10:
            T_cu, T_bulk = T_cu_new, T_bulk_new
            break
        T_cu, T_bulk = T_cu_new, T_bulk_new
    rho = rho_at(T_cu)
    P_cu = 2 * math.pi * rho * C ** 2 * L * lnr / lam
    R_turn = 2 * math.pi * rho / (d.d_plate * (1 - f_h) * lnr)
    R_copper = R_turn * n_turns
    R_tot = R_copper + R_contact
    V = I * R_tot
    P_contact = I ** 2 * R_contact
    P = P_cu + P_contact  # total electrical power

    # Per-ring hot spot (outlet end). Graded layouts need this; uniform reduces to inner ring.
    edges = lay["edges"]
    n_per = lay["n_per_row"]
    T_hot_rows = []
    dT_w_rows = []
    dT_film_rows = []
    dT_cond_rows = []
    q_flux_rows = []
    for i in range(lay["n_rows"]):
        n_i = int(n_per[i])
        r_in_i, r_out_i = float(edges[i]), float(edges[i + 1])
        r_mid = float(lay["r_rows"][i])
        cell_area = (r_out_i - r_in_i) * 2 * math.pi * r_mid / n_i
        b_cell = math.sqrt(cell_area / math.pi)
        a_hole = math.sqrt(lay["area_hole"] / math.pi)  # equivalent radius for conduction cell
        k_cond = K_CU if k_solid is None else k_solid
        T_h = T_cu
        for _ in range(n_iter):
            rho_h = rho_at(T_h)
            # heat per hole in this ring [W/m]
            q_line = (2 * math.pi / n_i) * rho_h * C ** 2 / lam * math.log(r_out_i / r_in_i)
            dT_w = q_line * L / (fl["m_dot"] * fl["cp"])
            q_flux = q_line / A_wet_cu_per_len
            dT_film = q_flux / fl["h"]
            q_v_cu = rho_h * (C / (lam * r_in_i)) ** 2
            dT_cond = thermal.annulus_conduction_dT(q_v_cu, a_hole, b_cell, k_cond)
            T_new = d.T_in + dT_w + dT_film + dT_cond
            if abs(T_new - T_h) < 1e-10:
                T_h = T_new
                break
            T_h = T_new
        T_hot_rows.append(T_h)
        dT_w_rows.append(dT_w)
        dT_film_rows.append(dT_film)
        dT_cond_rows.append(dT_cond)
        q_flux_rows.append(q_flux)

    # Overlap-sector contact heat: local rise over copper thermal mass of overlap wedge
    # (approximate; reported separately and folded into T_hot if R_c > 0).
    dT_contact = 0.0
    if P_contact > 0.0 and n_interfaces > 0:
        overlap_frac = float(d.overlap_deg) / 360.0
        # Steady: contact heat removed by neighbouring channel film over overlap arc length ~ overlap
        # Use a simple film balance on the overlap copper perimeter proxy.
        A_ov = A_wet_cu_per_len * n_holes * L * max(overlap_frac, 1e-6)
        dT_contact = (P_contact / n_interfaces) / (fl["h"] * A_ov / max(n_interfaces, 1))
        # milder: total contact heat / (h * overlap wet area)
        dT_contact = P_contact / (fl["h"] * A_wet_cu_per_len * n_holes * L * max(overlap_frac, 1e-6))

    i_hot = int(np.argmax(T_hot_rows))
    T_hot = float(T_hot_rows[i_hot]) + dT_contact
    dT_w = float(dT_w_rows[i_hot])
    dT_film = float(dT_film_rows[i_hot])
    dT_cond = float(dT_cond_rows[i_hot])
    q_flux = float(q_flux_rows[i_hot])
    T_wall_hot = d.T_in + dT_w + dT_film + dT_contact
    n0 = int(n_per[i_hot])

    # energy balance check: sum of cell heats over all rows == P_cu (at uniform rho)
    cell_heat = (2 * math.pi * rho * C ** 2 / lam * np.log(edges[1:] / edges[:-1])) * L
    energy_residual = (cell_heat.sum() - P_cu) / P_cu

    dp = fl["dp"]
    Q_tot = fl["Q"] * n_holes
    P_pump_hyd = dp * Q_tot
    P_pump = P_pump_hyd / d.eta_pump

    # mechanics (closed-form part); profile in mechanics.py
    sigma_hoop_est = C * d.B0 / lam   # E21 with Bz(R1) <= B0 (upper bound)

    dens_s = DENS_CU if dens is None else float(dens)
    cp_s = CP_CU if cp_solid is None else float(cp_solid)
    V_cu = math.pi * (R2 ** 2 - R1 ** 2) * L * lam
    C_th = V_cu * dens_s * cp_s
    R_th = 1.0 / (fl["h"] * A_wet_cu_per_len * n_holes * L)
    P20 = 2 * math.pi * rho_at(20.0) * C ** 2 * L * lnr / lam
    adiabatic_rate = P_cu / C_th

    # V headroom → max allowable R_c per interface at the supply ceiling
    V_cu_only = I * R_copper
    V_headroom = float(d.V_supply_max) - V_cu_only
    R_c_max = (V_headroom / I / n_interfaces) if (n_interfaces > 0 and I > 0) else float("inf")

    res = dict(d.p)
    res.update({
        "C_A_per_m": C, "NI_At": NI, "n_turns": n_turns, "I_A": I, "V_total_V": V,
        "V_per_turn_V": V / n_turns, "R_total_ohm": R_tot, "R_turn_ohm": R_turn,
        "R_copper_ohm": R_copper, "R_contact_ohm": R_contact, "n_interfaces": n_interfaces,
        "R_c_max_per_interface_ohm": R_c_max, "V_copper_V": V_cu_only,
        "V_headroom_V": V_headroom, "P_contact_W": P_contact, "P_copper_W": P_cu,
        "dT_contact_K": dT_contact, "overlap_deg": float(d.overlap_deg),
        "P_elec_W": P, "rho_cu_mean": rho, "T_cu_mean_C": T_cu,
        "J_cu_inner_A_per_mm2": C / (lam * R1) / 1e6, "J_cu_outer_A_per_mm2": C / (lam * R2) / 1e6,
        "fill_axial": lam_ax, "fill_total": lam, "hole_fraction": f_h,
        "n_holes": n_holes, "n_rows": lay["n_rows"], "holes_inner_row": int(n_per[0]),
        "holes_hot_row": n0, "hot_row_index": i_hot,
        "T_hot_rows_C": np.asarray(T_hot_rows, dtype=float),
        "hole_layout_mode": lay["mode"], "D_h_m": Dh, "elongated_aspect": lay["elongated_aspect"],
        "Re": fl["Re"], "Pr": fl["Pr"], "Nu": fl["Nu"], "Nu_gnielinski": fl["Nu_gnielinski"],
        "h_W_m2K": fl["h"], "h_gnielinski_W_m2K": fl["h_gnielinski"], "f_darcy": fl["f"],
        "f_darcy_smooth": fl["f_smooth"], "friction_multiplier": fl["friction_multiplier"],
        "Re_min_correlation": 5500.0,
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
