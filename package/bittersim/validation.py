"""Verification & validation suite (PDF section "Verification, Validation,
and Test Methods").  Every function returns plain numbers so that the
results tables in VALIDATION.md are generated, not typed."""
import math
import time
try:
    import resource  # POSIX only; absent on Windows
except ImportError:
    resource = None
import numpy as np
from .constants import MU_0
from . import fields, biot_savart as bs, inductance as ind
from .elliptic import loop_field_scipy, loop_field_numba, ellipke_agm, field_from_loops
from scipy.special import ellipk, ellipe


def v_elliptic():
    m = np.linspace(0, 0.999, 500)
    e = 0.0
    for mi in m:
        K, E = ellipke_agm(mi)
        e = max(e, abs(K - ellipk(mi)) / ellipk(mi), abs(E - ellipe(mi)) / ellipe(mi))
    rho = np.linspace(0.001, 0.2, 40); z = np.linspace(-0.2, 0.2, 40)
    RR, ZZ = np.meshgrid(rho, z)
    Br_s, Bz_s = loop_field_scipy(0.07, 100.0, RR.ravel(), ZZ.ravel())
    Br_n, Bz_n = field_from_loops(RR.ravel().copy(), ZZ.ravel().copy(), np.array([0.07]), np.array([0.0]), np.array([100.0]))
    Bs = np.hypot(Br_s, Bz_s)
    dB = np.max(np.hypot(Br_s - Br_n, Bz_s - Bz_n) / Bs)
    return {"agm_vs_scipy_max_rel": float(e), "numba_vs_scipy_loop_max_rel": float(dB)}


def v_single_loop():
    a, I = 0.1, 1000.0
    z = np.linspace(-0.3, 0.3, 61)
    num = np.array([loop_field_numba(a, I, 0.0, zi)[1] for zi in z])
    # off-axis-limit check: evaluate at tiny rho via elliptic branch
    num2 = np.array([loop_field_numba(a, I, 1e-7, zi)[1] for zi in z])
    ana = fields.B_loop_axis(a, I, z)
    # segment engine: polygon loop at centre
    rows = []
    for n in (36, 72, 144, 288):
        s0, s1, c = bs.polygon_loop(a, 0.0, n, I)
        B = bs.biot_savart(np.zeros((1, 3)), s0, s1, c)[0, 2]
        rows.append((n, abs(B - MU_0 * I / (2 * a)) / (MU_0 * I / (2 * a))))
    return {"axis_branch_max_rel": float(np.max(np.abs(num - ana) / ana)),
            "elliptic_branch_near_axis_max_rel": float(np.max(np.abs(num2 - ana) / ana)),
            "segments_center_error": rows}


def v_long_solenoid():
    out = []
    R, n, I = 0.05, 1000.0, 10.0
    for ratio in (10, 100, 1000):
        L = ratio * R
        coil = fields.CoilLoops(R, R, L, n * I, "thin", nz=16, z_panels=max(4, ratio // 2))
        B = coil.field(np.array([0.0]), np.array([0.0]))[1][0]
        Bi = fields.B_ideal(n, I)
        Bt = float(fields.B_thin_axis(R, L, n, I, 0.0))
        out.append({"L_over_R": ratio, "B_loops": float(B), "B_ideal": Bi, "B_E2": Bt,
                    "rel_vs_ideal": abs(B - Bi) / Bi, "rel_vs_E2": abs(B - Bt) / Bt})
    return out


def v_thick_closed_form():
    R1, R2, L, j = 0.05, 0.15, 0.4, 2e7
    coil = fields.CoilLoops(R1, R2, L, j, "uniform", nr=12, nz=16, z_panels=4)
    z = np.linspace(-0.3, 0.3, 61)
    Bn = coil.field(np.zeros_like(z), z)[1]
    Ba = fields.B_thick_axis(R1, R2, L, j, z)
    c_pdf = fields.B_thick_center(R1, R2, L, j)
    C = 3e5
    coilb = fields.CoilLoops(R1, R2, L, C, "bitter", nr=12, nz=16, z_panels=4)
    Bnb = coilb.field(np.zeros_like(z), z)[1]
    Bab = fields.B_bitter_axis(R1, R2, L, C, z)
    return {"uniform_loops_vs_E3_max_rel": float(np.max(np.abs(Bn - Ba) / np.abs(Ba))),
            "E3_center_vs_PDF_form_rel": float(abs(Ba[30] - c_pdf) / c_pdf),
            "bitter_loops_vs_E4_max_rel": float(np.max(np.abs(Bnb - Bab) / np.abs(Bab)))}


def v_segment_convergence(R1=0.05, R2=0.15, L=0.4, j=2e7, n_turns=40, radial=8):
    """PDF test 1 (uniform thick solenoid) with planar rings, and Richardson (E7).
    Radial/axial filament discretisation error is removed by comparing against the
    loop model with the same filaments, isolating the polygon (segment-length h) error;
    the absolute error vs E3 (all discretisation errors) is also reported."""
    ref_exact = float(fields.B_thick_center(R1, R2, L, j))
    rows = []
    vals = []
    for n_ang in (45, 90, 180, 360, 720):
        s0, s1, c = bs.generate_coil_segments(R1, R2, L, j, "uniform", n_turns, radial, n_ang, helical=False)
        B = bs.biot_savart(np.zeros((1, 3)), s0, s1, c)[0, 2]
        vals.append(B)
        rows.append((n_ang, abs(B - ref_exact) / ref_exact, int(s0.shape[0])))
    # same-filament reference: limit n_ang -> inf equals sum of exact loops
    dz = L / n_turns; dr = (R2 - R1) / radial
    rc = R1 + (np.arange(radial) + 0.5) * dr
    zc = -L / 2 + (np.arange(n_turns) + 0.5) * dz
    RR, ZZ = np.meshgrid(rc, zc, indexing="ij")
    Iloop = (j * dr * dz) * np.ones(RR.size)
    Bfil = field_from_loops(np.zeros(1), np.zeros(1), RR.ravel().copy(), ZZ.ravel().copy(), Iloop)[1][0]
    poly_err = [(rows[k][0], abs(vals[k] - Bfil) / Bfil) for k in range(len(vals))]
    p, Bext = bs.richardson_order(vals[2], vals[3], vals[4])
    return {"B_exact_E3": ref_exact, "rows_vs_E3": rows, "rows_polygon_only": poly_err,
            "richardson_p": p, "B_richardson": Bext, "B_filament_limit": float(Bfil),
            "filament_limit_vs_E3_rel": abs(Bfil - ref_exact) / ref_exact,
            "richardson_vs_E3_rel": abs(Bext - ref_exact) / ref_exact}


def v_bitter_segments(R1=0.05, R2=0.15, L=0.4, C=3e5):
    """PDF test 2: Bitter 1/r filaments (helical, as in the PDF generator, and planar) vs E4."""
    ana = float(fields.B_bitter_center(R1, R2, L, C))
    out = {"B_E4": ana}
    for hel in (False, True):
        s0, s1, c = bs.generate_coil_segments(R1, R2, L, C, "bitter", 80, 24, 360, helical=hel)
        B = bs.biot_savart(np.zeros((1, 3)), s0, s1, c)[0]
        out["helical" if hel else "planar"] = {"Bz": float(B[2]), "rel_err": abs(B[2] - ana) / ana,
                                               "B_transverse": float(math.hypot(B[0], B[1])),
                                               "n_segments": int(s0.shape[0])}
    return out


def v_inductance():
    rows = []
    for R, l, N in ((0.05, 0.05, 50), (0.05, 0.5, 200), (0.05, 5.0, 1000)):
        Ln = ind.sheet_inductance(R, l, N)
        La = ind.nagaoka_inductance(R, l, N)
        Llong = MU_0 * N ** 2 * math.pi * R ** 2 / l
        rows.append({"R": R, "len": l, "N": N, "L_numeric_H": Ln, "L_nagaoka_H": La,
                     "rel": abs(Ln - La) / La, "L_long_H": Llong, "nagaoka_over_long": La / Llong})
    # thick coil in the thin limit
    thin = []
    for t in (5e-3, 5e-4, 5e-5):
        Lt = ind.coil_inductance(0.05, 0.05 + t, 0.5, 200, "uniform", nr=10)
        La = ind.nagaoka_inductance(0.05 + t / 2, 0.5, 200)
        thin.append((t, abs(Lt - La) / La))
    # energy check for a long solenoid: W = L I^2/2 vs B^2/(2 mu0) * pi R^2 l
    R, l, N, I = 0.05, 5.0, 1000, 10.0
    W1 = 0.5 * ind.sheet_inductance(R, l, N) * I ** 2
    B = MU_0 * N / l * I
    W2 = B ** 2 / (2 * MU_0) * math.pi * R ** 2 * l
    return {"sheet": rows, "thick_thin_limit": thin,
            "energy_long_W_LI2": W1, "energy_long_W_B2": W2, "energy_ratio": W1 / W2}


def v_homogeneity_convergence(res):
    out = []
    for nr, nz, pan in ((6, 8, 2), (12, 16, 4), (16, 24, 6), (24, 32, 8)):
        coil = fields.CoilLoops(res["R1"], res["R2"], res["L"], res["C_A_per_m"], "bitter", nr, nz, pan)
        ppm, B0 = fields.homogeneity_ppm(coil, res["dsv"] / 2)
        out.append({"nr": nr, "nz_total": nz * pan, "ppm": ppm, "B0": B0})
    return out


def v_performance(n_seg_target=20000, n_obs=20000):
    s0, s1, c = bs.generate_coil_segments(0.05, 0.15, 0.4, 3e5, "bitter", 50, 8, 50, helical=True)
    obs = np.random.RandomState(0).uniform(-0.03, 0.03, (n_obs, 3))
    bs.biot_savart(obs[:10], s0, s1, c)  # JIT warm-up
    t = time.time()
    bs.biot_savart(obs, s0, s1, c, chunk=5000)
    dt = time.time() - t
    rss = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
           if resource is not None else float("nan"))
    return {"n_segments": int(s0.shape[0]), "n_obs": n_obs, "seconds": dt,
            "pair_evals_per_s": s0.shape[0] * n_obs / dt, "peak_rss_MB": rss}


def v_beta_benchmark():
    """BETA (Bates et al., RSI 89 054704, 2018) Tables I–II and measured ΔP.

    Compares our continuum formulas to the published analytic/experimental
    anchors. Tolerances are loose where geometry differs (elongated Florida
    holes, Vinokur packing, overlap fill); the resistance Eq 3 check is tight.
    """
    from .materials import rho_cu, water_props
    from . import thermal
    from . import inductance as ind

    # --- Table I geometry / electrical ---
    R1, R2 = 0.020, 0.06986
    L = 0.0805
    t_plate = 0.5e-3
    lam = 0.8113
    N_eff = 77.625
    T_avg = 44.5
    I = 1175.0
    B_target = 1.03
    table_I = {
        "B_T": 1.03, "R_ohm": 17.679e-3, "V_V": 20.773, "P_W": 24.41e3,
        "L_H": 196.86e-6, "I_A": 1175.0, "N": 77.625, "Nc": 81,
    }
    rho = float(rho_cu(T_avg))
    lnr = math.log(R2 / R1)
    # BETA Eq 3: R = 2 π N ρ / (λ t ln(r2/r1))
    R_eq3 = 2.0 * math.pi * N_eff * rho / (lam * t_plate * lnr)
    V_eq3 = I * R_eq3
    P_eq3 = I ** 2 * R_eq3
    # Inductance of thick Bitter winding (E9); N = effective turns
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=Warning)
        L_num = float(ind.coil_inductance(R1, R2, L, N_eff, "bitter", nr=6))
    # Field amplitude C from E4 for cross-check (not a Table I column)
    C = B_target / (MU_0 * (math.asinh(L / (2 * R1)) - math.asinh(L / (2 * R2))))

    # --- Table II hydraulics (5 rings × 18 elongated holes) ---
    # Approximate hydraulic diameter from 2.5 mm radial × 4.5° span at mid-radius.
    r_mid = 0.5 * (R1 + R2)
    a_rad, a_arc = 2.5e-3, r_mid * math.radians(4.5)
    area_h = a_rad * a_arc          # rectangular proxy for elongated hole
    perim = 2.0 * (a_rad + a_arc)
    Dh = 4.0 * area_h / perim
    n_holes = 5 * 18
    m_tot = 2.16                     # kg/s
    T_in = 5.0
    w = water_props(T_in + 5.0)      # bulk ~ mid-rise
    A_one = area_h
    v_mean = (m_tot / n_holes) / (w["rho"] * A_one)
    # Velocity range in Table II: 1.18–1.96 m/s; h = 7191–8089
    table_II = {"h_min": 7191.0, "h_max": 8089.0, "T_wall_inner": 53.7, "T_wall_outer": 36.1,
                "v_min": 1.18, "v_max": 1.96, "m_dot": 2.16}
    h_rows = []
    for v in (table_II["v_min"], v_mean, table_II["v_max"]):
        fl = thermal.channel_flow(v, Dh, L, T_in + 5.0, friction_multiplier=1.0)
        h_rows.append({"v": v, "h": fl["h"], "Re": fl["Re"], "f": fl["f"], "dp": fl["dp"]})
    # Measured ΔP 7.72 kPa; BETA analytic 9.69 kPa with 2 mm effective roughness.
    # Our stack-friction rule (×15 on f) brackets the measured drop at mean velocity.
    fl_smooth = thermal.channel_flow(v_mean, Dh, L, T_in + 5.0, friction_multiplier=1.0)
    fl_stack = thermal.channel_flow(v_mean, Dh, L, T_in + 5.0, friction_multiplier=15.0)
    dp_meas = 7.72e3
    dp_beta_analytic = 9.69e3

    return {
        "table_I": table_I,
        "R_eq3_ohm": R_eq3,
        "V_eq3_V": V_eq3,
        "P_eq3_W": P_eq3,
        "L_numeric_H": L_num,
        "C_A_per_m": C,
        "rho_at_Tavg": rho,
        "rel_R": abs(R_eq3 - table_I["R_ohm"]) / table_I["R_ohm"],
        "rel_V": abs(V_eq3 - table_I["V_V"]) / table_I["V_V"],
        "rel_P": abs(P_eq3 - table_I["P_W"]) / table_I["P_W"],
        "rel_L": abs(L_num - table_I["L_H"]) / table_I["L_H"],
        "table_II": table_II,
        "Dh_m": Dh,
        "v_mean_m_s": v_mean,
        "h_at_velocities": h_rows,
        "dp_smooth_Pa": fl_smooth["dp"],
        "dp_stack15_Pa": fl_stack["dp"],
        "dp_measured_Pa": dp_meas,
        "dp_beta_analytic_Pa": dp_beta_analytic,
        "notes": (
            "R/V/P from BETA Eq 3 at Table I T_avg; L from E9 bitter inductance. "
            "h from Dittus–Boelter on a rectangular Dh proxy for elongated holes. "
            "ΔP: smooth vs friction_multiplier=15 vs BETA measured 7.72 kPa / analytic 9.69 kPa."
        ),
    }


def run_all(opt_res):
    return {"elliptic": v_elliptic(), "single_loop": v_single_loop(), "long_solenoid": v_long_solenoid(),
            "thick": v_thick_closed_form(), "segments": v_segment_convergence(),
            "bitter_segments": v_bitter_segments(), "inductance": v_inductance(),
            "homogeneity_convergence": v_homogeneity_convergence(opt_res), "performance": v_performance(),
            "beta": v_beta_benchmark()}
