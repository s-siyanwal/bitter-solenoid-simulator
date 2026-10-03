"""Mesoscopic emulation of one copper RVE plus engineering scatter.

This is not molecular dynamics, not DFT, and not a particle model of the
2575 kg magnet. Every returned object carries model_grade. Heuristics are
labeled. The continuum evaluate_design path stays the reference and is not
replaced.
"""
import math
import numpy as np

from .constants import MU_0, ALPHA_CU
from .design import BitterDesign, evaluate_design
from . import catalog
from . import thermal

# CODATA-scale constants. Kept identical in docs/bittersim.js.
E_CHARGE = 1.60217662e-19       # C
M_ELECTRON = 9.1093837e-31      # kg
K_BOLTZMANN = 1.38064852e-23    # J/K

LABEL = "mesoscopic emulation, not molecular dynamics"

MODEL_GRADE = {
    "drude_tau": "established",
    "drift_velocity": "established",
    "johnson": "established",
    "shot_comparison": "established",
    "kohler": "approximate",
    "langevin_rve": "approximate",
    "contact_resistance": "approximate",
    "manufacturing_scatter": "approximate",
    "onb_bergles_rohsenow": "approximate",
    "antoine_tsat": "established",
    "hooge_1f": "speculative",
    "swissroll_snr": "speculative",
    "housing_convection": "approximate",
    "radial_rho": "approximate",
}


def saturation_temperature_C(p_pa):
    """Antoine inversion for water. P in Pa, T in C.

    log10(P_bar) = 5.1962 - 1730.63 / (T_C + 233.426), about 1-100 C,
    used here across the catalog window 0.07-0.2 MPa.
    """
    p_bar = float(p_pa) / 1.0e5
    if p_bar <= 0.0:
        raise ValueError("site pressure must be positive")
    return 1730.63 / (5.1962 - math.log10(p_bar)) - 233.426


def q_onb_bergles_rohsenow(p_pa, wall_minus_sat_C):
    """Bergles-Rohsenow (1964) water ONB heat flux [W/m^2].

    q = 1082 * p_bar^1.156 * [1.8 * dT_C] ^ (2.16 / p_bar^0.0234).
    Labeled approximate. Not a critical-heat-flux map.
    """
    p_bar = max(float(p_pa) / 1.0e5, 0.2)
    dts = float(wall_minus_sat_C)
    if dts <= 0.0:
        return 0.0
    return 1082.0 * (p_bar ** 1.156) * ((1.8 * dts) ** (2.16 / (p_bar ** 0.0234)))


def drude_rve(J, rho, T_C, B, n, d_plate, v_fermi, kohler=True, kohler_a=1.0):
    """Drude state at one current density. Kohler is optional and approximate.

    When kohler is off, J*E equals rho*J^2 exactly (the continuum cell heat).
    """
    rho = float(rho)
    J = float(J)
    n = float(n)
    if n <= 0.0 or rho <= 0.0:
        raise ValueError("n and rho must be positive")
    tau = M_ELECTRON / (n * E_CHARGE ** 2 * rho)
    # identity check rearranged: rho_from_tau must match the input rho
    rho_from_tau = M_ELECTRON / (n * E_CHARGE ** 2 * tau)
    v_d = J / (n * E_CHARGE)
    mean_free_path = float(v_fermi) * tau
    T_K = float(T_C) + 273.15
    if T_K < 0.0:
        T_K = 0.0
    v_thermal = math.sqrt(8.0 * K_BOLTZMANN * T_K / (math.pi * M_ELECTRON)) if T_K > 0.0 else 0.0
    t_plate = (float(d_plate) / v_d) if v_d != 0.0 else float("inf")
    t_mfp = tau  # lambda / v_F
    omega_c = E_CHARGE * float(B) / M_ELECTRON
    hall = omega_c * tau
    kohler_frac = (float(kohler_a) * hall ** 2) if kohler else 0.0
    rho_used = rho * (1.0 + kohler_frac)
    E_field = J * rho_used
    J_dot_E = J * E_field
    q_continuum = rho * J * J
    return {
        "label": LABEL,
        "n_m3": n,
        "tau_s": tau,
        "rho_ohm_m": rho,
        "rho_from_tau": rho_from_tau,
        "v_fermi_m_s": float(v_fermi),
        "v_thermal_m_s": v_thermal,
        "v_d_m_s": v_d,
        "mean_free_path_m": mean_free_path,
        "transit_plate_s": t_plate,
        "transit_mfp_s": t_mfp,
        "omega_c_tau": hall,
        "kohler_on": bool(kohler),
        "kohler_a": float(kohler_a),
        "kohler_drho_over_rho": kohler_frac,
        "model_grade_kohler": "approximate" if kohler else "off",
        "E_V_m": E_field,
        "J_dot_E_W_m3": J_dot_E,
        "q_continuum_W_m3": q_continuum,
        "power_density_rel_error": (J_dot_E - q_continuum) / q_continuum if q_continuum else 0.0,
        "model_grade": "established" if not kohler else "approximate",
    }


def langevin_ensemble(v_d, tau, T_C, n_carriers=4000, steps=600, seed=1):
    """1-D Langevin trajectory of notional carriers inside the RVE only.

    dv = -(v/tau) dt + (e/m) E dt + sqrt(2 kT /(m tau)) dW,
    with (e/m) E = v_d / tau from the Drude steady state.
    A fixed seed is required. This is a picture of the Drude process, not a
    sample of the magnet. At room temperature the thermal speed dwarfs v_d,
    so the standard error of a few thousand carriers does not resolve drift.
    """
    if seed is None:
        raise ValueError("seed is required")
    tau = float(tau)
    v_d = float(v_d)
    if tau <= 0.0:
        raise ValueError("tau must be positive")
    T_K = max(0.0, float(T_C) + 273.15)
    dt = tau / 20.0
    rs = np.random.RandomState(int(seed))
    sigma_v = math.sqrt(K_BOLTZMANN * T_K / M_ELECTRON) if T_K > 0.0 else 0.0
    v = rs.normal(0.0, sigma_v, size=int(n_carriers))
    noise = math.sqrt(2.0 * K_BOLTZMANN * T_K / (M_ELECTRON * tau) * dt) if T_K > 0.0 else 0.0
    for _ in range(int(steps)):
        kick = noise * rs.normal(size=v.shape[0]) if noise else 0.0
        v = v + (-v / tau) * dt + (v_d / tau) * dt + kick
    mean = float(np.mean(v))
    stderr = float(np.std(v, ddof=1) / math.sqrt(v.shape[0]))
    return {
        "label": LABEL,
        "model_grade": "approximate",
        "note": "Visualization of the Drude process in one RVE. Not a sample of the magnet.",
        "n_carriers": int(n_carriers),
        "steps": int(steps),
        "seed": int(seed),
        "mean_v_m_s": mean,
        "stderr_m_s": stderr,
        "v_d_m_s": v_d,
        "drift_resolved": bool(stderr < abs(v_d)) if v_d else False,
    }


def johnson_voltage(T_C, R_ohm, f_lo, f_hi):
    """Open-circuit Johnson-Nyquist voltage. S_v = 4 kT R, white over the band.

    The default band is a supply-sense bandwidth. A value at the Larmor
    frequency is a spectral density only. It is not the MRI noise floor.
    """
    if f_hi <= f_lo:
        raise ValueError("bandwidth upper bound must exceed the lower bound")
    T_K = float(T_C) + 273.15
    df = float(f_hi) - float(f_lo)
    S_v = 4.0 * K_BOLTZMANN * T_K * float(R_ohm)
    return {
        "model_grade": "established",
        "S_v_V2_per_Hz": S_v,
        "V_rms_V": math.sqrt(S_v * df),
        "f_lo_Hz": float(f_lo),
        "f_hi_Hz": float(f_hi),
        "note": "Supply-sense Johnson voltage. Coil Johnson noise is not the MRI noise floor.",
    }


def shot_vs_johnson(I_A, T_C, R_ohm, f_lo, f_hi):
    """Shot noise 2 e I is negligible next to Johnson current noise in a metal."""
    df = float(f_hi) - float(f_lo)
    T_K = float(T_C) + 273.15
    i_shot = math.sqrt(2.0 * E_CHARGE * float(I_A) * df)
    i_john = math.sqrt(4.0 * K_BOLTZMANN * T_K * df / float(R_ohm))
    return {
        "model_grade": "established",
        "I_shot_rms_A": i_shot,
        "I_johnson_rms_A": i_john,
        "shot_over_johnson": i_shot / i_john if i_john else float("inf"),
        "note": "Shot noise is negligible in this metal. The ratio is a comparison, not a measured current noise.",
    }


def hooge_current(I_A, n_carriers, alpha_H, f_lo, f_hi):
    """Optional Hooge-like 1/f current noise. Speculative. Default off."""
    if f_hi <= f_lo or n_carriers <= 0.0:
        raise ValueError("Hooge band and carrier count must be positive")
    var = (float(I_A) ** 2) * float(alpha_H) / float(n_carriers) * math.log(float(f_hi) / float(f_lo))
    return {
        "model_grade": "speculative",
        "alpha_H": float(alpha_H),
        "I_rms_A": math.sqrt(max(var, 0.0)),
        "note": "Speculative Hooge-like 1/f term. The coefficient is not a measurement on this coil.",
    }


def _pct(samples, nominal):
    a = np.asarray(samples, dtype=float)
    p05, p50, p95 = np.percentile(a, [5.0, 50.0, 95.0])
    return {
        "p05": float(p05),
        "p50": float(p50),
        "p95": float(p95),
        "mean": float(a.mean()),
        "bias_p50_minus_nominal": float(p50) - float(nominal),
    }


def _chip(name, value, limit, status, unit):
    return {"name": name, "value": float(value), "limit": float(limit), "status": status, "unit": unit}


def constraint_report(res, yield_pa, tsat_C, v_flow, aerated, insulator_max_C):
    """ok / warning / violation chips. Limit values are included."""
    chips = []

    def add_le(name, value, limit, warn_at, unit):
        if value > limit:
            status = "violation"
        elif value >= warn_at:
            status = "warning"
        else:
            status = "ok"
        chips.append(_chip(name, value, limit, status, unit))

    def add_ge(name, value, limit, warn_at, unit):
        if value < limit:
            status = "violation"
        elif value <= warn_at:
            status = "warning"
        else:
            status = "ok"
        chips.append(_chip(name, value, limit, status, unit))

    add_le("T_hot", res["T_hot_C"], 85.0, 70.0, "C")
    add_le("V_supply", res["V_total_V"], 8.0, 6.0, "V")
    if "homogeneity_ppm" in res:
        add_le("homogeneity", res["homogeneity_ppm"], 100.0, 90.0, "ppm")
    add_le("pressure", res["dp_Pa"], 5.0e5, 4.0e5, "Pa")
    add_ge("Reynolds", res["Re"], 1.0e4, 1.2e4, "1")
    add_ge("R1", res["R1"], 0.05, 0.052, "m")
    add_le("R2", res["R2"], 0.30, 0.29, "m")
    add_le("T_out_vs_Tsat_minus_10K", res["T_out_mixed_C"], tsat_C - 10.0, tsat_C - 15.0, "C")
    add_le("T_wall_vs_Tsat", res["T_wall_hot_C"], tsat_C, tsat_C - 20.0, "C")
    add_le("hoop_vs_0.6_yield", res["sigma_hoop_max_MPa"] * 1e6, 0.6 * yield_pa, 0.3 * yield_pa, "Pa")
    add_le("insulator_temperature", res["T_hot_C"], insulator_max_C, insulator_max_C - 30.0, "C")
    # velocity bands are warnings, encoded against the catalog edges
    if v_flow < 0.5:
        chips.append(_chip("velocity_fouling", v_flow, 0.5, "warning", "m/s"))
    elif (aerated and v_flow > 8.0) or v_flow > 15.0:
        lim = 8.0 if aerated and v_flow <= 15.0 else 15.0
        status = "violation" if v_flow > 15.0 else "warning"
        chips.append(_chip("velocity_erosion", v_flow, lim, status, "m/s"))
    else:
        chips.append(_chip("velocity", v_flow, 15.0, "ok", "m/s"))
    return chips


def velocity_warnings(v, aerated, Re):
    notes = []
    if v < 0.5:
        notes.append("Velocity below 0.5 m/s: fouling caution.")
    if Re < 1.0e4:
        notes.append("Re below 1e4: Dittus-Boelter is outside its validity range.")
    if aerated and v > 8.0:
        notes.append("Aerated water above about 8 m/s: copper erosion caution.")
    if v > 15.0:
        notes.append("Velocity above 15 m/s is outside this catalog.")
    return notes


def catalog_notes(conductor_id, coolant_id, insulator_id, housing, conductive_bore, radial_rho):
    notes = []
    cond = catalog.get_conductor(conductor_id)
    cool = catalog.get_coolant(coolant_id)
    ins = catalog.get_insulator(insulator_id)
    if conductor_id == "al_1350":
        notes.append("Warning: Al 1350 has lower strength, higher resistivity, and different water compatibility. Not recommended for this Bitter stack.")
    if conductor_id != "ofhc_cu":
        notes.append(cond["note"])
    notes.append(cool["note"])
    if cool.get("dielectric"):
        notes.append("Dielectric coolant, h will be lower.")
    if coolant_id == "water_glycol_30":
        notes.append("Water-glycol is electrically conductive. It is not a dielectric inside a Bitter hole.")
    notes.append(ins["note"])
    if insulator_id == "g10":
        notes.append("Warning: G-10 is a poor wet insulation choice.")
    notes.append(housing["note"])
    if housing["conductive"] and not conductive_bore:
        notes.append("Housing conductivity is ignored. Enable conductive bore to attach the shorted-turn warning.")
    if conductive_bore:
        if housing["conductive"]:
            notes.append("Conductive bore tube can act as a shorted turn and can support eddy currents. This is a warning, not a finite-element solution.")
        else:
            notes.append("Conductive bore was requested, but the selected bore wall is insulating. No shorted-turn path is added.")
    if radial_rho:
        notes.append("Radial rho(T) correction is reported as an approximate fractional field shift. It does not replace the locked continuum solution.")
    notes.append("Plate-to-plate contact resistance is sampled only in the emulation, not in the continuum resistance.")
    if coolant_id == "di_water":
        notes.append("ONB uses the Bergles-Rohsenow water correlation and is labeled approximate. There is no critical-heat-flux map.")
    return notes


def evaluate_with_catalog(design, conductor_id="ofhc_cu", coolant_id="di_water",
                          insulator_id="polyimide_kapton", housing_id="ss304",
                          temper="hard", T_amb=20.0, p_site_Pa=101325.0,
                          aerated=False, conductive_bore=False, radial_rho=False,
                          rho_scale=1.0, homogeneity=True):
    """Continuum evaluation plus catalog warnings. Default ids match evaluate_design."""
    cond = catalog.get_conductor(conductor_id)
    cool = catalog.get_coolant(coolant_id)
    ins = catalog.get_insulator(insulator_id)
    housing = catalog.get_housing(housing_id)
    yld = catalog.yield_pa(cond, temper)
    default_ids = (conductor_id == "ofhc_cu" and coolant_id == "di_water"
                   and rho_scale == 1.0 and not radial_rho)
    kwargs = {"homogeneity": homogeneity, "radial_rho": radial_rho}
    if not default_ids or conductor_id != "ofhc_cu":
        pass
    if conductor_id != "ofhc_cu":
        kwargs["rho20"] = cond["rho20"]
        kwargs["alpha"] = cond["alpha"]
        kwargs["k_solid"] = cond["k"]
        kwargs["dens"] = cond["density"]
        kwargs["cp_solid"] = cond["cp"]
    if coolant_id != "di_water":
        kwargs["fluid"] = (lambda T, cid=coolant_id: catalog.coolant_props(cid, T))
    if rho_scale != 1.0:
        kwargs["rho_scale"] = rho_scale
    # Default catalog ids call the untouched numeric path (rho_scale 1, no fluid).
    if conductor_id == "ofhc_cu" and coolant_id == "di_water" and rho_scale == 1.0:
        res = evaluate_design(design, homogeneity=homogeneity, radial_rho=radial_rho)
    else:
        res = evaluate_design(design, **kwargs)
    tsat = saturation_temperature_C(p_site_Pa)
    q_flux = res["q_flux_inner_W_m2"]
    dts = res["T_wall_hot_C"] - tsat
    q_onb = q_onb_bergles_rohsenow(p_site_Pa, dts) if coolant_id == "di_water" else None
    notes = list(res["warnings"])
    notes.extend(catalog_notes(conductor_id, coolant_id, insulator_id, housing, conductive_bore, radial_rho))
    notes.extend(velocity_warnings(design.v_flow, aerated, res["Re"]))
    if res["T_out_mixed_C"] > tsat - 10.0:
        notes.append("Mixed outlet is above Tsat - 10 K (Tsat %.2f C)." % tsat)
    if ins["t_min_m"] > design.d_ins:
        notes.append("Insulator thickness is below the catalog minimum of %.3g m." % ins["t_min_m"])
    if coolant_id != "di_water":
        notes.append("ONB correlation is the water Bergles-Rohsenow form and is not applied to this coolant.")
    # secondary housing term: reported, not subtracted from the water balance
    h_air = 8.0
    A_out = 2.0 * math.pi * design.R2 * design.L
    q_ext = h_air * A_out * max(res["T_cu_mean_C"] - float(T_amb), 0.0)
    t_wall = 5.0e-3
    mass_housing = housing["density"] * (2.0 * math.pi * design.R1 * t_wall * design.L)
    res = dict(res)
    res["warnings"] = notes
    res["catalog"] = {
        "conductor": conductor_id,
        "coolant": coolant_id,
        "insulator": insulator_id,
        "housing": housing_id,
        "temper": temper,
        "aerated": bool(aerated),
        "conductive_bore": bool(conductive_bore),
    }
    res["tsat_C"] = tsat
    res["q_onb_W_m2"] = q_onb
    res["onb_margin_K"] = -dts
    res["onb_model_grade"] = "approximate" if coolant_id == "di_water" else "not_applied"
    res["yield_Pa"] = yld
    res["housing_convection_W"] = q_ext
    res["housing_mass_kg"] = mass_housing
    res["housing_model_grade"] = "approximate"
    res["constraints"] = constraint_report(
        res, yld, tsat, design.v_flow, aerated, ins["t_max_C"])
    res["label"] = LABEL
    return res


def _fixed_current_B0(nominal, design):
    """B0 if the supply current stays at the nominal value while geometry changes."""
    n_turns = design.L / (design.d_plate + design.d_ins)
    NI = nominal["I_A"] * n_turns
    lnr = math.log(design.R2 / design.R1)
    C = NI / (design.L * lnr)
    return MU_0 * C * (math.asinh(design.L / (2.0 * design.R1)) - math.asinh(design.L / (2.0 * design.R2)))



def _eval_one(args):
    job, dr = args
    design = job["design"]
    R1, R2, L, d_plate, D_hole, T_in, v_flow, rho_lot, _rc = dr
    d_try = BitterDesign(
        R1=R1, R2=R2, L=L, d_plate=d_plate, d_ins=design.d_ins, D_hole=D_hole,
        v_flow=v_flow, pitch_factor=design.pitch_factor, T_in=T_in, B0=design.B0,
        dsv=design.dsv, K_minor=design.K_minor, eta_pump=design.eta_pump)
    B0_fixed_I = _fixed_current_B0({"I_A": job["nominal_I"]}, d_try)
    d_try = BitterDesign(
        R1=R1, R2=R2, L=L, d_plate=d_plate, d_ins=design.d_ins, D_hole=D_hole,
        v_flow=v_flow, pitch_factor=design.pitch_factor, T_in=T_in, B0=B0_fixed_I,
        dsv=design.dsv, K_minor=design.K_minor, eta_pump=design.eta_pump)
    ev = evaluate_with_catalog(
        d_try, job["conductor_id"], job["coolant_id"], job["insulator_id"], job["housing_id"], job["temper"],
        job["T_amb"], job["p_site_Pa"], job["aerated"], job["conductive_bore"], False,
        rho_scale=rho_lot, homogeneity=job["homogeneity"])
    keep = ("R_total_ohm", "V_total_V", "P_elec_W", "T_in", "T_hot_C", "B0_numeric_T", "homogeneity_ppm",
            "J_cu_inner_A_per_mm2", "rho_cu_mean", "T_cu_mean_C", "onb_margin_K", "sigma_hoop_max_MPa", "yield_Pa")
    return {k: ev[k] for k in keep if k in ev}, B0_fixed_I


def _map_realizations(job, draws, workers):
    """Serial map for workers <= 1, else a spawn-based process pool with 1 Numba thread per worker."""
    args = [(job, d) for d in draws]
    if workers is None or int(workers) <= 1 or len(draws) < 2:
        return [_eval_one(a) for a in args]
    import multiprocessing as mp
    ctx = mp.get_context("spawn")  # fork is unsafe once the OpenMP pool is live
    with ctx.Pool(int(workers), initializer=_single_thread) as pool:
        return pool.map(_eval_one, args, chunksize=max(1, len(args) // (4 * int(workers))))


def _single_thread():
    try:
        import numba
        numba.set_num_threads(1)
    except Exception:
        pass


def emulate(design=None, realizations=200, seed=None, conductor_id="ofhc_cu",
            coolant_id="di_water", insulator_id="polyimide_kapton", housing_id="ss304",
            temper="hard", T_amb=20.0, p_site_Pa=101325.0, aerated=False,
            conductive_bore=False, radial_rho=False, f_lo=1.0, f_hi=1.0e4,
            contact_median_ohm=1.0e-6, contact_sigma_ln=1.0, hooge=False,
            hooge_alpha=1.0e-3, kohler=True, kohler_a=1.0, homogeneity=True, workers=1):
    """Monte Carlo around a continuum design. Seed is required.

    workers > 1 evaluates the realizations in a process pool. All random draws happen first,
    in a fixed order, so the output is identical for any worker count.

    Scatter model (all approximate, stated in assumptions):
      geometry: normal, sigma 0.2 mm on R1 and R2, 0.5 mm on L
      hole diameter and plate thickness: uniform +/- 1 percent
      resistivity lot: uniform +/- 2 percent
      inlet sensor: normal, sigma 0.2 K
      inner-row flow: the single continuum velocity stands in for the inner
        row and is drawn uniformly within +/- 10 percent
      contact: log-normal per plate interface, median set by the caller
    """
    if seed is None:
        raise ValueError("seed is required")
    if realizations < 1:
        raise ValueError("realizations must be >= 1")
    if design is None:
        design = BitterDesign()
    cond = catalog.get_conductor(conductor_id)
    nominal = evaluate_with_catalog(
        design, conductor_id, coolant_id, insulator_id, housing_id, temper,
        T_amb, p_site_Pa, aerated, conductive_bore, radial_rho,
        homogeneity=homogeneity)
    J = nominal["J_cu_inner_A_per_mm2"] * 1.0e6
    drude = drude_rve(J, nominal["rho_cu_mean"], nominal["T_cu_mean_C"], nominal["B0"],
                      cond["n_m3"], design.d_plate, cond["v_fermi"],
                      kohler=False, kohler_a=kohler_a)
    drude_k = drude_rve(J, nominal["rho_cu_mean"], nominal["T_cu_mean_C"], nominal["B0"],
                        cond["n_m3"], design.d_plate, cond["v_fermi"],
                        kohler=kohler, kohler_a=kohler_a)
    langevin = langevin_ensemble(drude["v_d_m_s"], drude["tau_s"], nominal["T_cu_mean_C"], seed=seed)
    johnson = johnson_voltage(nominal["T_cu_mean_C"], nominal["R_total_ohm"], f_lo, f_hi)
    # Larmor note: density only, 1 Hz, not an MRI noise floor
    f_L = nominal["f_larmor_MHz"] * 1.0e6
    johnson_larmor = johnson_voltage(nominal["T_cu_mean_C"], nominal["R_total_ohm"], f_L, f_L + 1.0)
    johnson_larmor["note"] = ("Spectral density at the Larmor frequency, integrated over 1 Hz "
                               "only so a number exists. This is not the MRI noise floor.")
    shot = shot_vs_johnson(nominal["I_A"], nominal["T_cu_mean_C"], nominal["R_total_ohm"], f_lo, f_hi)
    volume = nominal["mass_cu_kg"] / cond["density"]
    n_carriers = cond["n_m3"] * volume
    hooge_rep = None
    if hooge:
        hooge_rep = hooge_current(nominal["I_A"], n_carriers, hooge_alpha, f_lo, f_hi)

    t, T = thermal.transient_lumped(
        nominal["P20_W"], ALPHA_CU if conductor_id == "ofhc_cu" else cond["alpha"],
        nominal["C_th_J_K"], nominal["R_th_K_W"], design.T_in,
        T0=nominal["T_cu_mean_C"], t_end=8000.0, n=8001, cooling=False)
    tsat = nominal["tsat_C"]
    hit = np.where(T >= tsat)[0]
    if hit.size:
        pump = {
            "reaches_tsat": True,
            "time_to_tsat_s": float(t[int(hit[0])]),
            "invalid_after": "Tsat",
            "model_grade": "approximate",
            "note": "Adiabatic lumped trajectory. Interpretation stops at Tsat. "
                    "The published time to 85 C is not a boiling model.",
        }
    else:
        pump = {
            "reaches_tsat": False,
            "time_to_tsat_s": None,
            "invalid_after": None,
            "model_grade": "approximate",
            "note": "Adiabatic lumped trajectory did not reach Tsat in 8000 s.",
        }

    rs = np.random.RandomState(int(seed))
    keys = {"B0": [], "V": [], "P": [], "T_hot": [], "ppm": [], "v_d": [], "mfp": [], "V_johnson": []}
    trip_hot = 0
    trip_onb = 0
    trip_stress = 0
    n_if_nom = max(1, int(round(nominal["n_turns"])) - 1)
    # phase 1: draw every realization serially (fixed RandomState order, so results do not depend
    # on the worker count); phase 2: evaluate in parallel; phase 3: aggregate in draw order.
    draws = []
    for _i in range(int(realizations)):
        R1 = design.R1 + rs.normal(0.0, 2.0e-4)
        R2 = design.R2 + rs.normal(0.0, 2.0e-4)
        L = design.L + rs.normal(0.0, 5.0e-4)
        R1 = max(0.03, R1)
        R2 = max(R1 + 0.02, R2)
        L = max(0.2, L)
        d_plate = design.d_plate * (1.0 + rs.uniform(-0.01, 0.01))
        D_hole = design.D_hole * (1.0 + rs.uniform(-0.01, 0.01))
        T_in = design.T_in + rs.normal(0.0, 0.2)
        v_flow = design.v_flow * (1.0 + rs.uniform(-0.10, 0.10))
        v_flow = max(0.05, v_flow)
        rho_lot = 1.0 + rs.uniform(-0.02, 0.02)
        n_if = max(1, int(round(L / (d_plate + design.d_ins))) - 1)
        r_each = rs.lognormal(math.log(contact_median_ohm), contact_sigma_ln, size=n_if)
        draws.append((R1, R2, L, d_plate, D_hole, T_in, v_flow, rho_lot, float(np.sum(r_each))))
    job = dict(design=design, nominal_I=nominal["I_A"], conductor_id=conductor_id, coolant_id=coolant_id,
               insulator_id=insulator_id, housing_id=housing_id, temper=temper, T_amb=T_amb,
               p_site_Pa=p_site_Pa, aerated=aerated, conductive_bore=conductive_bore, homogeneity=homogeneity)
    evals = _map_realizations(job, draws, workers)
    for (R1, R2, L, d_plate, D_hole, T_in, v_flow, rho_lot, r_contact), (ev, B0_fixed_I) in zip(draws, evals):
        scale = (ev["R_total_ohm"] + r_contact) / ev["R_total_ohm"]
        V = ev["V_total_V"] * scale
        P = ev["P_elec_W"] * scale
        T_hot = ev["T_in"] + (ev["T_hot_C"] - ev["T_in"]) * scale
        B0 = ev["B0_numeric_T"] if homogeneity else B0_fixed_I
        ppm = ev["homogeneity_ppm"] if homogeneity else float("nan")
        J_i = ev["J_cu_inner_A_per_mm2"] * 1.0e6
        dru = drude_rve(J_i, ev["rho_cu_mean"], ev["T_cu_mean_C"], B0, cond["n_m3"],
                        d_plate, cond["v_fermi"], kohler=kohler, kohler_a=kohler_a)
        john = johnson_voltage(ev["T_cu_mean_C"], ev["R_total_ohm"] + r_contact, f_lo, f_hi)
        keys["B0"].append(B0)
        keys["V"].append(V)
        keys["P"].append(P)
        keys["T_hot"].append(T_hot)
        keys["ppm"].append(ppm)
        keys["v_d"].append(dru["v_d_m_s"])
        keys["mfp"].append(dru["mean_free_path_m"])
        keys["V_johnson"].append(john["V_rms_V"])
        if T_hot > 85.0:
            trip_hot += 1
        if ev["onb_margin_K"] < 0.0:
            trip_onb += 1
        if ev["sigma_hoop_max_MPa"] * 1e6 > 0.6 * ev["yield_Pa"]:
            trip_stress += 1

    nreal = float(realizations)
    percentiles = {
        "B0_T": _pct(keys["B0"], nominal.get("B0_numeric_T", nominal["B0"])),
        "V_V": _pct(keys["V"], nominal["V_total_V"]),
        "P_W": _pct(keys["P"], nominal["P_elec_W"]),
        "T_hot_C": _pct(keys["T_hot"], nominal["T_hot_C"]),
        "ppm": _pct(keys["ppm"], nominal.get("homogeneity_ppm", float("nan"))),
        "v_d_m_s": _pct(keys["v_d"], drude["v_d_m_s"]),
        "mean_free_path_m": _pct(keys["mfp"], drude["mean_free_path_m"]),
        "V_johnson_V": _pct(keys["V_johnson"], johnson["V_rms_V"]),
    }
    return {
        "label": LABEL,
        "model_grade": dict(MODEL_GRADE),
        "seed": int(seed),
        "n_realizations": int(realizations),
        "nominal_summary": {
            "B0_T": nominal.get("B0_numeric_T", nominal["B0"]),
            "V_V": nominal["V_total_V"],
            "P_W": nominal["P_elec_W"],
            "T_hot_C": nominal["T_hot_C"],
            "ppm": nominal.get("homogeneity_ppm"),
            "R_ohm": nominal["R_total_ohm"],
            "I_A": nominal["I_A"],
            "v_d_m_s": drude["v_d_m_s"],
            "mean_free_path_m": drude["mean_free_path_m"],
            "tau_s": drude["tau_s"],
            "omega_c_tau": drude_k["omega_c_tau"],
            "V_johnson_V": johnson["V_rms_V"],
        },
        "drude_kohler_off": drude,
        "drude": drude_k,
        "langevin": langevin,
        "johnson": johnson,
        "johnson_larmor_1Hz": johnson_larmor,
        "shot_vs_johnson": shot,
        "hooge": hooge_rep,
        "pump_failure": pump,
        "thermal": {
            "tsat_C": tsat,
            "q_onb_W_m2": nominal["q_onb_W_m2"],
            "onb_margin_K": nominal["onb_margin_K"],
            "onb_correlation": "Bergles-Rohsenow 1964" if coolant_id == "di_water" else "not applied",
            "model_grade": nominal["onb_model_grade"],
        },
        "warnings": nominal["warnings"],
        "constraints": nominal["constraints"],
        "assumptions": {
            "hole_diameter_tolerance_rel": 0.01,
            "plate_thickness_tolerance_rel": 0.01,
            "radius_sigma_m": 2.0e-4,
            "length_sigma_m": 5.0e-4,
            "tin_sigma_K": 0.2,
            "rho_lot_half_width": 0.02,
            "flow_maldistribution_half_width": 0.10,
            "contact_median_ohm_per_interface": float(contact_median_ohm),
            "contact_log_sigma": float(contact_sigma_ln),
            "n_interfaces_nominal": n_if_nom,
            "contact_note": "Log-normal plate-to-plate resistance. Median is an assumption, not a measurement.",
            "flow_note": "The continuum model has one velocity. Perturbing it stands in for inner-row maldistribution.",
        },
        "percentiles": percentiles,
        "trip_probability": {
            "T_hot_above_85C": trip_hot / nreal,
            "onb_margin_below_0": trip_onb / nreal,
            "stress_above_0.6_yield": trip_stress / nreal,
        },
        "active_constraints": [c["name"] for c in nominal["constraints"] if c["status"] != "ok"],
    }
