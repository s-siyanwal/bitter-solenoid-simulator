"""Generate README.md and DESIGN_SUMMARY.md from results/*.json (no hand-typed numbers)."""
import json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = json.load(open(os.path.join(ROOT, "results", "results.json")))
V = json.load(open(os.path.join(ROOT, "results", "validation.json")))
o, b, x, env = D["optimal"], D["baseline_pdf_initial_guess"], D["extra"], D["env"]
sr = x["swissroll"]
envs = ", ".join("%s %s" % (k, env[k]) for k in ("python", "numpy", "scipy", "numba", "matplotlib"))

R = [  # (id, label, value string)
 ("R1", "Optimal inner radius R1", "%.2f mm" % (o["R1"] * 1e3)),
 ("R2", "Optimal outer radius R2", "%.2f mm" % (o["R2"] * 1e3)),
 ("R3", "Optimal stack length L", "%.1f mm" % (o["L"] * 1e3)),
 ("R4", "Plate thickness / insulator", "%.3f mm / %.2f mm (assumed)" % (o["d_plate"] * 1e3, o["d_ins"] * 1e3)),
 ("R5", "Number of plates (turns)", "%.1f" % o["n_turns"]),
 ("R6", "Cooling holes: diameter, pitch/D, count", "%.3f mm, %.3f, %d holes in %d rows" % (o["D_hole"] * 1e3, o["pitch_factor"], o["n_holes"], o["n_rows"])),
 ("R7", "Ampere-turns NI for Bz(0) = 0.5 T", "%.0f A" % o["NI_At"]),
 ("R8", "Current I", "%.1f A" % o["I_A"]),
 ("R9", "Supply voltage / voltage per plate V0", "%.3f V / %.3f mV" % (o["V_total_V"], o["V_per_turn_V"] * 1e3)),
 ("R10", "Resistance", "%.4f mOhm" % (o["R_total_ohm"] * 1e3)),
 ("R11", "Copper current density at R1 / R2", "%.3f / %.3f A/mm^2" % (o["J_cu_inner_A_per_mm2"], o["J_cu_outer_A_per_mm2"])),
 ("R12", "Electrical power (rho at mean Cu temperature)", "%.3f kW" % (o["P_elec_W"] / 1e3)),
 ("R13", "Pump shaft power (eta = 0.7)", "%.1f W" % o["P_pump_W"]),
 ("R14", "Isocentre field from exact loop model", "%.9f T" % o["B0_numeric_T"]),
 ("R15", "Homogeneity, peak-to-peak over 30 mm DSV", "%.2f ppm" % o["homogeneity_ppm"]),
 ("R16", "Water velocity / total flow", "%.3f m/s / %.1f L/min (%.2f kg/s)" % (o["v_flow"], o["flow_L_min"], o["m_dot_kg_s"])),
 ("R17", "Reynolds / Prandtl number", "%.0f / %.3f" % (o["Re"], o["Pr"])),
 ("R18", "h Dittus-Boelter / Gnielinski", "%.0f / %.0f W/m^2K" % (o["h_W_m2K"], o["h_gnielinski_W_m2K"])),
 ("R19", "Pressure drop (Darcy f = %.4f)" % o["f_darcy"], "%.2f kPa" % (o["dp_Pa"] / 1e3)),
 ("R20", "Mixed water outlet temperature (inlet 20 C)", "%.3f C" % o["T_out_mixed_C"]),
 ("R21", "Hottest channel: water rise / film dT / conduction dT", "%.3f / %.3f / %.3f K" % (o["dT_water_inner_K"], o["dT_film_inner_K"], o["dT_cond_inner_K"])),
 ("R22", "Maximum (hot-spot) copper temperature, at R1, outlet", "%.3f C" % o["T_hot_C"]),
 ("R23", "Mean copper temperature", "%.3f C" % o["T_cu_mean_C"]),
 ("R24", "Max hoop stress (E21, field profile) / simple bound C B0/lambda", "%.4f / %.4f MPa" % (x["hoop_stress_max_MPa_profile"], o["sigma_hoop_max_MPa"])),
 ("R25", "Axial compressive force on each half", "%.0f N" % abs(x["axial_compressive_force_N"])),
 ("R26", "Inductance / stored energy / L/R time constant", "%.4f mH / %.0f J / %.3f s" % (x["inductance_H"] * 1e3, x["stored_energy_J"], x["L_over_R_s"])),
 ("R27", "Lumped thermal time constant (63 %) / pump-failure time to 85 C: lumped mean / hot-spot adiabatic bound", "%.1f s / %.0f s / %.0f s" % (x["transient_tau63_s"], x["pump_failure_time_to_85C_s"], x["pump_failure_hotspot_adiabatic_s"])),
 ("R28", "Copper mass", "%.0f kg" % o["mass_cu_kg"]),
 ("R29", "Larmor frequency at 0.5 T", "%.6f MHz" % x["larmor_MHz"]),
 ("R30", "Swiss roll tuned exactly to f_L: mu_eff(f_L), Q", "%.3f + %.2fj, Q = %.1f" % (sr["mu_at_fL_when_tuned_to_fL"][0], sr["mu_at_fL_when_tuned_to_fL"][1], sr["Q"])),
 ("R31", "Best roll tuning f0/f_L, mu_eff(f_L), legacy heuristic SNR gain (E28, superseded by R49)", "%.5f, %.3f + %.3fj, %.3f" % (sr["best_f0_over_fL"], sr["mu_at_fL_best"][0], sr["mu_at_fL_best"][1], sr["snr_gain_best"])),
 ("R32", "Swiss roll mu_eff at DC", "%.1f + %.1fj (no static-field effect)" % tuple(sr["mu_at_DC"])),
 ("R33", "PDF initial guess (R2 = 0.15 m, L = 0.8 m, v = 2.5 m/s, 2 mm plates)", "P = %.3f kW, V = %.3f V (violates 8 V), %.2f ppm (violates 100 ppm), T_hot = %.2f C" % (b["P_elec_W"] / 1e3, b["V_total_V"], b["homogeneity_ppm"], b["T_hot_C"])),
 ("R34", "Energy-balance residual (sum of cell heats vs E31)", "%.1e" % o["energy_residual"]),
]
S = V["segments"]
VR = [
 ("R35", "Richardson order of convergence, segment engine (theory 2)", "%.6f" % S["richardson_p"]),
 ("R36", "Segment engine vs E3 (thick uniform), 360 angular steps", "%.3e relative" % S["rows_vs_E3"][3][1]),
 ("R37", "Bitter 1/r helical filaments vs E4 (691200 segments)", "%.3e relative" % V["bitter_segments"]["helical"]["rel_err"]),
 ("R38", "Loop model vs E3 / E4 closed forms (61 axial points)", "%.1e / %.1e" % (V["thick"]["uniform_loops_vs_E3_max_rel"], V["thick"]["bitter_loops_vs_E4_max_rel"])),
 ("R39", "Long-solenoid limit L/R = 1000 vs mu0 n I", "%.2e relative" % V["long_solenoid"][2]["rel_vs_ideal"]),
 ("R40", "Numerical sheet inductance vs Nagaoka (worst of 3)", "%.1e relative" % max(r["rel"] for r in V["inductance"]["sheet"])),
]
F = json.load(open(os.path.join(ROOT, "results", "fmri.json")))
s3, s4, hp, st, rq, rf = F["shim_30mm"], F["shim_40mm"], F["head_preset"], F["stability"], F["stability_requirements"], F["rf"]
cb, vb = st["current"]["budget_ppm"], st["voltage"]["budget_ppm"]
tol = F["tolerance_30mm"]["tesseral_ppm"]
FR = [
 ("R41", "Zonal Z2 / Z4 over 30 mm DSV (ppm at DSV radius, unshimmed)", "%.3f / %.4f" % (s3["zonal_ppm_unshimmed"][2], s3["zonal_ppm_unshimmed"][4])),
 ("R42", "Peak-to-peak over 30 mm / 40 mm DSV: unshimmed -> Z2+Z4 shim pairs", "%.2f -> %.3f ppm / %.2f -> %.3f ppm" % (s3["ppm_unshimmed"], s3["ppm_shimmed"], s4["ppm_unshimmed"], s4["ppm_shimmed"])),
 ("R43", "Shim pairs: radius, z positions, NI, power (J = 2 A/mm^2 assumed)", "%.1f mm, +-%.1f / +-%.1f mm, %.2f / %.2f A, %.3f W" % (s3["shim_radius_m"] * 1e3, s3["shim_z_m"][0] * 1e3, s3["shim_z_m"][1] * 1e3, s3["shim_NI_A"][0], s3["shim_NI_A"][1], s3["shim_power_W"])),
 ("R44", "Tesseral terms from 0.5 mm offset + 1 mrad tilt (assumed tolerance): A11 / B21", "%.3f / %.4f ppm" % (tol["A11"], tol["B21"])),
 ("R45", "Head preset (R1 = 190 mm, 200 mm DSV, <= 10 ppm after Z2/Z4; best on the old R2 <= 0.7 m grid, see R60 for the optimiser)", "R2 %.0f mm, L %.0f mm, P %.1f kW + pump %.2f kW, %.2f V, %.0f A, %.1f t Cu, %.0f L/min, %.0f -> %.2f ppm, T_hot %.2f C%s" % (hp["R2"] * 1e3, hp["L"] * 1e3, hp["P_elec_W"] / 1e3, hp["P_pump_W"] / 1e3, hp["V_total_V"], hp["I_A"], hp["mass_cu_kg"] / 1e3, hp["flow_L_min"], hp["ppm_unshimmed"], hp["ppm_shimmed"], hp["T_hot_C"], " (violates 8 V)" if hp["violates_8V_supply"] else "")),
 ("R46", "Field drift vs copper temperature at fixed current (numeric, isotropic expansion)", "%.3f ppm/K" % (F["alpha_B_per_K_numeric"] * 1e6)),
 ("R47", "B0 stability over 10 min, current-regulated PSU: ripple / drift / expansion / total", "%.3f / %.3f / %.3f / %.3f ppm (%.1f Hz at f_L); voltage-regulated total %.1f ppm" % (cb["psu_ripple"], cb["psu_drift"], cb["thermal_expansion"], cb["total"], st["current"]["budget_Hz"]["total"], vb["total"])),
 ("R48", "Max copper dT / water oscillation amplitude for 1 ppm (current / voltage mode)", "%.4f / %.2e K ; %.4f / %.2e K" % (rq["max_copper_dT_K_current_mode"], rq["max_copper_dT_K_voltage_mode"], rq["max_water_amp_K_current_mode"], rq["max_water_amp_K_voltage_mode"])),
 ("R49", "RF SNR, Swiss-roll slab vs same coil over air / vs coil on tissue (loss x50, best tuning f0/f_L)", "%.3f / %.3f (f0/f_L = %.4f, mu = %.2f + %.2fj)" % (rf["gain_vs_air"], rf["gain_vs_contact"], rf["detune"], rf["mu"][0], rf["mu"][1])),
 ("R50", "RF SNR vs air for loss multiplier 1 / 10 / 50", "%.3f / %.3f / %.3f" % (rf["gain_vs_air_by_loss_multiplier"]["1.0"], rf["gain_vs_air_by_loss_multiplier"]["10.0"], rf["gain_vs_air"])),
 ("R51", "RF resistances with slab: coil / tissue / slab", "%.4f / %.4f / %.4f ohm" % (rf["slab"]["R_coil"], rf["slab"]["R_tissue"], rf["slab"]["R_slab"])),
]

PJ = json.load(open(os.path.join(ROOT, "results", "particles.json")))
_pc = dict((r["key"], r) for r in PJ["comparison"])
_pb = PJ["benchmark"]
PR = [
 ("R52", "Particle emulation: speedup on %d threads (field / carriers / walkers), serial == parallel bit-identical" % PJ["threads"],
  "%.2fx / %.2fx / %.2fx, %s" % (_pb["field"]["speedup"], _pb["carriers"]["speedup"], _pb["walkers"]["speedup"], "yes" if all(v["bit_identical"] for v in _pb.values()) else "NO")),
 ("R53", "Particle vs continuum: B0, P (relative difference)", "%.2e / %.2e" % (_pc["B0"]["rel_err"], _pc["P"]["rel_err"])),
 ("R54", "Particle vs continuum: DSV homogeneity", "%.2f ppm vs %.2f ppm" % (_pc["ppm"]["particle"], _pc["ppm"]["continuum"])),
 ("R55", "Carrier (Bloch-Gruneisen) vs linear resistivity at mean Cu T", "%+.3f %% (+/- %.3f %% stat.)" % (100 * _pc["rho_mean"]["rel_err"], 100 * _pc["rho_mean"]["stderr_rel"])),
 ("R56", "Hot-spot copper T: particle / continuum", "%.3f / %.3f C" % (_pc["T_hot"]["particle"], _pc["T_hot"]["continuum"])),
]


def _g(v):
    return "%.6g" % v


def ptable():
    h = "| quantity | continuum | particle | rel. difference | stat. error (1 sigma) | expected |\n|---|---|---|---|---|---|\n"
    rows = []
    for r in PJ["comparison"]:
        se = "-" if r["stderr_rel"] is None else ("0 (deterministic)" if r["stderr_rel"] == 0 else "%.2e" % r["stderr_rel"])
        rows.append("| %s [%s] | %s | %s | %+.2e | %s | %s |" % (r["label"], r["unit"], _g(r["continuum"]), _g(r["particle"]), r["rel_err"], se, r["expect"]))
    return h + "\n".join(rows)


def btable():
    h = "| kernel | N | serial [s] | parallel [s] | speedup | identical |\n|---|---|---|---|---|---|\n"
    names = [("field", "P1 current elements (Biot-Savart, 404 DSV points)"), ("moments", "P1 moments (NI, P)"),
             ("carriers", "P2 carriers (Green-Kubo)"), ("walkers", "P3 heat walkers (Feynman-Kac)"),
             ("tolerance_mc_processes", "emulate() tolerance MC, 8 processes")]
    return h + "\n".join("| %s | %d | %.3f | %.3f | %.2fx | %s |" % (lab, _pb[k]["n"], _pb[k]["serial_s"], _pb[k]["parallel_s"], _pb[k]["speedup"], "yes" if _pb[k]["bit_identical"] else "NO") for k, lab in names)

HX = json.load(open(os.path.join(ROOT, "results", "helical.json")))
HD = json.load(open(os.path.join(ROOT, "results", "head.json")))
_h30, _h40 = HX["dsv"]["30mm"], HX["dsv"]["40mm"]


def _lad(a, i):
    return a["ladder_ppm"][i][1]


def _top(a, n=4):
    t = sorted(a["tesseral_ppm"].items(), key=lambda kv: -abs(kv[1]))[:n]
    return ", ".join("%s %+.2f" % (k, v) for k, v in t) or "none >= 0.01"


_LAD = [("+ retune Z1, Z2, Z4", ["Z"]), ("+ X, Y", ["X", "Y"]), ("+ ZX, ZY", ["ZX", "ZY"]), ("+ X2-Y2, XY", ["X2-Y2", "XY"]),
        ("+ all n = 3 tesseral", ["n=3"]), ("+ all n = 4 tesseral", ["n=4"])]


def _contrib(a):
    """Orders in the ladder prefix that reach the target whose step lowers p-p by >= 0.1 ppm."""
    out, L_ = [], a["ladder_ppm"]
    for i in range(2, len(L_)):
        if L_[i - 1][1] - L_[i][1] >= 0.1:
            out += dict(_LAD)[L_[i][0]]
        if L_[i][0] == a["ladder_reached"]:
            break
    return out


def htable():
    h = ("| DSV | path variant | p-p unshimmed | p-p after Z2/Z4 pairs | largest tesseral terms of abs(B) [ppm at r0] | "
         "max B_perp total / helical part [uT] | ideal shims needed (ladder) | terms >= 1 ppm |\n|---|---|---|---|---|---|---|---|\n")
    rows = []
    for key, d in (("30 mm", _h30), ("40 mm", _h40)):
        for kind in ("ideal", "uniform", "aligned", "rotating"):
            a = d[kind]
            need = "-" if kind == "ideal" else ("not reached with n <= 4" if a.get("ladder_orders_needed") is None else (", ".join(_contrib(a)) or "none"))
            rows.append("| %s | %s | %.2f | %.2f | %s | %.0f / %.0f | %s | %s |" % (
                key, kind, _lad(a, 0), _lad(a, 1), "-" if kind == "ideal" else _top(a), a["B_perp_max_uT"], a["B_perp_helical_max_uT"],
                need, "-" if kind == "ideal" else (", ".join(a["needed_shims"]) or "none")))
    return h + "\n".join(rows)


def ladtable():
    names = [n for n, _ in _h30["ideal"]["ladder_ppm"]]
    h = "| step | " + " | ".join("%s %s" % (k, dsv) for dsv in ("30 mm", "40 mm") for k in ("ideal", "uniform", "aligned", "rotating")) + " |\n|---|" + "---|" * 8 + "\n"
    rows = []
    for i, nm in enumerate(names):
        vals = [d[k]["ladder_ppm"][i][1] for d in (_h30, _h40) for k in ("ideal", "uniform", "aligned", "rotating")]
        rows.append("| %s | " % nm + " | ".join("%.2f" % v for v in vals) + " |")
    return h + "\n".join(rows)


def headtable():
    h = ("| run | seed | R2 [mm] | L [mm] | plate [mm] | hole [mm] | v [m/s] | pitch/D | P_elec + pump + shim [kW] | V [V] | I [A] | "
         "Cu [t] | p-p shimmed [ppm] | T_hot [C] | bounds hit | active constraints |\n|---|" + "---|" * 15 + "\n")
    lab = {"power_8V": "min power, V <= 8 V", "power_noV": "min power, no V limit", "min_voltage": "min supply voltage"}
    rows = []
    for k in ("power_8V", "power_noV", "min_voltage"):
        for r in HD["runs"][k]:
            x = r["x"]
            rows.append("| %s | %d | %.0f | %.0f | %.2f | %.2f | %.3f | %.2f | %.2f | %.3f | %.0f | %.1f | %.2f | %.2f | %s | %s |" % (
                lab[k], r["seed"], x["R2"] * 1e3, x["L"] * 1e3, x["d_plate"] * 1e3, x["D_hole"] * 1e3, x["v_flow"], x["pitch_factor"],
                r["P_total_W"] / 1e3, r["V_total_V"], r["I_A"], r["mass_cu_kg"] / 1e3, r["ppm_shimmed"], r["T_hot_C"],
                ", ".join(r["at_bound"]) or "none (interior)", ", ".join(r["active_constraints"]) or "none"))
    return h + "\n".join(rows)


_hb = min(HD["runs"]["power_8V"], key=lambda r: r["P_total_W"])
_hv = min(HD["runs"]["min_voltage"], key=lambda r: r["V_total_V"])
_spread = (max(r["P_total_W"] for r in HD["runs"]["power_8V"]) - _hb["P_total_W"]) / _hb["P_total_W"]
HR = [
 ("R57", "Helical path, 30 mm DSV, p-p after Z2/Z4 pairs: ideal / uniform helix / aligned slits / rotating slits",
  "%.2f / %.2f / %.2f / %.2f ppm" % (_lad(_h30["ideal"], 1), _lad(_h30["uniform"], 1), _lad(_h30["aligned"], 1), _lad(_h30["rotating"], 1))),
 ("R58", "Helical path, 40 mm DSV, p-p after Z2/Z4 pairs: ideal / uniform / aligned / rotating",
  "%.2f / %.2f / %.2f / %.2f ppm" % (_lad(_h40["ideal"], 1), _lad(_h40["uniform"], 1), _lad(_h40["aligned"], 1), _lad(_h40["rotating"], 1))),
 ("R59", "Aligned slits, 30 mm DSV: X / Y terms of abs(B); max transverse field (total, with return bus)",
  "%.2f / %.2f ppm; %.0f uT" % (_h30["aligned"]["tesseral_ppm"].get("A11", 0.0), _h30["aligned"]["tesseral_ppm"].get("B11", 0.0), _h30["aligned"]["B_perp_max_uT"])),
 ("R60", "Head preset optimum (DE, best of seeds 1-3, V <= 8 V): R2, L, plate, total power, V, Cu mass",
  "%.0f mm, %.0f mm, %.2f mm, %.2f kW, %.3f V, %.1f t (seed spread in power %.2f %%)" % (_hb["x"]["R2"] * 1e3, _hb["x"]["L"] * 1e3, _hb["x"]["d_plate"] * 1e3, _hb["P_total_W"] / 1e3, _hb["V_total_V"], _hb["mass_cu_kg"] / 1e3, 100 * _spread)),
 ("R61", "Head preset: minimum achievable supply voltage (other constraints kept)", "%.3f V (%s)" % (_hv["V_total_V"], ", ".join(_hv["at_bound"]) or "interior")),
]

rmap = dict((r[0], r[2]) for r in R + VR + FR + PR + HR)

def rtable(rows):
    return "| ID | quantity | value |\n|---|---|---|\n" + "\n".join("| %s | %s | %s |" % r for r in rows)

summary = open(os.path.join(ROOT, "examples", "summary_template.md")).read()
summary = summary.replace("{{RESULTS_TABLE}}", rtable(R + VR + FR + PR + HR)).replace("{{ENV}}", envs)
import sys
CHECK = "--check" in sys.argv
def _emit(name, text):
    path = os.path.join(ROOT, name)
    if CHECK:
        if open(path).read() != text:
            raise SystemExit("%s is stale: run examples/make_docs.py" % name)
    else:
        open(path, "w").write(text)
_emit("DESIGN_SUMMARY.md", summary)

readme = open(os.path.join(ROOT, "examples", "readme_template.md")).read()
readme = readme.replace("{{RESULTS_TABLE}}", rtable(R)).replace("{{FMRI_TABLE}}", rtable(FR)).replace("{{VALIDATION_TABLE}}", rtable(VR)).replace("{{ENV}}", envs)
_pp = PJ["production"]
readme = readme.replace("{{PARTICLES_ROWS}}", rtable(PR)).replace("{{PARTICLES_TABLE}}", ptable()).replace("{{PARTICLES_BENCH}}", btable())
readme = readme.replace("{{P_FIELD_N}}", str(_pp["field_N"])).replace("{{P_CARRIERS_N}}", str(_pp["carriers_N"])).replace("{{P_WALKERS_N}}", str(_pp["walkers_N"]))
readme = readme.replace("{{P_AXIS_INNER}}", "%.1e" % _pp["axis_rel_max_inner"]).replace("{{P_AXIS_ALL}}", "%.1e" % _pp["axis_rel_max_all"])
readme = readme.replace("{{P_MC_MS}}", "%.1f" % (1e3 * _pb["tolerance_mc_processes"]["serial_s"] / _pb["tolerance_mc_processes"]["n"]))
readme = readme.replace("{{P_B0_ERR}}", "%.1e" % _pc["B0"]["rel_err"]).replace("{{P_PPM_SAME}}", "%.2f" % _pp["ppm_loop_same_points"])
readme = readme.replace("{{P_RHO_MODEL}}", "%+.1e" % _pc["rho_bg_model"]["rel_err"]).replace("{{P_RHO_SE}}", "%.1e" % _pc["rho_mean"]["stderr_rel"]).replace("{{P_RHO_85}}", "%+.2f %%" % (100 * _pc["rho_85"]["rel_err"]))
readme = readme.replace("{{P_PPM_R15}}", "%.2f" % _pp["ppm_R15"]).replace("{{P_THETA}}", "%.0f" % _pp["theta_R"]).replace("{{P_RRR}}", "%.0f" % _pp["RRR"])
readme = readme.replace("{{HELICAL_TABLE}}", htable()).replace("{{LADDER_TABLE}}", ladtable()).replace("{{HEAD_TABLE}}", headtable()).replace("{{NEXT_ROWS}}", rtable(HR))
_ha = HX["assumptions"]
readme = readme.replace("{{HX_NT}}", str(_ha["n_turns"])).replace("{{HX_OV}}", "%.0f" % _ha["overlap_deg"]).replace("{{HX_BUS}}", "%.0f" % (_ha["bus_radius_m"] * 1e3)).replace("{{HX_NR}}", str(_ha["radial_filaments"]))
readme = readme.replace("{{HX_SCALE}}", "%.4f" % _h30["rotating"]["current_scale"]).replace("{{HX_RICH}}", "%.2f" % max(max(d[k]["tesseral_richardson_err_ppm"].values() or [0]) for d in (_h30, _h40) for k in ("uniform", "aligned", "rotating")))
readme = readme.replace("{{HX_RICH_AL}}", "%.3f" % max(max(d["aligned"]["tesseral_richardson_err_ppm"].values() or [0]) for d in (_h30, _h40)))
readme = readme.replace("{{HEAD_MIN_SECONDS}}", "%.0f" % (HD["seconds"] / 60.0))
readme = readme.replace("{{OPT_SECONDS}}", "%.0f" % D["optimiser"]["seconds"]).replace("{{NFEV}}", str(D["optimiser"]["de_nfev"]))
_emit("README.md", readme)
print("ok")
