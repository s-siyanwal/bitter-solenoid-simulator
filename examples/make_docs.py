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
 ("R27", "Lumped thermal time constant (63 %) / pump-failure time to 85 C", "%.1f s / %.0f s" % (x["transient_tau63_s"], x["pump_failure_time_to_85C_s"])),
 ("R28", "Copper mass", "%.0f kg" % o["mass_cu_kg"]),
 ("R29", "Larmor frequency at 0.5 T", "%.6f MHz" % x["larmor_MHz"]),
 ("R30", "Swiss roll tuned exactly to f_L: mu_eff(f_L), Q", "%.3f + %.2fj, Q = %.1f" % (sr["mu_at_fL_when_tuned_to_fL"][0], sr["mu_at_fL_when_tuned_to_fL"][1], sr["Q"])),
 ("R31", "Best roll tuning f0/f_L, mu_eff(f_L), heuristic SNR gain", "%.5f, %.3f + %.3fj, %.3f" % (sr["best_f0_over_fL"], sr["mu_at_fL_best"][0], sr["mu_at_fL_best"][1], sr["snr_gain_best"])),
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
rmap = dict((r[0], r[2]) for r in R + VR)

def rtable(rows):
    return "| ID | quantity | value |\n|---|---|---|\n" + "\n".join("| %s | %s | %s |" % r for r in rows)

summary = open(os.path.join(ROOT, "examples", "summary_template.md")).read()
summary = summary.replace("{{RESULTS_TABLE}}", rtable(R + VR)).replace("{{ENV}}", envs)
open(os.path.join(ROOT, "DESIGN_SUMMARY.md"), "w").write(summary)

readme = open(os.path.join(ROOT, "examples", "readme_template.md")).read()
readme = readme.replace("{{RESULTS_TABLE}}", rtable(R)).replace("{{VALIDATION_TABLE}}", rtable(VR)).replace("{{ENV}}", envs)
readme = readme.replace("{{OPT_SECONDS}}", "%.0f" % D["optimiser"]["seconds"]).replace("{{NFEV}}", str(D["optimiser"]["de_nfev"]))
open(os.path.join(ROOT, "README.md"), "w").write(readme)
print("ok")
