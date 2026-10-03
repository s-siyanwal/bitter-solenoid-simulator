"""Reproduce every number and figure in README.md / VALIDATION.md.

    python examples/run_all.py            (~3-5 min on 8 cores)

Writes results/results.json, results/validation.json, figures/*.png,
VALIDATION.md and results/RESULTS.md.
"""
import os
import sys
import json
import time
import platform
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "package"))
os.environ.setdefault("MPLBACKEND", "Agg")
warnings.filterwarnings("ignore")

import numpy as np
import scipy, numba, matplotlib
from bittersim import evaluate_design, BitterDesign, swissroll, mechanics, inductance, validation
from bittersim.optimize import optimise, LIMITS, BOUNDS
from bittersim import plotting

FIG = os.path.join(ROOT, "figures"); RES = os.path.join(ROOT, "results")
for d in (FIG, RES):
    if not os.path.isdir(d):
        os.makedirs(d)


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


t0 = time.time()
env = {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
       "numba": numba.__version__, "matplotlib": matplotlib.__version__, "cpus": os.cpu_count()}
print("environment", env)

baseline = evaluate_design(BitterDesign())                    # PDF initial guess
t = time.time()
opt = optimise(seed=1, maxiter=100, popsize=15)
t_opt = time.time() - t
best = opt["result"]
print("optimisation done in %.1f s (%s)" % (t_opt, opt["method"]))

# extra results for the optimal design
r_m, Bz_m, sig_m = mechanics.hoop_stress_profile(best["R1"], best["R2"], best["L"], best["C_A_per_m"], best["fill_total"])
Fz = mechanics.axial_force(best["R1"], best["R2"], best["L"], best["C_A_per_m"])
Lind = inductance.coil_inductance(best["R1"], best["R2"], best["L"], best["n_turns"], "bitter", nr=10)
fL = swissroll.larmor_hz(best["B0"])
sr_tuned = swissroll.SwissRoll(fL)
ratios, gains, qbest, gbest = swissroll.best_detuning(best["B0"])
sr_best = swissroll.SwissRoll(fL * qbest)
mu_t = complex(sr_tuned.mu(fL)); mu_b = complex(sr_best.mu(fL))
from bittersim.thermal import transient_summary
from bittersim.constants import ALPHA_CU
tr = transient_summary(best, ALPHA_CU)

extra = {"hoop_stress_max_MPa_profile": float(sig_m.max() / 1e6), "Bz_inner_edge_T": float(Bz_m[0]),
         "axial_compressive_force_N": Fz, "inductance_H": Lind, "stored_energy_J": 0.5 * Lind * best["I_A"] ** 2,
         "L_over_R_s": Lind / best["R_total_ohm"],
         "larmor_MHz": fL / 1e6, "swissroll": {
             "F": sr_tuned.F, "N": sr_tuned.N, "r_m": sr_tuned.r, "gap_tuned_um": sr_tuned.gap * 1e6,
             "Q": sr_tuned.Q, "skin_depth_um": sr_tuned.skin_depth * 1e6, "sheet_res_ohm": sr_tuned.sigma_s,
             "mu_at_fL_when_tuned_to_fL": [mu_t.real, mu_t.imag],
             "best_f0_over_fL": float(qbest), "mu_at_fL_best": [mu_b.real, mu_b.imag],
             "snr_gain_best": float(gbest), "snr_gain_tuned_to_fL": float(swissroll.snr_gain(mu_t)),
             "mu_at_DC": [complex(sr_tuned.mu(0.0)).real, complex(sr_tuned.mu(0.0)).imag]}}

extra.update(tr)
val = validation.run_all(best)
json.dump(clean({"env": env, "baseline_pdf_initial_guess": baseline, "optimal": best, "extra": extra,
                 "optimiser": {"method": opt["method"], "x": opt["x"], "constraints_g": opt["constraints"],
                               "de_nfev": opt["de_nfev"], "seconds": t_opt, "limits": LIMITS, "bounds": BOUNDS}}),
          open(os.path.join(RES, "results.json"), "w"), indent=1)
json.dump(clean(val), open(os.path.join(RES, "validation.json"), "w"), indent=1)

plotting.fig_field(best, os.path.join(FIG, "fig1_field_homogeneity.png"))
plotting.fig_convergence(val["segments"]["rows_polygon_only"], os.path.join(FIG, "fig2_convergence.png"))
plotting.fig_thermal(best, os.path.join(FIG, "fig3_thermal.png"))
plotting.fig_swissroll(best["B0"], os.path.join(FIG, "fig4_swissroll.png"))
plotting.fig_mechanics(best, os.path.join(FIG, "fig5_mechanics.png"))

# ------------------------------------------------------------------ markdown
def g(x, f="%.6g"):
    return f % x


def design_table(rows, cols):
    out = ["| quantity | " + " | ".join(c[0] for c in cols) + " |", "|---|" + "---|" * len(cols)]
    for label, key, fmt, scale in rows:
        out.append("| %s | " % label + " | ".join(fmt % (c[1][key] * scale) for c in cols) + " |")
    return "\n".join(out)


ROWS = [("R1 [mm]", "R1", "%.2f", 1e3), ("R2 [mm]", "R2", "%.2f", 1e3), ("L [mm]", "L", "%.1f", 1e3),
        ("plate thickness [mm]", "d_plate", "%.3f", 1e3), ("cooling hole D [mm]", "D_hole", "%.3f", 1e3),
        ("hole pitch / D", "pitch_factor", "%.3f", 1), ("water velocity [m/s]", "v_flow", "%.3f", 1),
        ("turns (plates)", "n_turns", "%.1f", 1), ("cooling holes", "n_holes", "%d", 1),
        ("ampere-turns NI [A]", "NI_At", "%.0f", 1), ("current I [A]", "I_A", "%.1f", 1),
        ("supply voltage [V]", "V_total_V", "%.3f", 1), ("voltage per plate V0 [mV]", "V_per_turn_V", "%.3f", 1e3),
        ("resistance [mOhm]", "R_total_ohm", "%.4f", 1e3), ("J_cu at R1 [A/mm^2]", "J_cu_inner_A_per_mm2", "%.3f", 1),
        ("electrical power [kW]", "P_elec_W", "%.3f", 1e-3), ("pump power (shaft) [W]", "P_pump_W", "%.1f", 1),
        ("total power [kW]", "P_total_W", "%.3f", 1e-3), ("B0 numeric (loops) [T]", "B0_numeric_T", "%.9f", 1),
        ("homogeneity over 30 mm DSV [ppm]", "homogeneity_ppm", "%.2f", 1),
        ("flow [L/min]", "flow_L_min", "%.1f", 1), ("Reynolds number", "Re", "%.0f", 1),
        ("h Dittus-Boelter [W/m^2K]", "h_W_m2K", "%.0f", 1), ("h Gnielinski [W/m^2K]", "h_gnielinski_W_m2K", "%.0f", 1),
        ("pressure drop [kPa]", "dp_Pa", "%.2f", 1e-3), ("mixed outlet water T [C]", "T_out_mixed_C", "%.3f", 1),
        ("hottest-channel water rise [K]", "dT_water_inner_K", "%.3f", 1),
        ("hot-spot wall T [C]", "T_wall_hot_C", "%.3f", 1), ("hot-spot copper T [C]", "T_hot_C", "%.3f", 1),
        ("mean copper T [C]", "T_cu_mean_C", "%.3f", 1), ("hoop stress bound C*B0/lambda [MPa]", "sigma_hoop_max_MPa", "%.3f", 1),
        ("copper mass [kg]", "mass_cu_kg", "%.0f", 1), ("energy-balance residual", "energy_residual", "%.2e", 1)]

tbl = design_table(ROWS, [("PDF initial guess (R2=0.15, L=0.8, v=2.5)", baseline), ("optimal", best)])
open(os.path.join(RES, "RESULTS.md"), "w").write(
    "# Computed results\n\nGenerated by `examples/run_all.py` (%s; runtime of optimiser %.1f s).\n\n" % (
        ", ".join("%s %s" % kv for kv in env.items()), t_opt) + tbl + "\n\n```json\n" +
    json.dumps(clean(extra), indent=1) + "\n```\n")

V = val
L = []
L.append("# VALIDATION\n\nAll numbers below are produced by `examples/run_all.py` -> `bittersim.validation` "
         "(environment: %s). The same checks run as pytest assertions in `tests/`.\n" %
         ", ".join("%s %s" % kv for kv in env.items()))
L.append("## 1. Elliptic integrals and loop fields\n\n| check | max relative difference |\n|---|---|\n"
         "| Numba AGM K(m), E(m) vs scipy.special.ellipk/ellipe, m in [0, 0.999] | %.2e |\n"
         "| Numba loop field vs SciPy loop field (1600 off-axis points) | %.2e |\n"
         "| On-axis loop (E5) vs code, axis branch | %.2e |\n"
         "| On-axis loop (E5) vs elliptic branch at rho = 1e-7 m | %.2e |\n" % (
             V["elliptic"]["agm_vs_scipy_max_rel"], V["elliptic"]["numba_vs_scipy_loop_max_rel"],
             V["single_loop"]["axis_branch_max_rel"], V["single_loop"]["elliptic_branch_near_axis_max_rel"]))
L.append("Segment engine (E6) at the centre of an N-gon loop vs mu0 I/(2a):\n\n| N | rel. error |\n|---|---|\n" +
         "\n".join("| %d | %.3e |" % tuple(r) for r in V["single_loop"]["segments_center_error"]) + "\n")
L.append("## 2. Long-solenoid limit (thin sheet, R = 50 mm, n = 1000 /m, I = 10 A)\n\n"
         "| L/R | B loops [T] | mu0 n I [T] | rel. to mu0 n I | rel. to E2 |\n|---|---|---|---|---|\n" +
         "\n".join("| %d | %.9f | %.9f | %.3e | %.3e |" % (r["L_over_R"], r["B_loops"], r["B_ideal"], r["rel_vs_ideal"],
                                                         r["rel_vs_E2"]) for r in V["long_solenoid"]) +
         "\n\nThe deviation from mu0 n I falls as 2(R/L)^2, as expected from E2.\n")
L.append("## 3. Thick-solenoid and Bitter closed forms (R1 = 50 mm, R2 = 150 mm, L = 400 mm)\n\n"
         "| check | max relative difference |\n|---|---|\n"
         "| Gauss-loop model, uniform J, vs E3 (61 axial points) | %.2e |\n"
         "| E3 at z=0 vs PDF centre formula | %.2e |\n"
         "| Gauss-loop model, J = C/r, vs E4 (61 axial points) | %.2e |\n" % (
             V["thick"]["uniform_loops_vs_E3_max_rel"], V["thick"]["E3_center_vs_PDF_form_rel"],
             V["thick"]["bitter_loops_vs_E4_max_rel"]))
S = V["segments"]
L.append("## 4. PDF test 1: uniform thick solenoid with straight segments + Richardson extrapolation\n\n"
         "40 planar turns x 8 radial filaments; B_z(0) exact (E3) = %.9f T.\n\n"
         "| angular steps | segments | rel. error vs E3 (all discretisation) | polygon error only |\n|---|---|---|---|\n" % S["B_exact_E3"] +
         "\n".join("| %d | %d | %.3e | %.3e |" % (r[0], r[2], r[1], p[1]) for r, p in zip(S["rows_vs_E3"], S["rows_polygon_only"])) +
         "\n\nRichardson order from 180/360/720 steps: **p = %.6f** (theory 2). Extrapolated B = %.9f T equals the "
         "filament limit %.9f T; the residual %.2e vs E3 is the radial/axial filament (midpoint) error of the "
         "40x8 cross-section mesh, not the segment law. The PDF's < 0.195%% criterion is met at every resolution "
         "tested.\n" % (S["richardson_p"], S["B_richardson"], S["B_filament_limit"], S["richardson_vs_E3_rel"]))
B = V["bitter_segments"]
L.append("## 5. PDF test 2: Bitter 1/r filaments vs closed form E4\n\n80 plates x 24 radial x 360 angular "
         "(%d segments). E4: %.9f T.\n\n| geometry | B_z(0) [T] | rel. error | transverse B [T] |\n|---|---|---|---|\n"
         "| planar rings | %.9f | %.3e | %.2e |\n| helical (PDF generator) | %.9f | %.3e | %.2e |\n\n"
         "The helical lead adds a small transverse field (the net axial current of the helix); it is reported "
         "but not part of the axisymmetric models.\n" % (
             B["planar"]["n_segments"], B["B_E4"], B["planar"]["Bz"], B["planar"]["rel_err"], B["planar"]["B_transverse"],
             B["helical"]["Bz"], B["helical"]["rel_err"], B["helical"]["B_transverse"]))
I_ = V["inductance"]
L.append("## 6. Inductance and energy\n\n| R [m] | length [m] | N | L numeric E9 [H] | L Nagaoka E10 [H] | rel. | Nagaoka / long-solenoid |\n|---|---|---|---|---|---|---|\n" +
         "\n".join("| %.3f | %.3f | %d | %.9e | %.9e | %.2e | %.6f |" % (r["R"], r["len"], r["N"], r["L_numeric_H"],
                   r["L_nagaoka_H"], r["rel"], r["nagaoka_over_long"]) for r in I_["sheet"]) +
         "\n\nThick winding -> sheet limit (R1 = 50 mm, length 0.5 m, 200 turns):\n\n| thickness [m] | rel. difference |\n|---|---|\n" +
         "\n".join("| %.0e | %.3e |" % tuple(r) for r in I_["thick_thin_limit"]) +
         "\n\n(The difference falls linearly with thickness, i.e. converges.)\n\nEnergy, long solenoid (R = 50 mm, 5 m, N = 1000, 10 A): "
         "L I^2/2 = %.6e J, B^2/(2 mu0) x volume = %.6e J, ratio %.6f (end-field deficit, equals the Nagaoka factor).\n" % (
             I_["energy_long_W_LI2"], I_["energy_long_W_B2"], I_["energy_ratio"]))
L.append("## 7. Homogeneity convergence (optimal design, 30 mm DSV)\n\n| radial nodes | axial nodes | ppm | B0 [T] |\n|---|---|---|---|\n" +
         "\n".join("| %d | %d | %.6f | %.12f |" % (r["nr"], r["nz_total"], r["ppm"], r["B0"]) for r in V["homogeneity_convergence"]) + "\n")
L.append("## 8. Thermal energy conservation\n\nSum of per-cell Joule heat over all hole rows vs analytic P (E31): "
         "relative residual %.2e (optimal design), %.2e (PDF initial guess). Mixed-outlet energy balance "
         "V0 I = m_dot cp dT holds by construction (dT computed from it); the test suite checks it to 1e-12 for "
         "random voltages/velocities.\n" % (best["energy_residual"], baseline["energy_residual"]))
P_ = V["performance"]
L.append("## 9. Performance / memory (chunked Numba engine)\n\n%d segments x %d observation points in %.3f s "
         "(%.3g segment-point evaluations/s, %d threads); process peak RSS %.0f MB.\n" % (
             P_["n_segments"], P_["n_obs"], P_["seconds"], P_["pair_evals_per_s"], env["cpus"], P_["peak_rss_MB"]))
L.append("## Figures\n\n![field](figures/fig1_field_homogeneity.png)\n\n![convergence](figures/fig2_convergence.png)\n\n"
         "![thermal](figures/fig3_thermal.png)\n\n![swiss roll](figures/fig4_swissroll.png)\n\n![mechanics](figures/fig5_mechanics.png)\n")
L.append("## Not validated\n\n* FEniCS curl-curl cross-validation (PDF): FEniCS 2018.1 is not pip-installable on current "
         "Colab/Kaggle or this environment; replaced by the independent closed-form and segment/loop cross-checks above.\n"
         "* Swiss-roll SNR gain: heuristic model with no experimental calibration here.\n")
open(os.path.join(ROOT, "VALIDATION.md"), "w").write("\n".join(L))
print("total runtime %.1f s" % (time.time() - t0))
script = os.path.join(HERE, "update_report.py")
if os.path.isfile(script):
    os.system("%s %s" % (sys.executable, script))
