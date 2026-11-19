"""Quick reproduction check: re-derive the published numbers and redraw figures 1-5.

    python examples/reproduce.py                 # check numbers, redraw figures into repro_out/
    python examples/reproduce.py --out DIR       # write the figures somewhere else
    python examples/reproduce.py --no-figures    # numbers only

The published design vector (results/results.json -> optimiser.x) is re-evaluated with the
current code. The optimiser is NOT rerun: the published optimum predates the switch to a soft
Reynolds-number penalty, so a fresh optimiser run lands on a different design (see README,
"Reproducing the results").

The quantities behind README rows R1-R51 are recomputed with the same calls that produced them
(examples/run_all.py for R1-R40, examples/run_fmri.py for R41-R51) and compared with the values
stored in results/results.json, results/validation.json and results/fmri.json. Nothing under
results/, figures/, docs/ or report/ is written.

Exit status: 0 if every quantity matches within its tolerance, 1 otherwise.
"""
from __future__ import print_function

import argparse
import json
import os
import sys
import time
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "package"))
os.environ.setdefault("MPLBACKEND", "Agg")
warnings.filterwarnings("ignore")

RTOL = 1e-6   # default relative tolerance; deterministic double-precision code, so this is loose


def load(name):
    with open(os.path.join(ROOT, "results", name)) as fh:
        return json.load(fh)


class Checker(object):
    """Collects (row, quantity, reproduced, stored) comparisons."""

    def __init__(self):
        self.rows = []

    def num(self, row, name, new, ref, rtol=RTOL, atol=0.0):
        new, ref = float(new), float(ref)
        err = abs(new - ref)
        ok = err <= atol + rtol * abs(ref)
        rel = err / abs(ref) if ref != 0.0 else (0.0 if err == 0.0 else float("inf"))
        self.rows.append((row, name, "%.10g" % new, "%.10g" % ref, rel, ok))

    def cplx(self, row, name, new, ref_pair, **kw):
        self.num(row, name + " (real)", complex(new).real, ref_pair[0], **kw)
        self.num(row, name + " (imag)", complex(new).imag, ref_pair[1], **kw)

    def same(self, row, name, new, ref):
        self.rows.append((row, name, repr(new), repr(ref), 0.0 if new == ref else float("inf"), new == ref))

    def report(self):
        print("%-4s %-58s %18s %18s %9s  %s" % ("row", "quantity", "reproduced", "stored", "rel.diff", "status"))
        for row, name, new, ref, rel, ok in self.rows:
            print("%-4s %-58s %18s %18s %9.1e  %s" % (row, name[:58], new, ref, rel, "ok" if ok else "MISMATCH"))
        bad = [r for r in self.rows if not r[5]]
        worst = max(r[4] for r in self.rows if r[5]) if len(bad) < len(self.rows) else float("nan")
        return bad, worst


def check_design(chk):
    """R1-R34: published optimum and PDF initial guess (as in examples/run_all.py)."""
    from bittersim import evaluate_design, BitterDesign, swissroll, mechanics, inductance
    from bittersim.optimize import x_to_design, constraints, correlation_penalty
    from bittersim.thermal import transient_summary
    from bittersim.constants import ALPHA_CU

    D = load("results.json")
    o, b, x = D["optimal"], D["baseline_pdf_initial_guess"], D["extra"]
    sr_ref = x["swissroll"]
    best = evaluate_design(x_to_design(D["optimiser"]["x"]))
    base = evaluate_design(BitterDesign())

    for row, keys in (("R1", ["R1"]), ("R2", ["R2"]), ("R3", ["L"]), ("R4", ["d_plate", "d_ins"]),
                      ("R5", ["n_turns"]), ("R6", ["D_hole", "pitch_factor"]), ("R7", ["NI_At"]),
                      ("R8", ["I_A"]), ("R9", ["V_total_V", "V_per_turn_V"]), ("R10", ["R_total_ohm"]),
                      ("R11", ["J_cu_inner_A_per_mm2", "J_cu_outer_A_per_mm2"]), ("R12", ["P_elec_W"]),
                      ("R13", ["P_pump_W"]), ("R14", ["B0_numeric_T"]), ("R15", ["homogeneity_ppm"]),
                      ("R16", ["v_flow", "flow_L_min", "m_dot_kg_s"]), ("R17", ["Re", "Pr"]),
                      ("R18", ["h_W_m2K", "h_gnielinski_W_m2K"]), ("R19", ["f_darcy", "dp_Pa"]),
                      ("R20", ["T_out_mixed_C"]),
                      ("R21", ["dT_water_inner_K", "dT_film_inner_K", "dT_cond_inner_K"]),
                      ("R22", ["T_hot_C"]), ("R23", ["T_cu_mean_C"]), ("R24", ["sigma_hoop_max_MPa"]),
                      ("R28", ["mass_cu_kg"])):
        for k in keys:
            chk.num(row, "optimum: " + k, best[k], o[k])
    chk.same("R6", "optimum: n_holes / n_rows", (best["n_holes"], best["n_rows"]), (o["n_holes"], o["n_rows"]))
    chk.num("R34", "optimum: energy-balance residual (round-off)", best["energy_residual"], o["energy_residual"],
            rtol=0.0, atol=1e-12)

    _, _, sig = mechanics.hoop_stress_profile(best["R1"], best["R2"], best["L"], best["C_A_per_m"], best["fill_total"])
    chk.num("R24", "max hoop stress, field profile [MPa]", sig.max() / 1e6, x["hoop_stress_max_MPa_profile"])
    chk.num("R25", "axial force on each half [N]",
            mechanics.axial_force(best["R1"], best["R2"], best["L"], best["C_A_per_m"]), x["axial_compressive_force_N"])
    L_ind = inductance.coil_inductance(best["R1"], best["R2"], best["L"], best["n_turns"], "bitter", nr=10)
    chk.num("R26", "inductance [H]", L_ind, x["inductance_H"])
    chk.num("R26", "stored energy [J]", 0.5 * L_ind * best["I_A"] ** 2, x["stored_energy_J"])
    chk.num("R26", "L/R time constant [s]", L_ind / best["R_total_ohm"], x["L_over_R_s"])

    tr = transient_summary(best, ALPHA_CU)
    chk.num("R27", "lumped thermal time constant (63 %) [s]", tr["transient_tau63_s"], x["transient_tau63_s"])
    # first point of a 1 s output grid at or above 85 C: allow one grid step
    chk.num("R27", "pump failure: lumped time to 85 C [s]", tr["pump_failure_time_to_85C_s"],
            x["pump_failure_time_to_85C_s"], atol=1.0)
    chk.num("R27", "pump failure: hot-spot adiabatic bound [s]", tr["pump_failure_hotspot_adiabatic_s"],
            x["pump_failure_hotspot_adiabatic_s"])

    fL = swissroll.larmor_hz(best["B0"])
    chk.num("R29", "Larmor frequency [MHz]", fL / 1e6, x["larmor_MHz"])
    tuned = swissroll.SwissRoll(fL)
    chk.cplx("R30", "Swiss roll tuned to f_L: mu_eff(f_L)", tuned.mu(fL), sr_ref["mu_at_fL_when_tuned_to_fL"])
    chk.num("R30", "Swiss roll Q", tuned.Q, sr_ref["Q"])
    _, _, q_best, g_best = swissroll.best_detuning(best["B0"])
    # argmax over a 3001-point grid: allow one grid step
    chk.num("R31", "best roll tuning f0/f_L", q_best, sr_ref["best_f0_over_fL"], atol=0.15 / 3000)
    chk.cplx("R31", "mu_eff(f_L) at best tuning", swissroll.SwissRoll(fL * q_best).mu(fL), sr_ref["mu_at_fL_best"])
    chk.num("R31", "legacy heuristic SNR gain (E28)", g_best, sr_ref["snr_gain_best"])
    chk.cplx("R32", "Swiss roll mu_eff at DC", tuned.mu(0.0), sr_ref["mu_at_DC"], atol=1e-15)

    for k in ("P_elec_W", "V_total_V", "homogeneity_ppm", "T_hot_C"):
        chk.num("R33", "PDF initial guess: " + k, base[k], b[k])

    g = constraints(best)
    chk.same("-", "optimum feasible under current hard constraints", bool((g >= -1e-9).all()), True)
    chk.same("-", "soft Re penalty at the optimum (Re >= 5500)", correlation_penalty(best), 0.0)
    return best, tr, L_ind / best["R_total_ohm"]


def check_validation(chk):
    """R35-R40 (as in bittersim.validation, used by examples/run_all.py)."""
    from bittersim import validation
    VJ = load("validation.json")
    seg = validation.v_segment_convergence()
    chk.num("R35", "Richardson order, segment engine", seg["richardson_p"], VJ["segments"]["richardson_p"])
    chk.num("R36", "segment engine vs E3, 360 angular steps", seg["rows_vs_E3"][3][1], VJ["segments"]["rows_vs_E3"][3][1])
    chk.num("R37", "Bitter helical filaments vs E4",
            validation.v_bitter_segments()["helical"]["rel_err"], VJ["bitter_segments"]["helical"]["rel_err"])
    th = validation.v_thick_closed_form()
    # these two are round-off-level agreements (~1e-13): check they stay at round-off level
    chk.num("R38", "loop model vs E3 (round-off level)", th["uniform_loops_vs_E3_max_rel"],
            VJ["thick"]["uniform_loops_vs_E3_max_rel"], rtol=0.0, atol=1e-12)
    chk.num("R38", "loop model vs E4 (round-off level)", th["bitter_loops_vs_E4_max_rel"],
            VJ["thick"]["bitter_loops_vs_E4_max_rel"], rtol=0.0, atol=1e-12)
    chk.num("R39", "long solenoid L/R = 1000 vs mu0 n I", validation.v_long_solenoid()[2]["rel_vs_ideal"],
            VJ["long_solenoid"][2]["rel_vs_ideal"])
    # quadrature-tolerance level (~1e-11)
    chk.num("R40", "sheet inductance vs Nagaoka (worst of 3)",
            max(r["rel"] for r in validation.v_inductance()["sheet"]),
            max(r["rel"] for r in VJ["inductance"]["sheet"]), rtol=0.0, atol=1e-9)
    return seg


def check_fmri(chk, best, tau, L_over_R):
    """R41-R51 (as in examples/run_fmri.py)."""
    from bittersim import harmonics as H, stability as S, rfsnr as Q
    F = load("fmri.json")
    geo = (best["R1"], best["R2"], best["L"], best["C_A_per_m"])

    s30, s40 = H.shimmed_homogeneity(*geo, dsv=0.03), H.shimmed_homogeneity(*geo, dsv=0.04)
    r30, r40 = F["shim_30mm"], F["shim_40mm"]
    chk.num("R41", "zonal Z2 over 30 mm DSV [ppm]", s30["zonal_ppm_unshimmed"][2], r30["zonal_ppm_unshimmed"][2])
    chk.num("R41", "zonal Z4 over 30 mm DSV [ppm]", s30["zonal_ppm_unshimmed"][4], r30["zonal_ppm_unshimmed"][4])
    for lab, new, ref in (("30", s30, r30), ("40", s40, r40)):
        chk.num("R42", "p-p over %s mm DSV, unshimmed [ppm]" % lab, new["ppm_unshimmed"], ref["ppm_unshimmed"])
        chk.num("R42", "p-p over %s mm DSV, Z2/Z4 shim pairs [ppm]" % lab, new["ppm_shimmed"], ref["ppm_shimmed"])
    chk.num("R43", "shim radius [m]", s30["shim_radius_m"], r30["shim_radius_m"])
    for i in (0, 1):
        chk.num("R43", "shim pair %d: z [m]" % (i + 1), s30["shim_z_m"][i], r30["shim_z_m"][i])
        chk.num("R43", "shim pair %d: NI [A]" % (i + 1), s30["shim_NI_A"][i], r30["shim_NI_A"][i])
    chk.num("R43", "shim power [W]", s30["shim_power_W"], r30["shim_power_W"])

    tol = H.tolerance_harmonics(*geo, dsv=0.03)["tesseral_ppm"]
    for k in ("A11", "B21"):
        chk.num("R44", "tolerance model, tesseral %s [ppm]" % k, tol[k], F["tolerance_30mm"]["tesseral_ppm"][k])

    hp, hr = H.head_preset(), F["head_preset"]
    for k in ("R2", "L", "P_elec_W", "P_pump_W", "V_total_V", "I_A", "mass_cu_kg", "flow_L_min",
              "ppm_unshimmed", "ppm_shimmed", "T_hot_C"):
        chk.num("R45", "head preset (grid): " + k, hp[k], hr[k])
    chk.same("R45", "head preset: violates 8 V supply", hp["violates_8V_supply"], hr["violates_8V_supply"])

    chk.num("R46", "field drift vs copper T at fixed current [1/K]",
            S.expansion_coefficient_numeric(best["R1"], best["R2"], best["L"], best["NI_At"]), F["alpha_B_per_K_numeric"])
    cur = S.simulate(best, tau, L_over_R, {"mode": "current"})
    vol = S.simulate(best, tau, L_over_R, {"mode": "voltage"})
    for k in ("psu_ripple", "psu_drift", "thermal_expansion", "total"):
        chk.num("R47", "current-regulated budget: %s [ppm]" % k, cur["budget_ppm"][k], F["stability"]["current"]["budget_ppm"][k])
    chk.num("R47", "current-regulated total [Hz]", cur["budget_Hz"]["total"], F["stability"]["current"]["budget_Hz"]["total"])
    chk.num("R47", "voltage-regulated total [ppm]", vol["budget_ppm"]["total"], F["stability"]["voltage"]["budget_ppm"]["total"])
    req = S.requirements(best, tau)
    for k in ("max_copper_dT_K_current_mode", "max_copper_dT_K_voltage_mode",
              "max_water_amp_K_current_mode", "max_water_amp_K_voltage_mode"):
        chk.num("R48", k, req[k], F["stability_requirements"][k])

    rf, rr = Q.compare(B0=best["B0"]), F["rf"]
    chk.num("R49", "RF SNR, slab vs same coil over air", rf["gain_vs_air"], rr["gain_vs_air"])
    chk.num("R49", "RF SNR, slab vs coil on tissue", rf["gain_vs_contact"], rr["gain_vs_contact"])
    chk.num("R49", "best roll tuning f0/f_L (RF model)", rf["detune"], rr["detune"])
    chk.cplx("R49", "mu_eff at best tuning", complex(*rf["mu"]), rr["mu"])
    for lm in (1.0, 10.0):
        chk.num("R50", "RF SNR vs air, loss multiplier %g" % lm, Q.compare(B0=best["B0"], loss_multiplier=lm)["gain_vs_air"],
                rr["gain_vs_air_by_loss_multiplier"][str(lm)])
    for k in ("R_coil", "R_tissue", "R_slab"):
        chk.num("R51", "RF resistance with slab: %s [ohm]" % k, rf["slab"][k], rr["slab"][k])


def draw_figures(best, seg, out):
    """Figures 1-5 of the README, drawn by the same bittersim.plotting calls as examples/run_all.py."""
    import matplotlib.pyplot as plt
    from bittersim import plotting
    if not os.path.isdir(out):
        os.makedirs(out)
    paths = []

    def save(name, fig):
        path = os.path.join(out, name)
        fig.savefig(path, dpi=130)
        plt.close(fig)
        paths.append(path)

    save("fig1_field_homogeneity.png", plotting.fig_field(best))
    save("fig2_convergence.png", plotting.fig_convergence(seg["rows_polygon_only"]))
    save("fig3_thermal.png", plotting.fig_thermal(best))
    save("fig4_swissroll.png", plotting.fig_swissroll(best["B0"]))
    save("fig5_mechanics.png", plotting.fig_mechanics(best)[0])
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.join(ROOT, "repro_out"), help="figure directory (default: repro_out/)")
    ap.add_argument("--no-figures", action="store_true", help="only check the numbers")
    a = ap.parse_args(argv)

    t0 = time.time()
    import numpy, scipy, numba, matplotlib
    print("environment: python %s, numpy %s, scipy %s, numba %s, matplotlib %s" % (
        sys.version.split()[0], numpy.__version__, scipy.__version__, numba.__version__, matplotlib.__version__))
    chk = Checker()
    best, tr, L_over_R = check_design(chk)
    seg = check_validation(chk)
    check_fmri(chk, best, tr["transient_tau_analytic_s"], L_over_R)
    print()
    bad, worst = chk.report()
    print()
    if bad:
        print("%d of %d quantities do NOT match the stored results." % (len(bad), len(chk.rows)))
    else:
        print("All %d quantities match the stored results (largest relative difference %.1e)." % (len(chk.rows), worst))
    if not a.no_figures:
        for p in draw_figures(best, seg, a.out):
            rel = os.path.relpath(p)
            print("wrote %s" % (p if rel.startswith("..") else rel))
    print("done in %.1f s" % (time.time() - t0))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
