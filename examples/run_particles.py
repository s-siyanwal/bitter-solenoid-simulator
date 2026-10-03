"""Particle-level emulation vs continuum, for the seed-1 optimum in results/results.json.

Writes results/particles.json, figures/fig8_particles.png (comparison) and
figures/fig9_particle_convergence.png (error vs particle count). Every README number in the
"Particle emulation" section is read from results/particles.json by examples/make_docs.py.
Usage: python examples/run_particles.py [--quick | --replot]  (--replot redraws fig9 from the JSON)
"""
from __future__ import division
import json, math, os, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "package"))
from bittersim import particles as Pt                                    # noqa: E402
from bittersim import harmonics as H                                     # noqa: E402
from bittersim.fields import B_bitter_axis                               # noqa: E402
from bittersim.thermal import annulus_conduction_dT                      # noqa: E402
from bittersim.constants import RHO_CU_20, ALPHA_CU, K_CU                # noqa: E402
from bittersim.design import BitterDesign                                # noqa: E402
from bittersim.emulation import emulate                                  # noqa: E402

QUICK = "--quick" in sys.argv


def timed(f, *a, **k):
    t = time.time()
    r = f(*a, **k)
    return r, time.time() - t


def rho_lin(T):
    return RHO_CU_20 * (1.0 + ALPHA_CU * (T - 20.0))


def ppm_of(b):
    B = np.linalg.norm(b, axis=1)
    return float((B.max() - B.min()) / B[-1] * 1e6)


def plot_convergence(out):
    """figures/fig9_particle_convergence.png from the convergence block of particles.json."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    conv = out["convergence"]
    fig, axc = plt.subplots(figsize=(7, 5))
    F = conv["field"]
    N = np.array([c["N"] for c in F], float)
    axc.loglog(N, [c["axis_rel_max_inner"] for c in F], "o-", label="axis |z|<=0.4L, max rel err")
    axc.loglog(N, [abs(c["ppm"] - out["production"]["ppm_loop_same_points"]) / 1e6 for c in F], "s-", label="DSV ppm error / 1e6")
    M = conv["moments"]
    axc.loglog([c["N"] for c in M], [max(c["P_rel"], 1e-17) for c in M], "^-", label="NI, P rel err")
    G = conv["carriers"]
    axc.loglog([c["N"] for c in G], [c["rho_rel"] for c in G], "d-", label="carrier rho rel err")
    W = conv["walkers"]
    axc.loglog([c["N"] for c in W], [c["dT_rel"] for c in W], "v-", label="walker dT rel err")
    nn = np.logspace(2, math.log10(4.0 * max(c["N"] for c in M)), 20)
    axc.loglog(nn, 1.0 / np.sqrt(nn), "k:", label="1/sqrt(N)")
    axc.loglog(nn, 1.0 / nn, "k--", label="1/N")
    axc.set_xlabel("particles N"); axc.set_ylabel("relative error vs continuum"); axc.legend(fontsize=8)
    axc.set_title("Particle emulation convergence")
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "fig9_particle_convergence.png"), dpi=110)


def main():
    import numba
    D = json.load(open(os.path.join(ROOT, "results", "results.json")))
    o = D["optimal"]
    R1, R2, L, C = o["R1"], o["R2"], o["L"], o["C_A_per_m"]
    design = BitterDesign(**dict((k, o[k]) for k in ("R1", "R2", "L", "d_plate", "d_ins", "D_hole", "v_flow",
                                                     "pitch_factor", "T_in", "B0", "dsv", "K_minor", "eta_pump")))
    lay = design.hole_layout()
    lam = o["fill_total"]
    n_turns = o["n_turns"]
    a_h = design.D_hole / 2.0
    cell_area = lay["dr_row"] * 2 * math.pi * lay["r_rows"][0] / lay["n_per_row"][0]
    b_h = math.sqrt(cell_area / math.pi)
    out = {"design": "seed-1 optimum (results/results.json)", "threads": int(getattr(numba, "get_num_threads", lambda: numba.config.NUMBA_NUM_THREADS)()),
           "cpu_count": os.cpu_count(), "quick": QUICK}

    # ---------------- continuum references ----------------
    loops = H.coil_loopset(R1, R2, L, C)
    ax, dsv = Pt.axis_and_dsv_points(L, o["dsv"])
    zin = np.abs(ax[:, 2]) <= 0.4 * L + 1e-12
    B_ax_E = B_bitter_axis(R1, R2, L, C, ax[:, 2])
    ppm_loop_same_pts = ppm_of(loops.field_xyz(dsv))

    # ---------------- warm-up (compile both builds) ----------------
    for par in (False, True):
        Pt.element_field(dsv[:2], R1, R2, L, C, 64, par)
        Pt.element_moments(R1, R2, L, C, 64, RHO_CU_20, lam, par)
        Pt.green_kubo_rho(20.0, 64, 1, par)
        Pt.walk_dT(1.0, a_h, b_h, K_CU, 64, 2, par)

    # ---------------- serial vs parallel benchmark + identity ----------------
    nb = {"field": 1 << (11 if QUICK else 14), "moments": 1 << 22, "carriers": 20000 if QUICK else 200000,
          "walkers": 4000 if QUICK else 40000}
    bench = {}
    jobs = {
        "field": lambda p: Pt.element_field(dsv, R1, R2, L, C, nb["field"], p),
        "moments": lambda p: np.array(sorted(Pt.element_moments(R1, R2, L, C, nb["moments"], RHO_CU_20, lam, p).items()))[:, 1].astype(float),
        "carriers": lambda p: np.array(Pt.green_kubo_rho(o["T_cu_mean_C"], nb["carriers"], 11, p)),
        "walkers": lambda p: np.array(Pt.walk_dT(1.0, a_h, b_h, K_CU, nb["walkers"], 12, p)),
    }
    for name, f in jobs.items():
        rs, ts = timed(f, False)
        rp, tp = timed(f, True)
        bench[name] = {"n": nb[name], "serial_s": ts, "parallel_s": tp, "speedup": ts / tp,
                       "bit_identical": bool(np.array_equal(np.asarray(rs), np.asarray(rp)))}
    nmc = 200 if QUICK else 2000
    em1, t1 = timed(emulate, design, realizations=nmc, seed=7, workers=1)
    em8, t8 = timed(emulate, design, realizations=nmc, seed=7, workers=8)
    bench["tolerance_mc_processes"] = {"n": nmc, "serial_s": t1, "parallel_s": t8, "speedup": t1 / t8,
                                       "bit_identical": json.dumps(em1, sort_keys=True, default=str) == json.dumps(em8, sort_keys=True, default=str)}
    out["benchmark"] = bench

    # ---------------- convergence sweeps ----------------
    conv = {"field": [], "moments": [], "carriers": [], "walkers": []}
    for e in range(10, (15 if QUICK else 19) + 1):
        N = 1 << e
        b, t = timed(Pt.element_field, np.vstack([ax, dsv]), R1, R2, L, C, N, True)
        bax, bd = b[:len(ax)], b[len(ax):]
        conv["field"].append({"N": N, "seconds": t,
                              "B0_rel": abs(bd[-1, 2] / o["B0_numeric_T"] - 1.0),
                              "axis_rel_max_inner": float(np.max(np.abs(bax[zin, 2] / B_ax_E[zin] - 1.0))),
                              "axis_rel_max_all": float(np.max(np.abs(bax[:, 2] / B_ax_E - 1.0))),
                              "ppm": ppm_of(bd)})
    NI_c = o["NI_At"]
    P_c = o["P_elec_W"]
    for e in range(8, 23, 2):
        m = Pt.element_moments(R1, R2, L, C, 1 << e, o["rho_cu_mean"], lam, True)
        conv["moments"].append({"N": 1 << e, "NI_rel": abs(m["NI"] / NI_c - 1), "P_rel": abs(m["P"] / P_c - 1)})
    for n in ([1000, 4000, 16000] if QUICK else [1000, 4000, 16000, 64000, 256000, 1000000]):
        r, rt, se = Pt.green_kubo_rho(o["T_cu_mean_C"], n, 21, True)
        conv["carriers"].append({"N": n, "rho_rel": abs(r / rt - 1), "stderr_rel": se})
    q_hot_c = rho_lin(o["T_hot_C"]) * (C / (lam * R1)) ** 2
    dT_E18 = annulus_conduction_dT(q_hot_c, a_h, b_h, K_CU)
    for n in ([100, 400, 1600] if QUICK else [100, 400, 1600, 6400, 25600, 100000]):
        w, se = Pt.walk_dT(q_hot_c, a_h, b_h, K_CU, n, 22, True)
        conv["walkers"].append({"N": n, "dT_rel": abs(w / dT_E18 - 1), "stderr_rel": se / dT_E18})
    out["convergence"] = conv

    # ---------------- production comparison ----------------
    f_last = conv["field"][-1]
    Nf = f_last["N"]
    b = Pt.element_field(np.vstack([ax, dsv]), R1, R2, L, C, Nf, True)
    bax, bd = b[:len(ax)], b[len(ax):]
    m = Pt.element_moments(R1, R2, L, C, 1 << 22, o["rho_cu_mean"], lam, True)
    I_p = m["NI"] / n_turns
    R_p = m["P"] / I_p ** 2
    ncar = 16000 if QUICK else 1000000
    rho_em, rho_bg, rho_se = Pt.green_kubo_rho(o["T_cu_mean_C"], ncar, 31, True)
    # Bloch-Gruneisen slope at 20 C (model property, finite difference)
    alpha_bg = (Pt.bloch_gruneisen_rho(20.5) - Pt.bloch_gruneisen_rho(19.5)) / Pt.RHO_CU_20 if hasattr(Pt, "RHO_CU_20") else \
        (Pt.bloch_gruneisen_rho(20.5) - Pt.bloch_gruneisen_rho(19.5)) / RHO_CU_20
    nw = 1600 if QUICK else 100000
    tau_fac, tau_se = Pt.walk_dT(1.0, a_h, b_h, K_CU, nw, 41, True)           # dT per unit q_v
    dT_walk = q_hot_c * tau_fac
    # emulated hot spot: water rise and film from the continuum correlations (scaled by the local
    # resistivity), conduction from the walkers, resistivity from the carrier model; fixed point in T
    wf = (o["dT_water_inner_K"] + o["dT_film_inner_K"]) / rho_lin(o["T_hot_C"])
    T = o["T_hot_C"]
    for _ in range(100):
        rh = Pt.bloch_gruneisen_rho(T)
        T_new = o["T_in"] + wf * rh + rh * (C / (lam * R1)) ** 2 * tau_fac
        if abs(T_new - T) < 1e-12:
            break
        T = T_new
    T_hot_em = T_new
    rho_hot_em, _, rho_hot_se = Pt.green_kubo_rho(T_hot_em, ncar, 32, True)
    T_hot_em_gk = o["T_in"] + (wf + (C / (lam * R1)) ** 2 * tau_fac) * rho_hot_em

    rows = []
    def row(key, label, unit, cont, emu, se_rel, expect):
        rows.append({"key": key, "label": label, "unit": unit, "continuum": cont, "particle": emu,
                     "rel_err": (emu - cont) / cont, "stderr_rel": se_rel, "expect": expect})
    row("B0", "Isocentre field B0", "T", o["B0_numeric_T"], float(bd[-1, 2]), 0.0, "agree")
    iz = int(np.argmin(np.abs(ax[:, 2] - 0.25 * L)))
    row("Bz_quarter", "Bz on axis at z = %.3f m (L/4)" % ax[iz, 2], "T", float(B_ax_E[iz]), float(bax[iz, 2]), None, "agree")
    iend = len(ax) - 1
    row("Bz_end", "Bz on axis at coil end z = L/2", "T", float(B_ax_E[iend]), float(bax[iend, 2]), None, "agree (slow MC)")
    row("ppm", "Homogeneity over 30 mm DSV", "ppm", ppm_loop_same_pts, ppm_of(bd), None, "agree (needs ~1e-7 rel. accuracy)")
    row("NI", "Ampere-turns NI", "A", NI_c, m["NI"], None, "agree")
    row("I", "Current I = NI / n_turns", "A", o["I_A"], I_p, None, "agree")
    row("R", "Resistance (rho at mean Cu T)", "ohm", o["R_total_ohm"], R_p, None, "agree")
    row("V", "Supply voltage", "V", o["V_total_V"], I_p * R_p, None, "agree")
    row("P", "Joule power (rho at mean Cu T)", "W", P_c, m["P"], None, "agree")
    row("rho_mean", "Resistivity at mean Cu T (%.3f C), carriers + Bloch-Gruneisen" % o["T_cu_mean_C"], "ohm m",
        o["rho_cu_mean"], rho_em, rho_se, "differ (model)")
    row("rho_bg_model", "Bloch-Gruneisen rho at mean Cu T (model value, no sampling noise)", "ohm m",
        o["rho_cu_mean"], rho_bg, 0.0, "differ (model)")
    row("rho_85", "rho at 85 C (trip limit): Bloch-Gruneisen vs linear", "ohm m", rho_lin(85.0), Pt.bloch_gruneisen_rho(85.0), 0.0, "differ (model)")
    row("P_bg", "Joule power with carrier resistivity", "W", P_c, P_c * rho_em / o["rho_cu_mean"], rho_se, "differ (model)")
    row("alpha", "d(rho)/dT / rho20 at 20 C", "1/K", ALPHA_CU, alpha_bg, 0.0, "differ (model)")
    row("dT_cond", "Hot-cell conduction rise (same q_v)", "K", dT_E18, dT_walk, tau_se / tau_fac, "agree")
    row("T_hot", "Hot-spot copper temperature (carrier rho + walkers)", "C", o["T_hot_C"], T_hot_em_gk, None, "differ slightly (model)")
    row("T_hot_bg", "Hot-spot copper temperature (Bloch-Gruneisen rho, walkers)", "C", o["T_hot_C"], T_hot_em, None, "differ slightly (model)")
    out["comparison"] = rows
    out["production"] = {"field_N": Nf, "moments_N": 1 << 22, "carriers_N": ncar, "walkers_N": nw,
                         "axis_rel_max_inner": float(np.max(np.abs(bax[zin, 2] / B_ax_E[zin] - 1.0))),
                         "axis_rel_max_all": float(np.max(np.abs(bax[:, 2] / B_ax_E - 1.0))),
                         "ppm_R15": o["homogeneity_ppm"], "ppm_loop_same_points": ppm_loop_same_pts,
                         "T_hot_em_fixed_point_bg": T_hot_em, "rho_hot_em": rho_hot_em, "rho_hot_bg": Pt.bloch_gruneisen_rho(T_hot_em),
                         "rho_hot_stderr_rel": rho_hot_se, "E18_at_continuum_qv": dT_E18,
                         "dT_cond_R21": o["dT_cond_inner_K"], "theta_R": Pt.THETA_R, "RRR": Pt.RRR}
    path = os.path.join(ROOT, "results", "particles.json")
    json.dump(out, open(path, "w"), indent=1, sort_keys=True)

    # ---------------- figures ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    zz = np.linspace(-0.5 * L, 0.5 * L, 201)
    axs[0].plot(zz, B_bitter_axis(R1, R2, L, C, zz), "k-", label="continuum (E4 axis)")
    axs[0].plot(ax[:, 2], bax[:, 2], "o", ms=4, label="particles, N=2^%d" % int(round(math.log(Nf, 2))))
    axs[0].set_xlabel("z [m]"); axs[0].set_ylabel("Bz on axis [T]"); axs[0].legend(); axs[0].set_title("On-axis field")
    keys = [r for r in rows if r["key"] in ("B0", "Bz_quarter", "ppm", "NI", "R", "P", "rho_mean", "rho_85", "dT_cond", "T_hot")]
    errs = [abs(r["rel_err"]) for r in keys]
    cols = ["tab:green" if r["expect"].startswith("agree") else "tab:orange" for r in keys]
    axs[1].barh(range(len(keys)), [max(e, 1e-16) for e in errs], color=cols)
    axs[1].set_yticks(range(len(keys))); axs[1].set_yticklabels([r["key"] for r in keys])
    axs[1].set_xscale("log"); axs[1].set_xlabel("|particle / continuum - 1|")
    axs[1].set_title("Relative difference (green: should agree, orange: model differs)")
    Ts = np.linspace(0, 100, 101)
    axs[2].plot(Ts, [rho_lin(t) * 1e8 for t in Ts], "k-", label="continuum linear, alpha=%.5f" % ALPHA_CU)
    axs[2].plot(Ts, [Pt.bloch_gruneisen_rho(t) * 1e8 for t in Ts], "--", label="Bloch-Gruneisen (carrier tau)")
    axs[2].errorbar([o["T_cu_mean_C"], T_hot_em], [rho_em * 1e8, rho_hot_em * 1e8],
                    yerr=[rho_em * rho_se * 1e8, rho_hot_em * rho_hot_se * 1e8], fmt="o", label="carrier ensemble")
    axs[2].set_xlabel("T [C]"); axs[2].set_ylabel("rho [1e-8 ohm m]"); axs[2].legend(fontsize=8); axs[2].set_title("Resistivity")
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "fig8_particles.png"), dpi=110)

    plot_convergence(out)
    print(json.dumps({"benchmark": bench, "comparison": [(r["key"], r["continuum"], r["particle"], r["rel_err"]) for r in rows]}, indent=1))


if __name__ == "__main__":
    if "--replot" in sys.argv:
        plot_convergence(json.load(open(os.path.join(ROOT, "results", "particles.json"))))
    else:
        main()
