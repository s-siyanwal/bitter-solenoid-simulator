"""Why the soft-Re optimum differs from the earlier hard-Re optimum, with evidence.

    python examples/run_re_investigation.py        (about 7 min on 4 CPUs: three DE runs)

Reads results/results.json (current optimum, soft Re penalty) and
results/archive/results_hard_re.json (earlier published optimum, hard Re >= 1e4) and writes
results/re_investigation.json and figures/fig12_re_investigation.png:

1. hard_re_rerun: the earlier formulation (verbatim logic of package/bittersim/optimize.py at
   commit b061e03: Re >= 1e4 as a hard constraint) run with the CURRENT evaluate_design and the
   same DE settings. If it returns the archived design vector, the move of the optimum is caused
   by the Re handling alone.
2. marginals: d(P_pump)/dv and d(P_elec)/dv at the current optimum (water velocity v only; it is
   the one design variable that does not change the discrete hole layout).
3. v_scan: P_elec + P_pump + soft Re penalty along v at the current optimum, all other variables
   fixed, and the lowest point of that line.
4. seeds: the current optimiser (examples/run_all.py settings) with seeds 2 and 3.
"""
from __future__ import print_function

import json
import os
import platform
import sys
import time
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "package"))
os.environ.setdefault("MPLBACKEND", "Agg")
warnings.filterwarnings("ignore")

import numpy as np
from scipy.optimize import differential_evolution, minimize
from bittersim import evaluate_design
from bittersim.optimize import BOUNDS, LIMITS, x_to_design, optimise, correlation_penalty, objective


def _eval(x):
    return evaluate_design(x_to_design(x))


# ---- earlier formulation: verbatim logic of optimize.py at commit b061e03 (hard Re >= 1e4) ----
HARD_LIMITS = dict(T_limit=85.0, V_max=8.0, ppm_max=100.0, dp_max=5e5, Re_min=1e4, R2_max=0.30)


def hard_constraints(res, lim=HARD_LIMITS):
    return np.array([
        (lim["T_limit"] - res["T_hot_C"]) / lim["T_limit"],
        (lim["V_max"] - res["V_total_V"]) / lim["V_max"],
        (lim["ppm_max"] - res["homogeneity_ppm"]) / lim["ppm_max"],
        (lim["dp_max"] - res["dp_Pa"]) / lim["dp_max"],
        (res["Re"] - lim["Re_min"]) / lim["Re_min"],
        (lim["R2_max"] - res["R2"]) / lim["R2_max"],
    ])


def hard_penalised(x, fixed=None):
    res = _eval(x)
    viol = np.minimum(hard_constraints(res), 0.0)
    return res["P_total_W"] / 1e3 + 1e4 * float(np.sum(viol ** 2)) + 1e2 * float(np.sum(-viol))


def hard_optimise(seed, maxiter, popsize):
    de = differential_evolution(hard_penalised, BOUNDS, args=({},), seed=seed, maxiter=maxiter,
                                popsize=popsize, tol=1e-8, polish=False)
    cons = [{"type": "ineq", "fun": (lambda x, i=i: hard_constraints(_eval(x))[i])} for i in range(6)]
    sl = minimize(lambda x: _eval(x)["P_total_W"] / 1e3, de.x, method="SLSQP", bounds=BOUNDS,
                  constraints=cons, options={"maxiter": 200, "ftol": 1e-10})
    if sl.success and np.all(hard_constraints(_eval(sl.x)) >= -1e-9) and \
            _eval(sl.x)["P_total_W"] < _eval(de.x)["P_total_W"]:
        return sl.x, de
    return de.x, de


def main():
    t0 = time.time()
    R = json.load(open(os.path.join(ROOT, "results", "results.json")))
    A = json.load(open(os.path.join(ROOT, "results", "archive", "results_hard_re.json")))
    opt = R["optimiser"]
    seed, maxiter, popsize = opt.get("seed", 1), opt.get("maxiter", 100), opt.get("popsize", 15)
    out = {"env": {"python": platform.python_version(), "numpy": np.__version__}, "cpus": os.cpu_count()}

    # 1. earlier formulation on the current code
    t = time.time()
    x_hard, de_hard = hard_optimise(seed, maxiter, popsize)
    x_arch = np.asarray(A["optimiser"]["x"])
    r_hard = _eval(x_hard)
    out["hard_re_rerun"] = {
        "x": [float(v) for v in x_hard], "x_equal_archive_bitwise": bool(np.array_equal(x_hard, x_arch)),
        "x_max_rel_diff_vs_archive": float(np.max(np.abs(x_hard / x_arch - 1.0))),
        "Re": float(r_hard["Re"]), "P_total_W": float(r_hard["P_total_W"]),
        "de_success": bool(de_hard.success), "de_nit": int(de_hard.nit), "seconds": time.time() - t}
    print("hard-Re rerun: bitwise equal to archive = %s" % out["hard_re_rerun"]["x_equal_archive_bitwise"], flush=True)

    # 2. marginals at the current optimum (water velocity only)
    x0 = np.asarray(opt["x"], dtype=float)
    def at_v(v):
        x = x0.copy()
        x[5] = v
        r = evaluate_design(x_to_design(x), homogeneity=False)
        return r, objective(r, LIMITS)
    h = 1e-4
    (rm, _), (rp, _) = at_v(x0[5] - h), at_v(x0[5] + h)
    out["marginals_W_per_m_s"] = {k: (rp[key] - rm[key]) / (2 * h) for k, key in
                                  (("dP_pump_dv", "P_pump_W"), ("dP_elec_dv", "P_elec_W"), ("dP_total_dv", "P_total_W"))}

    # 3. line scan in v (bounded below by the optimiser's lower bound)
    vs = np.linspace(BOUNDS[5][0], 1.25 * x0[5], 801)
    rows = []
    for v in vs:
        r, obj = at_v(v)
        rows.append({"v": float(v), "Re": float(r["Re"]), "P_elec_W": float(r["P_elec_W"]),
                     "P_pump_W": float(r["P_pump_W"]), "P_total_W": float(r["P_total_W"]),
                     "penalty_W": 1e3 * correlation_penalty(r, LIMITS), "objective_W": 1e3 * obj})
    r0, obj0 = at_v(x0[5])
    i = int(np.argmin([q["objective_W"] for q in rows]))
    out["v_scan"] = {"rows": rows, "optimum": {"v": float(x0[5]), "Re": float(r0["Re"]), "objective_W": 1e3 * obj0,
                                                "P_elec_W": float(r0["P_elec_W"]), "P_pump_W": float(r0["P_pump_W"])},
                     "line_min": rows[i], "gain_W": 1e3 * obj0 - rows[i]["objective_W"],
                     "gain_rel": (1e3 * obj0 - rows[i]["objective_W"]) / (1e3 * obj0), "Re_min": LIMITS["Re_min"]}
    print("v scan: line minimum at Re = %.0f, %.2f W below the optimum" % (rows[i]["Re"], out["v_scan"]["gain_W"]), flush=True)

    # 4. other seeds with the current optimiser
    out["seeds"] = [{"seed": seed, "from": "results/results.json", "Re": float(R["optimal"]["Re"]),
                     "P_total_W": float(R["optimal"]["P_total_W"]), "method": opt["method"]}]
    for s in (2, 3):
        t = time.time()
        o = optimise(seed=s, maxiter=maxiter, popsize=popsize)
        out["seeds"].append({"seed": s, "from": "this script", "Re": float(o["result"]["Re"]),
                             "P_total_W": float(o["result"]["P_total_W"]), "method": o["method"],
                             "de_success": o["de_success"], "x": [float(v) for v in o["x"]], "seconds": time.time() - t})
        print("seed %d: Re = %.0f, P_total = %.2f W (%s)" % (s, o["result"]["Re"], o["result"]["P_total_W"], o["method"]), flush=True)
    p = [q["P_total_W"] for q in out["seeds"]]
    out["seed_spread_W"] = max(p) - min(p)
    out["seed_spread_rel"] = (max(p) - min(p)) / min(p)
    out["seconds"] = time.time() - t0
    json.dump(out, open(os.path.join(ROOT, "results", "re_investigation.json"), "w"), indent=1)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    Re = np.array([q["Re"] for q in rows])
    for key, lab in (("P_pump_W", "pump power"), ("P_elec_W", "Joule (electrical) power"), ("objective_W", "P_elec + P_pump + Re penalty")):
        y = np.array([q[key] for q in rows])
        ax[0].plot(Re, y - y[np.argmin(np.abs(vs - x0[5]))], label=lab)
    ax[0].axvline(LIMITS["Re_min"], color="k", ls="--", lw=0.8, label="soft Re floor")
    ax[0].axvline(r0["Re"], color="tab:red", ls=":", lw=1.0, label="optimum (seed %d)" % seed)
    ax[0].set_xlabel("Reynolds number (water velocity varied, other variables fixed)")
    ax[0].set_ylabel("change from the optimum [W]")
    ax[0].set_ylim(-30, 100)
    ax[0].legend(fontsize=8)
    ax[0].set_title("Along the water velocity at the current optimum")
    ax[1].plot([A["optimal"]["Re"]], [A["optimal"]["P_total_W"]], "s", ms=9, label="earlier optimum (hard Re >= 1e4)")
    ax[1].plot([r_hard["Re"]], [r_hard["P_total_W"]], "x", ms=12, color="k", label="earlier formulation, current code")
    for q in out["seeds"]:
        ax[1].plot([q["Re"]], [q["P_total_W"]], "o", ms=8, label="soft Re, seed %d" % q["seed"])
    ax[1].axvline(LIMITS["Re_min"], color="k", ls="--", lw=0.8)
    ax[1].axvline(HARD_LIMITS["Re_min"], color="gray", ls="--", lw=0.8)
    ax[1].set_xlabel("Reynolds number at the optimum")
    ax[1].set_ylabel("P_elec + P_pump [W]")
    ax[1].legend(fontsize=8)
    ax[1].set_title("Optima found (dashed: soft Re floor %.0f, earlier hard floor %.0f)" % (LIMITS["Re_min"], HARD_LIMITS["Re_min"]))
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "fig12_re_investigation.png"), dpi=110)
    print("re investigation done in %.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
