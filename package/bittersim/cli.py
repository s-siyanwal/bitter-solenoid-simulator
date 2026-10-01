"""Command-line interface:  python -m bittersim <command> [options]

  evaluate   evaluate a design (defaults = PDF initial guess)
  optimise   run the constrained optimiser
  field      on-axis / off-axis field at given points
  swissroll  Swiss-roll mu_eff and heuristic SNR gain at the Larmor frequency
"""
import argparse
import json
import sys
import numpy as np
from .design import DEFAULTS, BitterDesign, evaluate_design


def _design_args(p):
    for k, v in DEFAULTS.items():
        if isinstance(v, float):
            p.add_argument("--" + k.replace("_", "-"), type=float, default=v, dest=k)


def _fmt(res):
    keys = ["R1", "R2", "L", "d_plate", "D_hole", "v_flow", "pitch_factor", "I_A", "V_total_V", "NI_At",
            "n_turns", "P_elec_W", "P_pump_W", "homogeneity_ppm", "T_hot_C", "T_out_mixed_C", "flow_L_min",
            "dp_Pa", "Re", "h_W_m2K", "n_holes", "sigma_hoop_max_MPa", "mass_cu_kg"]
    return "\n".join("%-22s %.6g" % (k, res[k]) for k in keys if k in res)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="bittersim", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    pe = sub.add_parser("evaluate"); _design_args(pe); pe.add_argument("--json", action="store_true")
    po = sub.add_parser("optimise"); po.add_argument("--maxiter", type=int, default=60)
    po.add_argument("--seed", type=int, default=1); po.add_argument("--json", action="store_true")
    pf = sub.add_parser("field"); _design_args(pf)
    pf.add_argument("--rho", type=float, nargs="+", default=[0.0])
    pf.add_argument("--z", type=float, nargs="+", default=[0.0])
    ps = sub.add_parser("swissroll"); ps.add_argument("--B0", type=float, default=0.5)
    a = ap.parse_args(argv)
    if a.cmd is None:
        ap.print_help(); return 1
    if a.cmd in ("evaluate", "field"):
        kw = {k: getattr(a, k) for k, v in DEFAULTS.items() if isinstance(v, float)}
        res = evaluate_design(BitterDesign(**kw))
        if a.cmd == "evaluate":
            print(json.dumps(res, indent=1, default=float) if a.json else _fmt(res))
        else:
            from .fields import CoilLoops
            coil = CoilLoops(res["R1"], res["R2"], res["L"], res["C_A_per_m"], "bitter")
            rho, z = np.broadcast_arrays(np.array(a.rho), np.array(a.z))
            Br, Bz = coil.field(rho, z)
            for i in range(rho.size):
                print("rho=%.4f z=%.4f  Br=%.9f T  Bz=%.9f T" % (rho.flat[i], z.flat[i], Br[i], Bz[i]))
    elif a.cmd == "optimise":
        from .optimize import optimise
        o = optimise(seed=a.seed, maxiter=a.maxiter)
        print(json.dumps(o["result"], indent=1, default=float) if a.json else
              "method: %s\n%s\nconstraints g>=0: %s" % (o["method"], _fmt(o["result"]), np.round(o["constraints"], 5)))
    elif a.cmd == "swissroll":
        from . import swissroll as sr
        fL = sr.larmor_hz(a.B0)
        r = sr.SwissRoll(fL)
        ratios, gains, q, g = sr.best_detuning(a.B0)
        rb = sr.SwissRoll(fL * q)
        print("Larmor f = %.6f MHz" % (fL / 1e6))
        print("tuned at f_L : mu_eff = %s  Q = %.1f  gap = %.2f um" % (complex(r.mu(fL)), r.Q, r.gap * 1e6))
        print("best tuning  : f0/fL = %.4f  mu_eff = %s  heuristic SNR gain = %.3f" % (q, complex(rb.mu(fL)), g))
        print("mu_eff at DC : %s (no effect on static B0)" % complex(r.mu(0.0)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
