"""Command-line interface:  python -m bittersim <command> [options]

  evaluate   evaluate a design (defaults = PDF initial guess)
  optimise   run the constrained optimiser
  field      on-axis / off-axis field at given points
  swissroll  Swiss-roll mu_eff and heuristic SNR gain at the Larmor frequency
  emulate    mesoscopic emulation (not molecular dynamics) around a design
  catalog    print the material allow-list
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


def _catalog_args(p):
    p.add_argument("--catalog-conductor", default="ofhc_cu")
    p.add_argument("--catalog-coolant", default="di_water")
    p.add_argument("--catalog-insulator", default="polyimide_kapton")
    p.add_argument("--catalog-housing", default="ss304")
    p.add_argument("--temper", default="hard", choices=["annealed", "hard"])
    p.add_argument("--T-amb", dest="T_amb", type=float, default=20.0)
    p.add_argument("--p-site", dest="p_site_Pa", type=float, default=101325.0)
    p.add_argument("--aeration", default="deaerated", choices=["aerated", "deaerated"])
    p.add_argument("--conductive-bore", action="store_true")
    p.add_argument("--radial-rho", action="store_true")


def _design_from_args(a):
    kw = {}
    for k, v in DEFAULTS.items():
        if isinstance(v, float):
            kw[k] = getattr(a, k)
    return BitterDesign(**kw)


def _fmt(res):
    keys = ["R1", "R2", "L", "d_plate", "D_hole", "v_flow", "pitch_factor", "I_A", "V_total_V", "NI_At",
            "n_turns", "P_elec_W", "P_pump_W", "homogeneity_ppm", "T_hot_C", "T_out_mixed_C", "flow_L_min",
            "dp_Pa", "Re", "h_W_m2K", "n_holes", "sigma_hoop_max_MPa", "mass_cu_kg"]
    return "\n".join("%-22s %.6g" % (k, res[k]) for k in keys if k in res)


def _print_warnings(res):
    print("warnings:")
    for w in res.get("warnings", []):
        print("- %s" % w)


def _json(obj):
    print(json.dumps(obj, indent=1, default=float))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="bittersim", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    pe = sub.add_parser("evaluate")
    _design_args(pe)
    _catalog_args(pe)
    pe.add_argument("--json", action="store_true")
    po = sub.add_parser("optimise")
    po.add_argument("--maxiter", type=int, default=60)
    po.add_argument("--seed", type=int, default=1)
    po.add_argument("--json", action="store_true")
    pf = sub.add_parser("field")
    _design_args(pf)
    pf.add_argument("--rho", type=float, nargs="+", default=[0.0])
    pf.add_argument("--z", type=float, nargs="+", default=[0.0])
    ps = sub.add_parser("swissroll")
    ps.add_argument("--B0", type=float, default=0.5)
    pm = sub.add_parser("emulate")
    _design_args(pm)
    _catalog_args(pm)
    pm.add_argument("--realizations", type=int, default=200)
    pm.add_argument("--seed", type=int, required=True)
    pm.add_argument("--bandwidth-lo", type=float, default=1.0)
    pm.add_argument("--bandwidth-hi", type=float, default=1.0e4)
    pm.add_argument("--contact-median", type=float, default=1.0e-6)
    pm.add_argument("--hooge", action="store_true")
    pm.add_argument("--json", action="store_true")
    sub.add_parser("catalog")
    a = ap.parse_args(argv)
    if a.cmd is None:
        ap.print_help()
        return 1
    if a.cmd in ("evaluate", "field"):
        design = _design_from_args(a)
        if a.cmd == "field":
            res = evaluate_design(design)
            from .fields import CoilLoops
            coil = CoilLoops(res["R1"], res["R2"], res["L"], res["C_A_per_m"], "bitter")
            rho, z = np.broadcast_arrays(np.array(a.rho), np.array(a.z))
            Br, Bz = coil.field(rho, z)
            for i in range(rho.size):
                print("rho=%.4f z=%.4f  Br=%.9f T  Bz=%.9f T" % (rho.flat[i], z.flat[i], Br[i], Bz[i]))
        else:
            from .emulation import evaluate_with_catalog
            res = evaluate_with_catalog(
                design,
                conductor_id=a.catalog_conductor,
                coolant_id=a.catalog_coolant,
                insulator_id=a.catalog_insulator,
                housing_id=a.catalog_housing,
                temper=a.temper,
                T_amb=a.T_amb,
                p_site_Pa=a.p_site_Pa,
                aerated=(a.aeration == "aerated"),
                conductive_bore=a.conductive_bore,
                radial_rho=a.radial_rho,
            )
            if a.json:
                _json(res)
            else:
                print(_fmt(res))
                _print_warnings(res)
                if res.get("constraints"):
                    print("constraints:")
                    for c in res["constraints"]:
                        print("- %s %s  value %.6g  limit %.6g %s" % (
                            c["status"], c["name"], c["value"], c["limit"], c["unit"]))
    elif a.cmd == "optimise":
        from .optimize import optimise
        o = optimise(seed=a.seed, maxiter=a.maxiter)
        if a.json:
            _json(o["result"])
        else:
            print("method: %s\n%s\nconstraints g>=0: %s" % (
                o["method"], _fmt(o["result"]), np.round(o["constraints"], 5)))
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
        print("warning: Swiss-roll SNR gain is a heuristic, not an MRI SNR measurement.")
    elif a.cmd == "emulate":
        from .emulation import emulate
        out = emulate(
            _design_from_args(a),
            realizations=a.realizations,
            seed=a.seed,
            conductor_id=a.catalog_conductor,
            coolant_id=a.catalog_coolant,
            insulator_id=a.catalog_insulator,
            housing_id=a.catalog_housing,
            temper=a.temper,
            T_amb=a.T_amb,
            p_site_Pa=a.p_site_Pa,
            aerated=(a.aeration == "aerated"),
            conductive_bore=a.conductive_bore,
            radial_rho=a.radial_rho,
            f_lo=a.bandwidth_lo,
            f_hi=a.bandwidth_hi,
            contact_median_ohm=a.contact_median,
            hooge=a.hooge,
        )
        if a.json:
            _json(out)
        else:
            print(out["label"])
            print("seed %d   realizations %d" % (out["seed"], out["n_realizations"]))
            print("%-22s %12s %12s %12s %12s" % ("quantity", "nominal", "p05", "p50", "p95"))
            nom = out["nominal_summary"]
            mapping = [
                ("B0_T", "B0_T"),
                ("V_V", "V_V"),
                ("P_W", "P_W"),
                ("T_hot_C", "T_hot_C"),
                ("ppm", "ppm"),
                ("v_d_m_s", "v_d_m_s"),
                ("mean_free_path_m", "mean_free_path_m"),
                ("V_johnson_V", "V_johnson_V"),
            ]
            for key, nk in mapping:
                p = out["percentiles"][key]
                base = nom.get(nk)
                base_s = "%.6g" % base if base is not None else "-"
                print("%-22s %12s %12.6g %12.6g %12.6g" % (key, base_s, p["p05"], p["p50"], p["p95"]))
            tr = out["trip_probability"]
            print("trip T_hot>85C %.3f   ONB margin<0 %.3f   stress>0.6 yield %.3f" % (
                tr["T_hot_above_85C"], tr["onb_margin_below_0"], tr["stress_above_0.6_yield"]))
            print("active: %s" % (", ".join(out["active_constraints"]) or "none"))
    elif a.cmd == "catalog":
        from .catalog import describe
        print(describe())
    return 0


if __name__ == "__main__":
    sys.exit(main())
