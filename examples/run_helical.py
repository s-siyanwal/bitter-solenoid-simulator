"""Helical Bitter current path (plate slit / overlap) -> non-axisymmetric field errors.

Segment-engine (biot_savart.py) model of three path variants for the seed-1 optimum, over the
30 mm and 40 mm DSVs: spherical harmonics of |B| before / after the Z2/Z4 shim pairs, transverse
field, the shim ladder, and the extra shim orders needed. Radial discretisation error is estimated
from 48 vs 96 radial filaments (Richardson, order 2 observed).
Writes results/helical.json and figures/fig10_helical.png.   python examples/run_helical.py (~3 min)
"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "package"))
os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
from bittersim import helical as Hx

QUICK = "--quick" in sys.argv
NR, NR_COARSE = (12, 6) if QUICK else (96, 48)
t0 = time.time()
o = json.load(open(os.path.join(ROOT, "results", "results.json"), encoding="utf-8"))["optimal"]
g = (o["R1"], o["R2"], o["L"], o["C_A_per_m"])
nt = int(round(o["n_turns"]))
out = {"assumptions": {"n_turns": nt, "overlap_deg": Hx.OVERLAP_DEG, "bus_radius_m": o["R2"] + Hx.BUS_GAP,
                       "radial_filaments": NR, "segment_deg": 2.5, "shim_threshold_ppm": 1.0,
                       "variants": {"uniform": "uniform-pitch helix", "aligned": "flat plates, slits in one column",
                                    "rotating": "flat plates, slit advances by the overlap angle per plate"}},
       "quick": QUICK, "dsv": {}}
for dsv in (0.03, 0.04):
    key = "%dmm" % round(dsv * 1e3)
    out["dsv"][key] = {"ideal": Hx.analyse(None, *g, nt, dsv)}
    for kind in ("uniform", "aligned", "rotating"):
        a = Hx.analyse(kind, *g, nt, dsv, radial_steps=NR)
        b = Hx.analyse(kind, *g, nt, dsv, radial_steps=NR_COARSE)
        # Richardson (order 2): extrapolated value and |extrapolation - fine| as the discretisation error
        a["tesseral_richardson_err_ppm"] = dict((k, abs(v - b["tesseral_ppm"].get(k, 0.0)) / 3.0) for k, v in a["tesseral_ppm"].items())
        a["ladder_coarse_ppm"] = b["ladder_ppm"]
        # shim orders needed: the shortest cumulative ladder prefix that brings the helical coil to
        # within 1 ppm of the ideal coil's Z2/Z4-shimmed peak-to-peak (ASSUMPTION: 1 ppm budget)
        target = out["dsv"][key]["ideal"]["ladder_ppm"][1][1] + 1.0
        a["ladder_target_ppm"] = target
        reached = None
        for i, (nm, v) in enumerate(a["ladder_ppm"][1:]):
            if v <= target:
                reached = i
                break
        a["ladder_reached"] = None if reached is None else a["ladder_ppm"][reached + 1][0]
        orders = []
        for nm, keys in Hx.LADDER[1:(len(Hx.LADDER) if reached is None else reached + 1)]:
            orders += [Hx.SHIM_NAMES.get(k, k) for k in keys if k not in ("A20", "A40")]
        a["ladder_orders_needed"] = orders if reached is not None else None
        out["dsv"][key][kind] = a
        print(key, kind, [round(v, 2) for _, v in a["ladder_ppm"]], a["needed_shims"], a["ladder_orders_needed"], flush=True)
out["seconds"] = time.time() - t0
json.dump(out, open(os.path.join(ROOT, "results", "helical.json"), "w", encoding="utf-8"), indent=1, default=float)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
terms = ["A11", "B11", "A21", "B21", "A22", "B22", "A31", "B31", "A41", "B41"]
for i, kind in enumerate(("uniform", "aligned", "rotating")):
    t = out["dsv"]["30mm"][kind]["tesseral_ppm"]
    ax[0].bar(np.arange(len(terms)) + 0.27 * i, [max(abs(t.get(k, 0.0)), 1e-3) for k in terms], 0.25, label=kind)
ax[0].set_yscale("log"); ax[0].axhline(1.0, color="k", lw=0.6, ls="--")
ax[0].set_xticks(np.arange(len(terms)) + 0.27)
ax[0].set_xticklabels(["%s\n(%s)" % (k, Hx.SHIM_NAMES[k]) for k in terms], fontsize=8)
ax[0].set_ylabel("|coefficient| at r0 [ppm of B0]"); ax[0].legend(); ax[0].set_title("Tesseral terms of |B|, 30 mm DSV (dashed: 1 ppm)")
for kind, mk in (("ideal", "k-"), ("uniform", "o-"), ("aligned", "s-"), ("rotating", "^-")):
    for key, ls in (("30mm", "-"), ("40mm", ":")):
        L_ = out["dsv"][key][kind]["ladder_ppm"]
        ax[1].semilogy(range(len(L_)), [max(v, 1e-3) for _, v in L_], mk[:-1] + ls, label="%s %s" % (kind, key))
labels = [n for n, _ in out["dsv"]["30mm"]["ideal"]["ladder_ppm"]]
ax[1].set_xticks(range(len(labels))); ax[1].set_xticklabels([l.replace(" (currents solved on the ideal coil)", "") for l in labels], rotation=30, ha="right", fontsize=7)
ax[1].set_ylabel("peak-to-peak |B| over DSV [ppm]"); ax[1].legend(fontsize=7, ncol=2); ax[1].set_title("Shim ladder (ideal shims of each order, cumulative)")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "fig10_helical.png"), dpi=110)
print("helical done in %.0f s" % out["seconds"])
