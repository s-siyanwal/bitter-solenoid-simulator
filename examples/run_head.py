"""Head-bore preset (R1 = 190 mm, 200 mm DSV, <= 10 ppm after Z2/Z4 shims): differential evolution
with seeds 1-3, (a) with the 8 V supply limit, (b) without it, (c) minimum achievable supply voltage.
Writes results/head.json and figures/fig11_head.png.   python examples/run_head.py  (~20 min)"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "package"))
os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
from bittersim import harmonics as H

QUICK = "--quick" in sys.argv
kw = dict(maxiter=3, popsize=4) if QUICK else {}
t0 = time.time()
out = {"assumptions": {"head": H.HEAD, "bounds": dict(zip(H.HEAD_NAMES, H.HEAD_BOUNDS)),
                       "objective": "P_elec + P_pump + shim power", "dp_max_Pa": 5e5, "Re_min": 1e4, "T_limit_C": 85.0},
       "quick": QUICK, "runs": {}}
for label, v_max, what in (("power_8V", 8.0, "power"), ("power_noV", None, "power"), ("min_voltage", None, "voltage")):
    out["runs"][label] = [H.head_optimise(seed=s, v_max=v_max, what=what, **kw) for s in (1, 2, 3)]
    print(label, [(r["seed"], round(r["P_total_W"]), round(r["V_total_V"], 3), r["at_bound"]) for r in out["runs"][label]], flush=True)
out["seconds"] = time.time() - t0
json.dump(out, open(os.path.join(ROOT, "results", "head.json"), "w"), indent=1, default=float)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
for label, mk in (("power_8V", "o"), ("power_noV", "s")):
    rr = out["runs"][label]
    ax[0].plot([r["x"]["R2"] for r in rr], [r["x"]["L"] for r in rr], mk, ms=9, label=label)
lo, hi = zip(*H.HEAD_BOUNDS[:2])
ax[0].add_patch(plt.Rectangle((lo[0], lo[1]), hi[0] - lo[0], hi[1] - lo[1], fill=False, ls="--"))
ax[0].set_xlabel("R2 [m]"); ax[0].set_ylabel("L [m]"); ax[0].legend(); ax[0].set_title("Head optimum per seed (dashed: bounds)")
names = ["power_8V", "power_noV", "min_voltage"]
for i, label in enumerate(names):
    rr = out["runs"][label]
    ax[1].bar(np.arange(3) + 0.27 * i, [r["P_total_W"] / 1e3 for r in rr], 0.25, label=label)
ax[1].set_xticks(np.arange(3) + 0.27); ax[1].set_xticklabels(["seed 1", "seed 2", "seed 3"])
ax[1].set_ylabel("total power [kW]"); ax[1].legend(); ax[1].set_title("Total power per seed")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "fig11_head.png"), dpi=110)
print("head done in %.0f s" % out["seconds"])
