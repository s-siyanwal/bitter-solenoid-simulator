"""Can one joint resistance explain the DC-resistance misses across campaigns?

Diagnostic only: for each built magnet, R_measured - R_layers is divided by the
number of current-carrying interfaces. A physical joint term would give similar
per-joint values (for similar joint construction). Nothing is fitted or written
back into the model. Output: results/joint_diagnostic.json.
"""
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "package"))
sys.path.insert(0, os.path.join(ROOT, "tests"))
from bittersim import stack                       # noqa: E402
from bittersim.constants import RHO_CU_20        # noqa: E402
import test_arxiv_checks as tac                  # noqa: E402  (shares the tested layer tables)

rows = []

# EPFL bulk-machined spiral (arXiv:1901.08791): monolithic, no plate joints.
rk = 32e-3 + (np.arange(31) + 0.5) * 40e-3 / 31
epfl = [stack.Layer(r - 0.46e-3, r + 0.46e-3, 0.0, 22e-3, profile="uniform") for r in rk]
rows.append(("EPFL spiral", epfl, 0, 10.4e-3, 1.0e-3, "monolithic; gap = leads/terminals/path"))

# Claw-ZS (arXiv:2607.02813): 106 half-layers, 105 Cu conducting spacers (2 faces each).
claw = tac._claw_stack()
t, zc, _ = tac._claw_layers()
gaps = np.genfromtxt(os.path.join(ROOT, "tests", "data", "claw_zs_layers.csv"), delimiter=",",
                     comments="#", skip_header=5, filling_values=0.0)[:-1, 2] * 25.4e-3
r_spacers = float(np.sum(RHO_CU_20 * gaps / (tac.R_OUT - tac.R_IN) ** 2))   # Radia footprint
rows.append(("Claw-ZS", claw, 2 * 105, 5.3e-3 - r_spacers, 0.2e-3, "spacer bulk %.1f uOhm removed" % (r_spacers * 1e6)))

# JQI curvature and anti-bias pairs (arXiv:2008.11181): one interface per layer transition.
for kind, Rm, s in (("curv", 9.2e-3, 0.6e-3), ("bias", 13.0e-3, 0.9e-3)):
    lay = stack.mirrored(tac._jqi_layers(kind))
    rows.append(("JQI " + kind, lay, len(lay) - 2, Rm, s, "brass rho 6.6e-8 assumed"))

out = []
print("%-14s %9s %9s %9s %6s %14s" % ("campaign", "R_meas", "R_layers", "gap", "joints", "per joint"))
for name, layers, nj, Rm, sR, note in rows:
    d = stack.implied_extra_resistance(Rm, layers, nj)
    Rl = stack.stack_resistance(layers)["R_layers"]
    pj = d["per_joint"] if nj else float("nan")
    print("%-14s %7.2f mO %7.2f mO %7.2f mO %6d %8.1f +- %5.1f uO" % (
        name, Rm * 1e3, Rl * 1e3, d["R_gap"] * 1e3, nj, pj * 1e6, (sR / nj if nj else float("nan")) * 1e6))
    out.append({"campaign": name, "R_measured_ohm": Rm, "sigma_ohm": sR, "R_layers_ohm": Rl,
                "R_gap_ohm": d["R_gap"], "n_joints": nj, "per_joint_ohm": pj, "note": note})

with open(os.path.join(ROOT, "results", "joint_diagnostic.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, allow_nan=True)
