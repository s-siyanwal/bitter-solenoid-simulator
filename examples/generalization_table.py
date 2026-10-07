"""Cross-campaign scorecard of the bare solver, per channel, recomputed live.

Every row comes from the same layer tables the tests use. Convention:
rel = (measured - solver) / solver. Campaigns marked held_out are scored only.
Output: results/GENERALIZATION.md and results/generalization.json.
"""
import json
import math
import os
import sys
import warnings

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "package"), os.path.join(ROOT, "tests")]
warnings.filterwarnings("ignore")
from bittersim import fields, inductance, stack     # noqa: E402
import test_arxiv_checks as tac                      # noqa: E402

rows = []
SUMMARY = """
## Reading (2026-10-07)
- **Field closes wherever the winding geometry is known:** EPFL -0.3 %, Claw-ZS -1.1 %, and Bates -2 to -3 % (in bitter-residual-ml). For JQI the shape closes (B''/B0 +0.1 %) but the amplitude is 9 % high, which points to the as-built turn count in the design notebook rather than to physics.
- **AC inductance now closes:** the PEEC model brings Claw-ZS L(f) to within 1 % (it was -39 % with DC L) and cuts Olsen's 10 kHz miss from -45 % to -19 %. The remaining JQI and EPFL L gaps have unknown measurement conditions and turn counts.
- **DC resistance does not generalize:** the layer model is always low (EPFL +24 %, Claw-ZS +45 %, JQI curvature +88 %, Olsen x4.8), except for the JQI anti-bias coil (-0.8 %). The extra is coil- and setup-specific (results/joint_diagnostic.json), not one joint constant.
- **Implication for the residual corrector:** field (B) is the only channel whose bare residual is small and smooth across geometries. R and L stay out of any corrector. Olsen and Claw-ZS remain held out; their rows here are transfer scores of physics changes, not fits.
"""


def add(campaign, split, channel, measured, solver, model, note=""):
    rel = (measured - solver) / solver
    rows.append(dict(campaign=campaign, split=split, channel=channel, measured=float(measured),
                     solver=float(solver), rel=float(rel), model=model, flag=bool(abs(rel) > 0.05), note=note))


# EPFL bulk-machined spiral (arXiv:1901.08791)
E = tac.EPFL
j = E["N"] / ((E["R2"] - E["R1"]) * E["H"])
add("EPFL spiral", "published_bench", "B at 52.2 mm [G/A]", 1.28,
    fields.B_thick_axis(E["R1"], E["R2"], E["H"], j, 52.2e-3) * 1e4, "E3")
rk = E["R1"] + (np.arange(E["N"]) + 0.5) * (E["R2"] - E["R1"]) / E["N"]
ep = [stack.Layer(r - 0.46e-3, r + 0.46e-3, 0.0, E["H"], profile="uniform") for r in rk]
add("EPFL spiral", "published_bench", "R_dc [mOhm]", 10.4, stack.stack_resistance(ep)["R_layers"] * 1e3,
    "layers only", "no joints; leads/terminals")
add("EPFL spiral", "published_bench", "L [uH]", 116.0, stack.stack_inductance(ep, 1, 8) * 1e6,
    "DC filaments", "measurement frequency not stated")

# Claw-ZS (held out)
z, Bm = tac._claw_measured()
zz = np.linspace(0.0, 0.1, 2001)
Bs = tac._claw_axis(zz)
add("Claw-ZS", "held_out", "B0 peak 10 A [mT]", Bm.max() * 1e3, Bs.max() * 1e3, "E4 half-layers")
add("Claw-ZS", "held_out", "R_dc [mOhm]", 5.3, stack.stack_resistance(tac._claw_stack())["R_layers"] * 1e3,
    "layers only", "per-joint 7.8 uOhm if joints carried it")
f, Zs, ph, Zm = tac._claw_ac()
w = 2 * np.pi * f
for fk in (20.0, 2000.0):
    k = int(np.argmin(np.abs(f - fk)))
    add("Claw-ZS", "held_out", "L at %g Hz [uH] (x10 reading)" % fk, 10 * Zs[k] * np.sin(ph[k]) / w[k] * 1e6,
        Zm[k].imag / w[k] * 1e6, "PEEC AC")
Ldc = stack.stack_inductance([ly._replace(phi=2 * math.pi, sign=0.5) for ly in tac._claw_stack()], 8, 1)
k = int(np.argmin(np.abs(f - 2000.0)))
add("Claw-ZS", "held_out", "L at 2000 Hz [uH] (x10 reading)", 10 * Zs[k] * np.sin(ph[k]) / w[k] * 1e6,
    Ldc * 1e6, "DC (before)", "old comparison")

# JQI Bitter Ioffe-Pritchard (published_bench, design-notebook layer table)
m = tac.JQI["measured_150A"]
c0, c2, b0, b1 = tac._jqi_coeffs()
add("JQI curvature", "published_bench", "B0 [uT/A]", abs(m["curv_B0_uT_per_A"]), c0, "E4 stack")
add("JQI curvature", "published_bench", "B''/B0 [1/cm^2]", m["curv_B2_uT_per_cm2_A"] / m["curv_B0_uT_per_A"], c2 / c0, "E4 stack")
add("JQI anti-bias", "published_bench", "B0 Helmholtz [uT/A]", m["bias_B0_helmholtz_uT_per_A"], b0, "E4 stack")
for kind, lab in (("curv", "JQI curvature"), ("bias", "JQI anti-bias")):
    add(lab, "published_bench", "R_dc pair [mOhm]", m["R_mOhm"][kind],
        stack.stack_resistance(stack.mirrored(tac._jqi_layers(kind)))["R_layers"] * 1e3, "layers only")
    add(lab, "published_bench", "L one coil [uH]", m["L_uH"][kind],
        stack.stack_impedance(tac._jqi_layers(kind), [100.0], 12, 2)[0].imag / (2 * math.pi * 100) * 1e6, "PEEC 100 Hz")

# Olsen Bitter-ZS (held out)
lay = tac._olsen_layers()
fo = np.array(tac.OLSEN["Z_magnitude"]["f_Hz"]); Zo = np.array(tac.OLSEN["Z_magnitude"]["Z_ohm"])
Rdc = stack.stack_resistance(lay)["R_layers"]
Rx = tac.OLSEN["R_dc_measured_ohm"] - Rdc
k = int(np.argmin(np.abs(fo - 10000.0)))
add("Olsen Bitter-ZS", "held_out", "abs(Z) at 10 kHz [mOhm]", Zo[k] * 1e3,
    abs(stack.stack_impedance(lay, [fo[k]], 8, 2)[0] + Rx) * 1e3, "PEEC AC + DC contact")
add("Olsen Bitter-ZS", "held_out", "abs(Z) at 10 kHz [mOhm]", Zo[k] * 1e3,
    abs(Rdc + Rx + 1j * 2 * np.pi * fo[k] * stack.stack_inductance(lay, 8, 2)) * 1e3, "DC L (before)")
add("Olsen Bitter-ZS", "held_out", "R_dc [mOhm]", 26.5, Rdc * 1e3, "layers only", "faulty contact (paper)")

with open(os.path.join(ROOT, "results", "generalization.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1)
L = ["# Bare-solver scorecard across campaigns (generated by examples/generalization_table.py)", "",
     "rel = (measured - solver)/solver; |rel| > 5 % is flagged. held_out rows are scores only. Nothing in this table is fitted.", "",
     "| campaign | split | channel | measured | solver | rel | model | flag | note |", "|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    L.append("| %s | %s | %s | %.4g | %.4g | %+.1f %% | %s | %s | %s |" % (
        r["campaign"], r["split"], r["channel"], r["measured"], r["solver"], 100 * r["rel"], r["model"],
        "**>5 %**" if r["flag"] else "", r["note"]))
with open(os.path.join(ROOT, "results", "GENERALIZATION.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n" + SUMMARY)
print("\n".join(L[4:]))
