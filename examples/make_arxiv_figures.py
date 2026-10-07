"""Figures for the measured-data comparisons (figures/fig12-fig18, PNG + PDF).

Measured data is drawn in ink (black markers); the current model is blue and
solid, the superseded model orange and dashed (palette validated for CVD and
contrast). Every number is recomputed from the tested layer tables.
"""
import json
import math
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "package"), os.path.join(ROOT, "tests")]
warnings.filterwarnings("ignore")
from bittersim import fields, stack     # noqa: E402
import test_arxiv_checks as tac          # noqa: E402

FIG = os.path.join(ROOT, "figures")
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
NEW, OLD = "#2a78d6", "#eb6834"          # current model / superseded model
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "font.size": 9, "axes.titlesize": 10, "axes.titlecolor": INK,
    "legend.frameon": False, "legend.fontsize": 8, "lines.linewidth": 1.6, "lines.markersize": 5,
})


def meas(ax, x, y, label="measured", **kw):
    ax.plot(x, y, "o", color=INK, mfc=SURF, mew=1.2, label=label, zorder=5, **kw)


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, "%s.%s" % (name, ext)), dpi=200)
    plt.close(fig)
    print("wrote figures/%s.{png,pdf}" % name)


# ---- fig12: Claw-ZS L(f) and phase(f) ------------------------------------------
f, Zs, ph, Zm = tac._claw_ac()
w = 2 * np.pi * f
Ldc = stack.stack_inductance([ly._replace(phi=2 * math.pi, sign=0.5) for ly in tac._claw_stack()], 8, 1)
R0 = Zm[0].real
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
band = f >= 2
meas(a1, f[band], 10 * Zs[band] * np.sin(ph[band]) / w[band] * 1e6, label="measured (x10 reading)")
a1.plot(f, Zm.imag / w * 1e6, "-", color=NEW, label="PEEC AC model")
a1.axhline(Ldc * 1e6, color=OLD, ls="--", label="DC inductance (before)")
a1.set_xscale("log"); a1.set_xlabel("frequency [Hz]"); a1.set_ylabel("apparent L [uH]")
a1.set_title("Claw-ZS: inductance vs frequency"); a1.legend(loc="lower left")
meas(a2, f, np.degrees(ph))
a2.plot(f, np.degrees(np.angle(Zm)), "-", color=NEW, label="PEEC AC model")
a2.plot(f, np.degrees(np.arctan2(w * Ldc, R0)), "--", color=OLD, label="DC inductance (before)")
a2.set_xscale("log"); a2.set_xlabel("frequency [Hz]"); a2.set_ylabel("phase of Z [deg]")
a2.set_title("Claw-ZS: impedance phase (scale-free)"); a2.legend(loc="upper left")
save(fig, "fig12_clawzs_impedance")

# ---- fig13: Olsen |Z(f)| --------------------------------------------------------
lay = tac._olsen_layers()
fo = np.array(tac.OLSEN["Z_magnitude"]["f_Hz"]); Zo = np.array(tac.OLSEN["Z_magnitude"]["Z_ohm"])
Rdc = stack.stack_resistance(lay)["R_layers"]; Rx = tac.OLSEN["R_dc_measured_ohm"] - Rdc
Zac = np.abs(stack.stack_impedance(lay, fo, 8, 2) + Rx)
Zdc = np.abs(Rdc + Rx + 1j * 2 * np.pi * fo * stack.stack_inductance(lay, 8, 2))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
meas(a1, fo, Zo * 1e3)
a1.plot(fo, Zac * 1e3, "-", color=NEW, label="PEEC AC model")
a1.plot(fo, Zdc * 1e3, "--", color=OLD, label="DC inductance (before)")
a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlabel("frequency [Hz]"); a1.set_ylabel("|Z| [mOhm]")
a1.set_title("Olsen Bitter-ZS (held out): |Z(f)|"); a1.legend(loc="upper left")
a2.axhspan(-5, 5, color=GRID, alpha=0.6, lw=0)
a2.axhline(0, color=AXIS, lw=0.8)
a2.plot(fo, (Zac / Zo - 1) * 100, "-o", color=NEW, mfc=SURF, label="PEEC AC model")
a2.plot(fo, (Zdc / Zo - 1) * 100, "--s", color=OLD, mfc=SURF, label="DC inductance (before)")
a2.set_xscale("log"); a2.set_xlabel("frequency [Hz]"); a2.set_ylabel("model / measured - 1 [%]")
a2.set_title("error; grey band = +-5 %"); a2.legend(loc="upper left")
save(fig, "fig13_olsen_impedance")

# ---- fig14: Claw-ZS B(z) --------------------------------------------------------
z, Bm = tac._claw_measured()
zz = np.linspace(z.min(), z.max(), 800)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.0, 4.2), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.3]))
meas(a1, z * 1e3, Bm * 1e3, label="measured (Fig. 4a, 10 A)")
a1.plot(zz * 1e3, tac._claw_axis(zz) * 1e3, "-", color=NEW, label="bare solver, 106 half-layers (E4)")
a1.set_ylabel("Bz on axis [mT]"); a1.set_title("Claw-ZS (held out): axial field, nothing fitted"); a1.legend(loc="upper right")
a2.axhline(0, color=AXIS, lw=0.8)
a2.plot(z * 1e3, (Bm - tac._claw_axis(z)) * 1e3, "o", color=INK, mfc=SURF, mew=1.0)
a2.set_ylabel("meas - solver\n[mT]"); a2.set_xlabel("z [mm]")
save(fig, "fig14_clawzs_bz")

# ---- fig15: JQI ratios ----------------------------------------------------------
m = tac.JQI["measured_150A"]
c0, c2, b0, b1 = tac._jqi_coeffs()
items = [
    ("curvature B0", c0 / abs(m["curv_B0_uT_per_A"])),
    ("curvature B''", c2 / abs(m["curv_B2_uT_per_cm2_A"])),
    ("anti-bias B0", b0 / m["bias_B0_helmholtz_uT_per_A"]),
    ("anti-bias B'", b1 / abs(m["bias_Bp_antihelmholtz_uT_per_cm_A"])),
    ("curvature B''/B0 (shape)", (c2 / c0) / (m["curv_B2_uT_per_cm2_A"] / m["curv_B0_uT_per_A"])),
    ("anti-bias B'/B0 (shape)", (b1 / b0) / abs(m["bias_Bp_antihelmholtz_uT_per_cm_A"] / m["bias_B0_helmholtz_uT_per_A"])),
]
for kind, lab in (("curv", "curvature"), ("bias", "anti-bias")):
    items.append(("%s R (pair)" % lab, stack.stack_resistance(stack.mirrored(tac._jqi_layers(kind)))["R_layers"] * 1e3 / m["R_mOhm"][kind]))
    items.append(("%s L (one coil, 100 Hz)" % lab, stack.stack_impedance(tac._jqi_layers(kind), [100.0], 12, 2)[0].imag / (2 * math.pi * 100) * 1e6 / m["L_uH"][kind]))
fig, ax = plt.subplots(figsize=(6.4, 3.6))
y = np.arange(len(items))[::-1]
ax.axvspan(0.95, 1.05, color=GRID, alpha=0.6, lw=0)
ax.axvline(1.0, color=AXIS, lw=0.8)
for yi, (lab, r) in zip(y, items):
    ax.plot([1.0, r], [yi, yi], color=NEW, lw=1.2)
    ax.plot(r, yi, "o", color=NEW, ms=6)
    ax.text(r + (0.02 if r >= 1 else -0.02), yi, "%.3f" % r, va="center", ha="left" if r >= 1 else "right", color=INK2, fontsize=8)
ax.set_yticks(y); ax.set_yticklabels([lab for lab, _ in items])
ax.set_xlim(0.38, 1.38)
ax.set_xlabel("solver / measured  (grey band = +-5 %)")
ax.set_title("JQI Bitter Ioffe-Pritchard: shape closes, amplitude 9-12 % high")
save(fig, "fig15_jqi_ratios")

# ---- fig16: scorecard -----------------------------------------------------------
rows = json.load(open(os.path.join(ROOT, "results", "generalization.json"), encoding="utf-8"))
def kind(ch):
    if ch.startswith(("L", "abs(Z)")):
        return "inductance / impedance"
    return "resistance" if ch.startswith("R") else "field"
groups = [(k, [r for r in rows if kind(r["channel"]) == k]) for k in ("field", "inductance / impedance", "resistance")]
fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.0), sharex=True,
                         gridspec_kw=dict(height_ratios=[len(g) for _, g in groups]))
XMAX = 130
for ax, (k, sel) in zip(axes, groups):
    yy = np.arange(len(sel))[::-1]
    ax.axvspan(-5, 5, color=GRID, alpha=0.8, lw=0)
    ax.axvline(0, color=AXIS, lw=0.8)
    for yi, r in zip(yy, sel):
        before = "before" in r["model"]
        v = 100 * r["rel"]
        c = OLD if before else NEW
        ax.plot(min(v, XMAX), yi, "s" if before else "o", color=c, mfc=SURF if before else c, mew=1.4, ms=6)
        lab = "%+.1f %%" % v + (" (off scale)" if v > XMAX else "")
        ax.text(min(v, XMAX) + (-4 if v > XMAX else 4), yi + 0.28, lab, ha="right" if v > XMAX else "left",
                va="center", fontsize=7, color=INK2)
    ax.set_yticks(yy)
    ax.set_yticklabels(["%s: %s%s" % (r["campaign"], r["channel"], "  [before]" if "before" in r["model"] else "") for r in sel], fontsize=7.5)
    ax.set_ylim(-0.7, len(sel) - 0.3)
    ax.set_title(k, loc="left", fontsize=9)
axes[-1].set_xlim(-60, XMAX + 5)
axes[-1].set_xlabel("(measured - solver) / solver  [%]     grey band = +-5 %")
fig.suptitle("Bare solver across campaigns: filled blue = current model, open orange = superseded model", fontsize=9, color=INK)
save(fig, "fig16_scorecard")

# ---- fig17: joint diagnostic ----------------------------------------------------
jd = [r for r in json.load(open(os.path.join(ROOT, "results", "joint_diagnostic.json"), encoding="utf-8")) if r["n_joints"]]
fig, ax = plt.subplots(figsize=(5.6, 2.8))
yy = np.arange(len(jd))[::-1]
ax.axvline(0, color=AXIS, lw=0.8)
for yi, r in zip(yy, jd):
    v, e = r["per_joint_ohm"] * 1e6, r["sigma_ohm"] / r["n_joints"] * 1e6
    ax.errorbar(v, yi, xerr=e, fmt="o", color=NEW, ecolor=NEW, capsize=3, ms=6)
    ax.text(v + e + 6, yi, "%.0f +- %.0f uOhm (%d joints)" % (v, e, r["n_joints"]), va="center", fontsize=8, color=INK2)
ax.set_yticks(yy); ax.set_yticklabels([r["campaign"] for r in jd])
ax.set_xlim(-60, 420); ax.set_xlabel("measured-R gap per interface [uOhm]")
ax.set_title("One joint resistance does not explain all coils")
save(fig, "fig17_joint_diagnostic")

# ---- fig18: RRE 441 mm stripwound bound -----------------------------------------
I, N, r1 = 17275.0, 240, 0.2205
r2 = r1 + 15 * (16e-3 + 0.017e-3)
gaps = np.linspace(0, 6e-3, 61)
B = [fields.B_thick_center(r1, r2, 16 * 43e-3 + 15 * g, N * I / ((r2 - r1) * (16 * 43e-3 + 15 * g))) for g in gaps]
fig, ax = plt.subplots(figsize=(5.2, 2.8))
ax.plot(gaps * 1e3, B, "-", color=NEW, label="bare solver (E3), 240 turns, 17,275 A")
ax.axhline(5.0, color=INK, ls=":", lw=1.2, label="paper: 50 kG (rounded)")
ax.set_xlabel("unstated gap between pancakes [mm]"); ax.set_ylabel("B0 [T]")
ax.set_title("RRE 441 mm stripwound solenoid (1972): bound check only"); ax.legend(loc="upper right")
save(fig, "fig18_rre_stripwound_bound")
