"""fMRI-usefulness layer: harmonics + shims, head preset, B0 stability, RF SNR.

    python examples/run_fmri.py        (~1-2 min)

Reads the published optimum from results/results.json (it does not rerun the optimiser) and
writes results/fmri.json and figures/fig7_fmri.png. Every default input marked ASSUMPTION in
harmonics.py / stability.py / rfsnr.py is copied into the JSON under "assumptions".
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "package"))
os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
from bittersim import harmonics as H, stability as S, rfsnr as Q

t0 = time.time()
D = json.load(open(os.path.join(ROOT, "results", "results.json"), encoding="utf-8"))
o, x = D["optimal"], D["extra"]
geo = (o["R1"], o["R2"], o["L"], o["C_A_per_m"])

shim30 = H.shimmed_homogeneity(*geo, dsv=0.03)
shim40 = H.shimmed_homogeneity(*geo, dsv=0.04)
coil = H.coil_loopset(*geo)
sh_coef, _, sh_rms = H.sh_fit(coil.field_xyz, 0.015, nmax=8)
tol = H.tolerance_harmonics(*geo, dsv=0.03)
head = H.head_preset()
head_grid = head.pop("grid")

tau, LR = x["transient_tau_analytic_s"], x["L_over_R_s"]
alpha_num = S.expansion_coefficient_numeric(o["R1"], o["R2"], o["L"], o["NI_At"])
stab = {}
for mode in ("current", "voltage"):
    r = S.simulate(o, tau, LR, {"mode": mode})
    stab[mode] = {"budget_ppm": r["budget_ppm"], "budget_Hz": r["budget_Hz"], "meets_target": r["meets_target"]}
    if mode == "current":
        ts = r
req = S.requirements(o, tau)

cmp_ = Q.compare(B0=o["B0"])
cmp_lm = {str(lm): Q.compare(B0=o["B0"], loss_multiplier=lm)["gain_vs_air"] for lm in (1.0, 10.0)}
rd, rg = Q.best_detune(o["B0"])[2:]

out = {
    "assumptions": {"shim_J_A_per_m2": H.J_SHIM, "head": H.HEAD, "stability_spec": S.DEFAULT_SPEC,
                    "alpha_L_cu": S.ALPHA_L_CU, "rf": Q.DEFAULTS, "swissroll_loss_multiplier": 50.0},
    "shim_30mm": shim30, "shim_40mm": shim40,
    "sh_fit_30mm": {"coef_ppm": {k: v for k, v in sh_coef.items() if abs(v) > 1e-3}, "rms_resid_ppm": sh_rms},
    "tolerance_30mm": tol, "head_preset": head, "head_grid": head_grid,
    "stability": stab, "stability_requirements": req, "alpha_B_per_K_numeric": alpha_num,
    "rf": {"detune": cmp_["detune"], "mu": cmp_["mu"], "gain_vs_air": cmp_["gain_vs_air"],
           "gain_vs_contact": cmp_["gain_vs_contact"], "slab": cmp_["slab"], "air": cmp_["air"],
           "contact": cmp_["contact"], "gain_vs_air_by_loss_multiplier": cmp_lm},
    "seconds": time.time() - t0,
}
json.dump(out, open(os.path.join(ROOT, "results", "fmri.json"), "w", encoding="utf-8"), indent=1, default=float)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
zz = np.linspace(-0.02, 0.02, 201)
b0 = coil.field_xyz(np.c_[0 * zz, 0 * zz, zz])[:, 2]
a_s, zs, I = shim30["shim_radius_m"], shim30["shim_z_m"], shim30["shim_NI_A"]
tot = (coil + H.pair_set(a_s, zs[0], I[0]) + H.pair_set(a_s, zs[1], I[1])).field_xyz(np.c_[0 * zz, 0 * zz, zz])[:, 2]
ax[0].plot(zz * 1e3, (b0 / b0[100] - 1) * 1e6, label="unshimmed"); ax[0].plot(zz * 1e3, (tot / tot[100] - 1) * 1e6, label="Z2/Z4 shim pairs")
ax[0].set_xlabel("z [mm]"); ax[0].set_ylabel("Bz deviation on axis [ppm]"); ax[0].legend(); ax[0].set_title("Zonal shimming (30 mm DSV)")
ax[1].plot(ts["t"], ts["slow_ppm"]); ax[1].set_xlabel("t [s]"); ax[1].set_ylabel("dB/B [ppm] (slow terms)")
ax[1].set_title("B0 drift, current-regulated PSU (ripple +-%.1f ppm not drawn)" % ts["ripple_amp_ppm"], fontsize=9)
ax[2].plot(rd, rg); ax[2].axhline(1, color="k", lw=0.5); ax[2].set_xlabel("roll f0 / f_L"); ax[2].set_ylabel("SNR vs same coil over air")
ax[2].set_title("Quasi-static RF model (loss x50)")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "fig7_fmri.png"), dpi=110)
print("fmri done in %.1f s" % out["seconds"])
