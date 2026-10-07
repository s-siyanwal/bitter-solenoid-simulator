"""Cross-code cases: bittersim vs Elmer (and GetDP where verified) vs measurement.

Usage:  python crosscode/run_cases.py [A] [B] [C]   -> crosscode/results/<case>.json
Work directories go to $XC_WORK (default: crosscode/_work, git-ignored).
Every FEM model uses the same layer table as the corresponding bittersim test; nothing
is tuned to the measurement.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [os.path.join(ROOT, "package"), os.path.join(ROOT, "tests"), HERE]
import fem_axi as F                      # noqa: E402
from bittersim import fields             # noqa: E402

WORK = os.environ.get("XC_WORK", os.path.join(HERE, "_work"))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def save(name, obj):
    with open(os.path.join(OUT, name + ".json"), "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1)
    print("wrote crosscode/results/%s.json" % name)


# ---------------------------------------------------------------- case A: EPFL spiral
def case_A():
    R1, R2, H, N = 0.032, 0.072, 0.022, 31
    J = N * 1.0 / ((R2 - R1) * H)                 # 1 A drive -> field per ampere
    reg = [dict(r1=R1, r2=R2, z1=-H / 2, z2=H / 2, kind="js", J=J)]
    zs = [0.0, 0.0261, 0.0522]
    rows = []
    for lc in (0.002, 0.001, 0.0005):              # mesh refinement study
        t = time.time()
        out = F.run_elmer(os.path.join(WORK, "A_%g" % lc), reg, [(0.0, z) for z in zs],
                          Rbox=1.5, Zbox=1.5, lc_reg=lc, lc_air=0.05, lc_axis=lc)
        rows.append(dict(lc=lc, mesh=out["nodes_elements"], seconds=time.time() - t,
                         Bz=[b[3] for b in out["b"]]))
    gd = []
    for lc in (0.002, 0.001, 0.0005):
        t = time.time()
        out = F.run_getdp(os.path.join(WORK, "Ag_%g" % lc), reg, [(0.0, z) for z in zs],
                          Rbox=1.5, Zbox=1.5, lc_reg=lc, lc_air=0.05, lc_axis=lc)
        gd.append(dict(lc=lc, mesh=out["nodes_elements"], seconds=time.time() - t,
                       Bz=[abs(b[3]) for b in out["b"]]))
    exact = [float(fields.B_thick_axis(R1, R2, H, J, z)) for z in zs]
    save("A_epfl", dict(z=zs, exact_E3=exact, elmer=rows, getdp=gd, measured_G_per_A_at_52mm=1.28,
                        getdp_G_per_A_at_52mm=gd[-1]["Bz"][2] * 1e4,
                        bittersim_G_per_A_at_52mm=exact[2] * 1e4,
                        elmer_G_per_A_at_52mm=rows[-1]["Bz"][2] * 1e4))


# ---------------------------------------------------------------- case B: Claw-ZS profile
def case_B():
    import test_arxiv_checks as tac
    I = tac.I_CLAW
    regs = []
    for ly in tac._claw_stack():                   # 106 half-layers, smeared to 1/2 current
        C = 0.5 * I / (ly.t * math.log(ly.r_out / ly.r_in))
        regs.append(dict(r1=ly.r_in, r2=ly.r_out, z1=ly.z - ly.t / 2, z2=ly.z + ly.t / 2, kind="bitter", C=C))
    z, Bm = tac._claw_measured()
    probes = [(0.0, float(zz)) for zz in z]
    t = time.time()
    out = F.run_elmer(os.path.join(WORK, "B"), regs, probes, Rbox=2.0, Zbox=2.0,
                      lc_reg=0.0005, lc_air=0.08, lc_axis=0.001)
    Be = np.array([b[3] for b in out["b"]])
    Bs = tac._claw_axis(z)
    rms = lambda d: float(np.sqrt(np.mean(d ** 2)))
    save("B_clawzs", dict(z_m=z.tolist(), measured_T=Bm.tolist(), bittersim_T=Bs.tolist(),
                          elmer_T=Be.tolist(), mesh=out["nodes_elements"], seconds=time.time() - t,
                          rms_meas_minus_bittersim_mT=rms(Bm - Bs) * 1e3,
                          rms_meas_minus_elmer_mT=rms(Bm - Be) * 1e3,
                          rms_elmer_minus_bittersim_mT=rms(Be - Bs) * 1e3,
                          peak_meas_mT=float(Bm.max() * 1e3), peak_bittersim_mT=float(Bs.max() * 1e3),
                          peak_elmer_mT=float(Be.max() * 1e3)))


# ---------------------------------------------------------------- case C: NASA MDF coil 1
# disk ID/OD from drawing RES/MAG-12-01 are DIAMETERS (bore 184.2 mm dia); radii = /2
NASA_COILS = [(d1 / 2, d2 / 2) for d1, d2 in
              [(190.3e-3, 213.1e-3), (214.7e-3, 239.5e-3), (241.1e-3, 268.7e-3), (270.3e-3, 303.7e-3)]]
NASA_I_3KW = 31.86
NASA_GAP = 72e-3
SHELLS = {"far_legacy": (165.688e-3, 170.80e-3), "close_pdf_1026": (152.8e-3, 157.8e-3)}


def nasa_regions(t_cu, N, shell, mur, insul=0.0381e-3):
    pitch = t_cu + insul
    Lh = (N / 2) * pitch                           # each half-coil (split at the midplane gap)
    regs = []
    for R1, R2 in NASA_COILS:
        C = (N / 2) * NASA_I_3KW / (Lh * math.log(R2 / R1))
        for sgn in (1, -1):
            z1, z2 = sorted((sgn * NASA_GAP / 2, sgn * (NASA_GAP / 2 + Lh)))
            regs.append(dict(r1=R1, r2=R2, z1=z1, z2=z2, kind="bitter", C=C))
    if mur > 1:
        ri, ro = SHELLS[shell]
        regs.append(dict(r1=ri, r2=ro, z1=-0.109, z2=0.109, kind="iron", mur=mur))
        for sgn in (1, -1):                        # annular heads, 7 mm, r 94..172 mm, outboard
            z1, z2 = sorted((sgn * 0.109, sgn * 0.116))
            regs.append(dict(r1=0.094, r2=0.172, z1=z1, z2=z2, kind="iron", mur=mur))
    return regs, Lh


def case_C():
    rows = []
    meas = 136.9e-3
    ladder = [("air core, t=1.0, N=113", 1.0e-3, 113, "far_legacy", 1.0, 50.84e-3),
              ("steel mu_r=200, far shell", 1.0e-3, 113, "far_legacy", 200.0, 85.52e-3),
              ("steel mu_r=2000, far shell", 1.0e-3, 113, "far_legacy", 2000.0, 93.65e-3),
              ("steel mu_r=2000, as-built shell", 1.0e-3, 113, "close_pdf_1026", 2000.0, 100.1e-3),
              ("thin disks t=0.20, N=140, as-built shell", 0.20e-3, 140, "close_pdf_1026", 2000.0, 132.4e-3),
              ("thin disks t=0.20, N=150, far shell", 0.20e-3, 150, "far_legacy", 2000.0, 137.2e-3)]
    for k, (lab, t, N, sh, mur, prior) in enumerate(ladder):
        regs, Lh = nasa_regions(t, N, sh, mur)
        tic = time.time()
        out = F.run_elmer(os.path.join(WORK, "C%d" % k), regs, [(0.0, 0.0)], Rbox=1.5, Zbox=1.5,
                          lc_reg=0.002, lc_air=0.08, lc_axis=0.002)
        Bz = abs(out["b"][0][3])
        og = F.run_getdp(os.path.join(WORK, "Cg%d" % k), regs, [(0.0, 0.0)], Rbox=1.5, Zbox=1.5,
                         lc_reg=0.002, lc_air=0.08, lc_axis=0.002)
        Bg = abs(og["b"][0][3])
        rows.append(dict(case=lab, getdp_peak_mT=Bg * 1e3, t_cu_mm=t * 1e3, N=N, shell=sh, mu_r=mur, half_length_mm=Lh * 1e3,
                         elmer_peak_mT=Bz * 1e3, earlier_inhouse_fem_mT=prior * 1e3,
                         rel_meas_vs_elmer=(meas - Bz) / Bz, rel_meas_vs_inhouse=(meas - prior) / prior,
                         mesh=out["nodes_elements"], seconds=time.time() - tic))
        print("%-42s getdp %.1f  elmer %.1f mT (in-house %.1f)  meas-sol %+.1f %%" % (lab, Bg * 1e3, Bz * 1e3, prior * 1e3, 100 * rows[-1]["rel_meas_vs_elmer"]))
    # air-core cross-check against the bittersim closed form (E4 superposition)
    t, N = 1.0e-3, 113
    regs, Lh = nasa_regions(t, N, "far_legacy", 1.0)
    e4 = sum(fields.B_bitter_axis(g["r1"], g["r2"], g["z2"] - g["z1"], g["C"], -(g["z1"] + g["z2"]) / 2) for g in regs)
    save("C_nasa", dict(measured_peak_3kW_mT=meas * 1e3, ladder=rows, bittersim_air_core_E4_mT=float(e4) * 1e3))


if __name__ == "__main__":
    todo = sys.argv[1:] or ["A", "B", "C"]
    for c in todo:
        {"A": case_A, "B": case_B, "C": case_C}[c]()
