"""Matplotlib figures (API compatible with Matplotlib 2.2)."""
import math
import numpy as np
import matplotlib.pyplot as plt
from .constants import MU_0
from . import fields, thermal, swissroll, mechanics
from .materials import rho_cu, water_props
from .constants import ALPHA_CU


def fig_field(res, path=None):
    R1, R2, L, C = res["R1"], res["R2"], res["L"], res["C_A_per_m"]
    z = np.linspace(-L, L, 801)
    NI = res["NI_At"]
    j_uni = NI / (L * (R2 - R1))
    Rm = 0.5 * (R1 + R2)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].plot(z * 100, fields.B_bitter_axis(R1, R2, L, C, z), label="Bitter J~1/r (E4)", lw=2)
    ax[0].plot(z * 100, fields.B_thick_axis(R1, R2, L, j_uni, z), "--", label="uniform thick, same NI (E3)")
    ax[0].plot(z * 100, fields.B_thin_axis(Rm, L, NI / L, 1.0, z), ":", label="thin at R_mean, same NI (E2)")
    ax[0].axhline(MU_0 * NI / L, color="k", lw=0.8, label="ideal mu0 n I (E1)")
    ax[0].set_xlabel("z [cm]"); ax[0].set_ylabel("B_z on axis [T]"); ax[0].legend(fontsize=8)
    ax[0].set_title("On-axis field, optimal design")
    coil = fields.CoilLoops(R1, R2, L, C, "bitter", nr=12, nz=16, z_panels=4)
    rr = np.linspace(0, 0.9 * R1, 61); zz = np.linspace(-0.06, 0.06, 121)
    RR, ZZ = np.meshgrid(rr, zz)
    Br, Bz = coil.field(RR.ravel(), ZZ.ravel())
    B = np.sqrt(Br ** 2 + Bz ** 2).reshape(RR.shape)
    ppm = (B / B[60, 0] - 1) * 1e6
    lev = [-1000, -300, -100, -30, -10, 0, 10, 30, 100, 300, 1000]
    cs = ax[1].contour(ZZ * 100, RR * 100, ppm, levels=lev, cmap="coolwarm")
    ax[1].clabel(cs, fmt="%d", fontsize=7)
    t = np.linspace(0, np.pi, 100); rd = res["dsv"] / 2 * 100
    ax[1].plot(rd * np.cos(t), rd * np.sin(t), "k--", lw=1, label="DSV %.0f mm" % (res["dsv"] * 1000))
    ax[1].set_xlabel("z [cm]"); ax[1].set_ylabel("r [cm]"); ax[1].set_aspect("equal")
    ax[1].set_title("|B| deviation [ppm] (exact loop fields)"); ax[1].legend(fontsize=8)
    fig.tight_layout()
    if path: fig.savefig(path, dpi=130)
    return fig


def fig_convergence(rows, path=None):
    n = np.array([r[0] for r in rows], dtype=float); e = np.array([r[1] for r in rows])
    fig, ax = plt.subplots(figsize=(6, 4.4))
    ax.loglog(n, e, "o-", label="Biot-Savart segments vs E3")
    ax.loglog(n, e[0] * (n[0] / n) ** 2, "k--", label="slope -2 reference")
    ax.set_xlabel("angular steps per turn"); ax.set_ylabel("relative error of B_z(0)")
    ax.set_title("Thick uniform solenoid: mesh convergence"); ax.legend(); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    if path: fig.savefig(path, dpi=130)
    return fig


def fig_thermal(res, path=None):
    fl_cp = 4184.0
    m_hole = res["m_dot_kg_s"] / res["n_holes"]
    q_line = res["dT_water_inner_K"] * m_hole * fl_cp / res["L"]
    perim = math.pi * res["D_hole"] * res["fill_axial"]
    z, Tf, Ts = thermal.axial_profile(res["T_in"], q_line, m_hole, fl_cp, res["h_W_m2K"], perim, res["L"],
                                      res["dT_cond_inner_K"])
    Tw = res["T_in"] + 0.5 * res["dT_water_mixed_K"]
    t, T1 = thermal.transient_lumped(res["P20_W"], ALPHA_CU, res["C_th_J_K"], res["R_th_K_W"], Tw, t_end=60.0)
    t2, T2 = thermal.transient_lumped(res["P20_W"], ALPHA_CU, res["C_th_J_K"], res["R_th_K_W"], Tw,
                                      t_end=900.0, cooling=False)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    ax[0].plot(z * 100, Tf, label="water, innermost row (E16)")
    ax[0].plot(z * 100, Ts - res["dT_cond_inner_K"], label="wall (E17)")
    ax[0].plot(z * 100, Ts, label="copper hot spot (E17+E18)")
    ax[0].set_xlabel("z from inlet [cm]"); ax[0].set_ylabel("T [C]"); ax[0].legend()
    ax[0].set_title("Axial temperature, hottest channel")
    ax[1].plot(t2 / 60, T2, "r", label="pump failure (adiabatic)")
    ax[1].plot(t / 60, T1, "b", label="normal cooling")
    ax[1].axhline(85, color="k", ls="--", lw=0.8, label="T_limit 85 C")
    ax[1].set_xlabel("time [min]"); ax[1].set_ylabel("mean copper T [C]")
    ax[1].set_title("Lumped transient (E19)"); ax[1].legend()
    fig.tight_layout()
    if path: fig.savefig(path, dpi=130)
    return fig


def fig_swissroll(B0=0.5, path=None):
    fL = swissroll.larmor_hz(B0)
    sr = swissroll.SwissRoll(fL)
    f = np.linspace(0.97 * fL, 1.03 * fL, 2001)
    mu = sr.mu(f)
    ratios, gains, qbest, gbest = swissroll.best_detuning(B0)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    ax[0].plot(f / 1e6, mu.real, label="mu' (tuned to f_L)")
    ax[0].plot(f / 1e6, mu.imag, label="mu''")
    ax[0].axvline(fL / 1e6, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel("f [MHz]"); ax[0].set_ylabel("mu_eff (E25)"); ax[0].legend()
    ax[0].set_title("Swiss roll mu_eff (RF only; mu_eff(0) = 1)")
    ax[1].plot(ratios, gains)
    ax[1].axvline(qbest, color="r", ls="--", lw=0.8)
    ax[1].set_xlabel("roll resonance / Larmor frequency"); ax[1].set_ylabel("heuristic SNR gain (E28)")
    ax[1].set_title("SNR gain vs tuning (illustrative model)")
    fig.tight_layout()
    if path: fig.savefig(path, dpi=130)
    return fig


def fig_mechanics(res, path=None):
    r, Bz, sig = mechanics.hoop_stress_profile(res["R1"], res["R2"], res["L"], res["C_A_per_m"],
                                               res["fill_total"])
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ax[0].plot(r * 100, res["C_A_per_m"] / (res["fill_total"] * r) / 1e6)
    ax[0].set_xlabel("r [cm]"); ax[0].set_ylabel("J_cu [A/mm^2]"); ax[0].set_title("Bitter current density J ~ 1/r")
    ax[1].plot(r * 100, sig / 1e6, label="hoop stress (E21)")
    ax2 = ax[1].twinx(); ax2.plot(r * 100, Bz, "g--", label="B_z(r, z=0)"); ax2.set_ylabel("B_z [T]")
    ax[1].set_xlabel("r [cm]"); ax[1].set_ylabel("sigma_theta [MPa]"); ax[1].set_title("Midplane hoop stress")
    fig.tight_layout()
    if path: fig.savefig(path, dpi=130)
    return fig, (r, Bz, sig)
