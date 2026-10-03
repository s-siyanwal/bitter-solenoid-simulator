"""Tests for the fMRI layer: harmonics/shims, stability, quasi-static RF SNR."""
import math
import numpy as np
import pytest
from scipy.special import j1
from scipy.integrate import quad
from bittersim import harmonics as H, stability as S, rfsnr as Q
from bittersim.constants import MU_0, ALPHA_CU
from bittersim.fields import B_bitter_center

GEO = (0.05, 0.3, 1.1, 230440.0)          # close to the seed-1 optimum


# ------------------------------------------------------------------ (a) harmonics
def test_axis_zonal_matches_full_sh_fit():
    zon, B0 = H.zonal_axis(*GEO, r0=0.015)
    assert B0 == pytest.approx(B_bitter_center(*GEO), rel=1e-12)
    coef, B0f, rms = H.sh_fit(H.coil_loopset(*GEO).field_xyz, 0.015, nmax=8)
    assert coef["A20"] == pytest.approx(zon[2], abs=0.02)        # ppm; loop quadrature floor ~0.01 ppm
    assert coef["A40"] == pytest.approx(zon[4], abs=0.02)
    assert abs(coef["A11"]) < 1e-6 and abs(coef["A10"]) < 1e-6   # axisymmetric, mid-plane symmetric
    assert rms < 1e-2


def test_sh_fit_recovers_known_harmonics():
    # synthetic field Bz = 1 + 2e-6 (x/r0) + 5e-6 P2 (r/r0)^2 ; B from gradient of a harmonic potential
    r0 = 0.02
    def f(p):
        x, y, z = p[:, 0], p[:, 1], p[:, 2]
        Bz = 1 + 2e-6 * x / r0 + 5e-6 * (z ** 2 - 0.5 * (x ** 2 + y ** 2)) / r0 ** 2
        return np.column_stack([0 * x, 0 * x, Bz])
    c, B0, rms = H.sh_fit(f, r0, nmax=4)
    assert c["A11"] == pytest.approx(2.0, rel=1e-9)
    assert c["A20"] == pytest.approx(5.0, rel=1e-9)


def test_shims_null_z2_z4_and_improve_ppm():
    s = H.shimmed_homogeneity(*GEO, dsv=0.03)
    assert abs(s["zonal_ppm_shimmed"][2]) < 1e-6 and abs(s["zonal_ppm_shimmed"][4]) < 1e-6
    assert s["ppm_shimmed"] < s["ppm_unshimmed"] / 50.0
    assert s["shim_radius_m"] < GEO[0]
    # shim power formula
    assert s["shim_power_W"] == pytest.approx(
        sum(2 * 1.68e-8 * H.J_SHIM * abs(I) * 2 * math.pi * s["shim_radius_m"] for I in s["shim_NI_A"]), rel=1e-12)


def test_shimmed_ppm_is_converged_in_quadrature():
    a = H.shimmed_homogeneity(*GEO, dsv=0.03)["ppm_shimmed"]
    b = H.shimmed_homogeneity(*GEO, dsv=0.03, nr=24, nz=64, z_panels=8)["ppm_shimmed"]
    assert a == pytest.approx(b, abs=0.02)


def test_tolerance_offset_gives_expected_a11():
    # a lateral shift dx of a field with zonal A20 creates A11 = A20 * dx / r0 (to first order)
    r0, dx = 0.015, 0.5e-3
    t = H.tolerance_harmonics(*GEO, dsv=2 * r0, dx=dx, tilt=0.0)
    zon, _ = H.zonal_axis(*GEO, r0=r0)
    assert t["tesseral_ppm"]["A11"] == pytest.approx(zon[2] * dx / r0, rel=2e-3)


def test_head_preset_reports_honest_numbers():
    h = H.head_preset(R2_grid=[0.7], L_grid=[2.6])
    assert h["feasible"] and h["R1"] == 0.19 and h["dsv_m"] == 0.2
    assert h["ppm_shimmed"] <= H.HEAD["ppm_target"] < h["ppm_unshimmed"]
    assert "copper_cost_usd" not in h                      # no built-in prices
    h2 = H.head_preset(R2_grid=[0.7], L_grid=[2.6], copper_usd_per_kg=1.0)
    assert h2["copper_cost_usd"] == pytest.approx(h2["mass_cu_kg"])


# ------------------------------------------------------------------ (b) stability
def test_expansion_coefficient_is_minus_alpha_L():
    a = S.expansion_coefficient_numeric(0.05, 0.3, 1.1, 4.5e5)
    assert a == pytest.approx(-S.ALPHA_L_CU, rel=1e-4)        # ~ -16.5 ppm/K


def test_first_order_lag_step_response():
    t = np.linspace(0, 500, 5001)
    y = S.first_order_lag(t, np.ones_like(t), 50.0, 0.0)
    assert np.interp(50.0, t, y) == pytest.approx(1 - math.exp(-1), abs=1e-3)


def _res():
    return {"B0": 0.5, "T_in": 20.0, "T_cu_mean_C": 20.9}


def test_stability_budget_components():
    r = S.simulate(_res(), tau=54.3, L_over_R=0.657, spec={"water_amp_K": 0.0, "water_drift_K_per_h": 0.0})
    b = r["budget_ppm"]
    assert b["psu_ripple"] == pytest.approx(2.0) and b["thermal_expansion"] == pytest.approx(0.0, abs=1e-12)
    assert b["psu_drift"] == pytest.approx(2.0 * 600 / 3600)
    assert r["budget_Hz"]["total"] == pytest.approx(b["total"] * r["f_larmor_Hz"] / 1e6)
    # pure slow water drift, long run: copper follows and B falls by alpha_L per K
    r2 = S.simulate(_res(), tau=1.0, L_over_R=0.657, spec={"ripple_ppm": 0, "drift_ppm_per_h": 0,
                    "water_amp_K": 0, "water_drift_K_per_h": 3.6, "duration_s": 1000})
    assert r2["budget_ppm"]["thermal_expansion"] == pytest.approx(16.5 * 0.999, rel=2e-3)


def test_voltage_mode_resistance_drift_dominates():
    r = S.simulate(_res(), tau=54.3, L_over_R=0.657, spec={"mode": "voltage"})
    c = S.simulate(_res(), tau=54.3, L_over_R=0.657, spec={"mode": "current"})
    assert r["budget_ppm"]["resistance_drift"] > 100 * c["budget_ppm"]["thermal_expansion"]
    q = S.requirements(_res(), tau=54.3)
    assert q["max_copper_dT_K_current_mode"] == pytest.approx(1e-6 / 16.5e-6)


# ------------------------------------------------------------------ (c) RF model
def test_slab_transfer_matches_direct_4x4_solve():
    k, mu, t, z1 = np.array([10.0, 50.0]), 5.7 + 0.4j, 0.05, 0.002
    tau, A, Bp, q = Q.slab_transfer(k, mu, t, z1)
    for i in range(2):
        kk, qq = k[i], q[i]
        E = np.exp(-qq * t); inc = np.exp(-kk * z1)
        M = np.array([[1, -1, -E, 0], [kk, mu * qq, -mu * qq * E, 0],
                      [0, E, 1, -1], [0, -mu * qq * E, mu * qq, kk]], dtype=complex)
        rhs = np.array([-inc, kk * inc, 0, 0], dtype=complex)
        r, a_, b_, tt = np.linalg.solve(M, rhs)
        assert tt == pytest.approx(tau[i], rel=1e-10)


def test_unit_mu_slab_is_air_and_b1_matches_loop_formula():
    s, a = Q.evaluate(1.0 + 0j), Q.evaluate(None)
    assert s["B1_T_per_A"] == pytest.approx(a["B1_T_per_A"], rel=1e-12) and s["R_slab"] == 0.0
    p = Q.DEFAULTS; z = p["gap"] + p["t_slab"] + p["depth"]
    assert a["B1_T_per_A"] == pytest.approx(MU_0 * p["a"] ** 2 / (2 * (p["a"] ** 2 + z ** 2) ** 1.5), rel=1e-4)


def test_tissue_loss_half_space_closed_form():
    a, s = Q.DEFAULTS["a"], Q.DEFAULTS["sigma"]
    r = Q.evaluate(None, {"gap": 0.0, "t_slab": 0.0})
    w = 2 * math.pi * r["f_Hz"]
    assert quad(lambda x: j1(x) ** 2 / x ** 2, 0, 4000, limit=4000)[0] == pytest.approx(4 / (3 * math.pi), rel=1e-6)
    assert r["R_tissue"] == pytest.approx(s * w ** 2 * MU_0 ** 2 * a ** 3 / 3.0, rel=1e-3)


def test_lossless_high_mu_slab_guides_flux_and_losses_cut_gain():
    lossless = Q.evaluate(50.0 + 0j)["B1_T_per_A"]
    assert lossless > 2 * Q.evaluate(None)["B1_T_per_A"]
    g1 = Q.compare(loss_multiplier=1.0, detune=1.03)["gain_vs_air"]
    g50 = Q.compare(loss_multiplier=50.0, detune=1.03)["gain_vs_air"]
    assert g1 > g50
