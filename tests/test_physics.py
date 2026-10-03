import math
import numpy as np
import pytest
from bittersim import BitterDesign, evaluate_design, inductance as ind, swissroll, thermal, mechanics
from bittersim.constants import MU_0
from bittersim.materials import water_props, rho_cu


def test_inductance_nagaoka():
    for R, l, N in ((0.05, 0.05, 50), (0.05, 0.5, 200), (0.05, 5.0, 1000)):
        assert ind.sheet_inductance(R, l, N) == pytest.approx(ind.nagaoka_inductance(R, l, N), rel=1e-8)


def test_inductance_long_limit_and_energy():
    R, l, N, I = 0.05, 50.0, 10000, 1.0
    L = ind.nagaoka_inductance(R, l, N)
    Llong = MU_0 * N ** 2 * math.pi * R ** 2 / l
    assert L / Llong == pytest.approx(1 - 8 * R / (3 * math.pi * l), rel=1e-5)
    B = MU_0 * N / l * I
    assert 0.5 * L * I ** 2 == pytest.approx(B ** 2 / (2 * MU_0) * math.pi * R ** 2 * l, rel=2e-3)


def test_thick_inductance_thin_limit():
    t = 5e-5
    Lt = ind.coil_inductance(0.05, 0.05 + t, 0.5, 200, "uniform", nr=6)
    assert Lt == pytest.approx(ind.nagaoka_inductance(0.05 + t / 2, 0.5, 200), rel=1e-3)


def test_design_hits_target_field():
    r = evaluate_design()
    assert r["B0_numeric_T"] == pytest.approx(0.5, rel=1e-8)


def test_bitter_power_formula_matches_pdf():
    # E31 with lambda = 1: P = V0^2 L ln(R2/R1)/(2 pi rho), C = V0/(2 pi rho)
    R1, R2, L, V0, rho = 0.05, 0.15, 0.8, 0.05, 1.68e-8
    C = V0 / (2 * math.pi * rho)
    P_numeric = 0.0
    r = np.linspace(R1, R2, 200001)
    q = rho * (C / r) ** 2
    P_numeric = np.trapz(q * 2 * math.pi * r, r) * L if hasattr(np, "trapz") else np.trapezoid(q * 2 * math.pi * r, r) * L
    assert P_numeric == pytest.approx(V0 ** 2 * L * math.log(R2 / R1) / (2 * math.pi * rho), rel=1e-8)
    # E30: plate current
    d = 2e-3
    I_plate = np.sum((C / r[:-1] + C / r[1:]) / 2 * np.diff(r)) * d
    assert I_plate == pytest.approx(V0 * d * math.log(R2 / R1) / (2 * math.pi * rho), rel=1e-8)


@pytest.mark.parametrize("v,dp", [(1.2, 1.5e-3), (2.0, 2e-3), (4.0, 3e-3)])
def test_energy_balance(v, dp):
    r = evaluate_design(BitterDesign(v_flow=v, d_plate=dp), homogeneity=False)
    # V I = m_dot cp (T_out - T_in)
    lhs = r["V_total_V"] * r["I_A"]
    rhs = r["m_dot_kg_s"] * 4184.0 * (r["T_out_mixed_C"] - r["T_in"])
    assert lhs == pytest.approx(rhs, rel=1e-12)
    assert lhs == pytest.approx(r["P_elec_W"], rel=1e-12)
    assert abs(r["energy_residual"]) < 1e-12


def test_water_props_match_pdf_constants():
    w = water_props(20.0)
    assert w["rho"] == pytest.approx(998.2, rel=1e-3)
    assert w["nu"] == pytest.approx(1.004e-6, rel=1e-2)
    assert w["Pr"] == pytest.approx(7.0, rel=0.02)
    assert w["k"] == pytest.approx(0.598, rel=0.01)


def test_correlations():
    assert thermal.nusselt_dittus_boelter(1e4, 7.0) == pytest.approx(0.023 * 1e4 ** 0.8 * 7 ** 0.4)
    # Gnielinski and Dittus-Boelter agree within ~20 % for Re 1e4-1e5, Pr ~ 7
    for Re in (1e4, 3e4, 1e5):
        a = thermal.nusselt_dittus_boelter(Re, 7.0); b = thermal.nusselt_gnielinski(Re, 7.0)
        assert abs(a - b) / b < 0.2


def test_resistivity_temperature():
    assert rho_cu(20.0) == pytest.approx(1.68e-8)
    assert rho_cu(85.0) / rho_cu(20.0) == pytest.approx(1 + 3.93e-3 * 65)


def test_hotspot_above_mean():
    r = evaluate_design(homogeneity=False)
    assert r["T_hot_C"] > r["T_cu_mean_C"] > r["T_in"]
    assert r["T_hot_C"] > r["T_out_mixed_C"]


def test_swissroll_static_and_resonance():
    fL = swissroll.larmor_hz(0.5)
    assert fL / 1e6 == pytest.approx(21.2887, rel=1e-4)
    sr = swissroll.SwissRoll(fL)
    assert complex(sr.mu(0.0)) == 1.0 + 0j            # no effect on static B0
    mu = complex(sr.mu(fL))
    assert mu.real == pytest.approx(1.0, abs=1e-9)   # at exact resonance mu' = 1, mu'' = F Q
    assert mu.imag == pytest.approx(sr.F * sr.Q, rel=1e-9)
    # high-frequency limit 1 - F
    assert complex(sr.mu(1e12)).real == pytest.approx(1 - sr.F, rel=1e-3)
    # Pendry geometry mapping reproduces the requested resonance
    from bittersim.constants import C_LIGHT
    w0 = math.sqrt(sr.gap * C_LIGHT ** 2 / (2 * math.pi ** 2 * sr.eps_r * sr.r ** 3 * (sr.N - 1)))
    assert w0 == pytest.approx(2 * math.pi * fL, rel=1e-12)


def test_snr_gain_unity_without_metamaterial():
    assert float(swissroll.snr_gain(1.0 + 0j)) == pytest.approx(1.0)


def test_axial_force_compressive():
    r = evaluate_design(homogeneity=False)
    assert mechanics.axial_force(r["R1"], r["R2"], r["L"], r["C_A_per_m"]) < 0


def test_transient_steady_state():
    t, T = thermal.transient_lumped(1000.0, 0.0, 1e4, 0.01, 25.0, t_end=2000.0)
    assert T[-1] == pytest.approx(25.0 + 1000.0 * 0.01, rel=1e-6)


def test_lumped_time_constant_against_true_steady_state():
    """R27/E19: tau = C_th/(1/R_th - alpha P20); the 63 % point must use the true steady
    state (the published 29.8 s used T(60 s) as the 'steady state'; analytic is ~54.3 s)."""
    import math
    from bittersim.constants import ALPHA_CU
    P20, C_th, R_th = 11662.47, 991509.5, 5.46278e-5
    res = {"P20_W": P20, "C_th_J_K": C_th, "R_th_K_W": R_th, "T_in": 20.0, "dT_water_mixed_K": 0.5224,
           "T_hot_C": 25.95, "J_cu_inner_A_per_mm2": 4.856}
    tau = thermal.lumped_time_constant(P20, ALPHA_CU, C_th, R_th)
    assert tau == pytest.approx(C_th / (1.0 / R_th - ALPHA_CU * P20), rel=1e-12)
    s = thermal.transient_summary(res, ALPHA_CU)
    assert s["transient_tau63_s"] == pytest.approx(tau, rel=2e-3)
    assert s["transient_tau63_s"] == pytest.approx(54.3, rel=3e-3)
    assert thermal.lumped_time_constant(2.0 / ALPHA_CU, ALPHA_CU, 1.0, 1.0) == float("inf")


def test_hotspot_adiabatic_bound_is_conservative():
    """J ~ 1/r: local heating at R1 is far faster than the stack-average lumped estimate."""
    from bittersim.constants import ALPHA_CU, RHO_CU_20, DENS_CU, CP_CU
    J = 4.856e6
    t = thermal.adiabatic_time_to(85.0, 25.95, RHO_CU_20, ALPHA_CU, J, DENS_CU, CP_CU)
    # closed form vs direct integration of dens cp dT/dt = rho(T) J^2
    from scipy.integrate import quad
    t_num = quad(lambda T: DENS_CU * CP_CU / (RHO_CU_20 * (1 + ALPHA_CU * (T - 20.0)) * J ** 2), 25.95, 85.0)[0]
    assert t == pytest.approx(t_num, rel=1e-9)
    assert t < 4922.0 / 5.0
    assert thermal.adiabatic_time_to(20.0, 25.0, RHO_CU_20, ALPHA_CU, J, DENS_CU, CP_CU) == 0.0
