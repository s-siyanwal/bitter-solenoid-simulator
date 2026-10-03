"""Mesoscopic emulation. The published continuum numbers stay put."""
import math
import pytest
from bittersim import BitterDesign, evaluate_design
from bittersim.catalog import coolant_props, get_conductor
from bittersim.emulation import (
    E_CHARGE, M_ELECTRON, K_BOLTZMANN, drude_rve, johnson_voltage,
    langevin_ensemble, emulate, saturation_temperature_C,
)
from bittersim.materials import water_props
from bittersim.cli import main


def test_drude_identities():
    rho = 1.68e-8
    n = 8.47e28
    J = 4.856e6
    d = drude_rve(J, rho, 26.0, 0.5, n, 0.006, 1.57e6, kohler=False)
    assert d["rho_from_tau"] == pytest.approx(rho, rel=1e-12)
    assert rho == pytest.approx(M_ELECTRON / (n * E_CHARGE ** 2 * d["tau_s"]), rel=1e-12)
    assert d["v_d_m_s"] == pytest.approx(J / (n * E_CHARGE), rel=1e-12)
    assert abs(d["power_density_rel_error"]) < 0.01
    assert d["J_dot_E_W_m3"] == pytest.approx(d["q_continuum_W_m3"], rel=1e-12)


def test_hall_angle_small_at_half_tesla():
    c = get_conductor("ofhc_cu")
    d = drude_rve(5e6, c["rho20"], 26.0, 0.5, c["n_m3"], 0.006, c["v_fermi"], kohler=True)
    assert d["omega_c_tau"] < 1e-2
    assert d["model_grade_kohler"] == "approximate"


def test_johnson_integral():
    T_C, R, flo, fhi = 20.0, 0.0212356, 1.0, 1.0e4
    j = johnson_voltage(T_C, R, flo, fhi)
    Sv = 4.0 * K_BOLTZMANN * (T_C + 273.15) * R
    assert j["S_v_V2_per_Hz"] == pytest.approx(Sv, rel=1e-12)
    assert j["V_rms_V"] == pytest.approx(math.sqrt(Sv * (fhi - flo)), rel=1e-12)


def test_catalog_rejects_unknown():
    with pytest.raises(ValueError) as exc:
        get_conductor("unobtainium")
    assert "ofhc_cu" in str(exc.value)
    with pytest.raises(ValueError) as exc2:
        coolant_props("beer", 20.0)
    assert "di_water" in str(exc2.value)


def test_glycol_lowers_reynolds():
    w = water_props(20.0)
    g = coolant_props("water_glycol_30", 20.0)
    v, D = 1.7, 0.006
    Re_w = w["rho"] * v * D / w["mu"]
    Re_g = g["rho"] * v * D / g["mu"]
    assert Re_g < Re_w
    assert g["k"] < w["k"]


def test_default_evaluate_matches_results_md():
    r = evaluate_design()
    # Printed precision in results/RESULTS.md for the PDF initial guess.
    assert r["P_elec_W"] / 1000.0 == pytest.approx(17.974, abs=5e-4)
    assert r["V_total_V"] == pytest.approx(19.537, abs=5e-4)
    assert r["T_hot_C"] == pytest.approx(21.153, abs=5e-4)
    assert "warnings" in r
    assert any("Florida-Bitter" in w for w in r["warnings"])
    assert any("not stress-limited" in w for w in r["warnings"])


def test_default_catalog_does_not_move_numbers():
    plain = evaluate_design()
    from bittersim.emulation import evaluate_with_catalog
    cat = evaluate_with_catalog(BitterDesign())
    for key in ("P_elec_W", "V_total_V", "T_hot_C", "Re", "homogeneity_ppm", "mass_cu_kg", "C_th_J_K"):
        assert cat[key] == pytest.approx(plain[key], rel=0, abs=0)


def test_alloy_mass_and_glycol_follow_catalog():
    from bittersim.emulation import evaluate_with_catalog
    d = BitterDesign()
    cu = evaluate_with_catalog(d, homogeneity=False)
    al = evaluate_with_catalog(d, conductor_id="al_1350", homogeneity=False)
    gly = evaluate_with_catalog(d, coolant_id="water_glycol_30", homogeneity=False)
    assert al["P_elec_W"] > cu["P_elec_W"]
    assert al["mass_cu_kg"] / cu["mass_cu_kg"] == pytest.approx(2705.0 / 8960.0, rel=1e-12)
    assert al["C_th_J_K"] / cu["C_th_J_K"] == pytest.approx((2705.0 * 900.0) / (8960.0 * 385.0), rel=1e-12)
    assert gly["Re"] < cu["Re"]


def test_radial_rho_does_not_change_power():
    a = evaluate_design()
    b = evaluate_design(radial_rho=True)
    assert b["P_elec_W"] == a["P_elec_W"]
    assert b["V_total_V"] == a["V_total_V"]
    assert "radial_rho_dB_over_B" in b
    assert b["radial_rho_model_grade"] == "approximate"


def test_langevin_cold_drift():
    out = langevin_ensemble(1.0e-3, 1.0e-14, -273.15, n_carriers=800, steps=500, seed=2)
    assert out["mean_v_m_s"] == pytest.approx(1.0e-3, rel=0.05)
    assert out["model_grade"] == "approximate"
    assert "not a sample of the magnet" in out["note"].lower() or "Not a sample" in out["note"]


def test_tsat_near_100c():
    assert saturation_temperature_C(101325.0) == pytest.approx(100.0, abs=0.5)


def test_unknown_cli_flag():
    with pytest.raises(SystemExit):
        main(["evaluate", "--not-a-real-flag"])


def test_emulate_seed_required_and_percentiles():
    with pytest.raises(ValueError):
        emulate(BitterDesign(), realizations=4, seed=None, homogeneity=False)
    out = emulate(BitterDesign(), realizations=8, seed=1, homogeneity=False)
    assert out["label"].startswith("mesoscopic emulation")
    assert "molecular dynamics" in out["label"]
    assert out["model_grade"]["hooge_1f"] == "speculative"
    assert out["model_grade"]["johnson"] == "established"
    p = out["percentiles"]["T_hot_C"]
    assert p["p05"] <= p["p50"] <= p["p95"]
    assert out["hooge"] is None
    assert "not the MRI noise floor" in out["johnson"]["note"]
    again = emulate(BitterDesign(), realizations=8, seed=1, homogeneity=False)
    assert again["percentiles"]["V_V"]["p50"] == out["percentiles"]["V_V"]["p50"]
