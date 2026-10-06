"""Particle emulation: serial == parallel bit for bit, and convergence to the continuum."""
import json, math, os
import numpy as np
import pytest
from bittersim import particles as Pt
from bittersim.fields import B_bitter_axis
from bittersim.thermal import annulus_conduction_dT
from bittersim.constants import RHO_CU_20, ALPHA_CU, K_CU, MU_0

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R1, R2, L, C = 0.05, 0.30, 1.10, 2.3e5


def test_serial_and_parallel_are_bit_identical():
    ax, dsv = Pt.axis_and_dsv_points(L, 0.03, n_axis=5, n_dsv=20)
    pts = np.vstack([ax, dsv])
    assert np.array_equal(Pt.element_field(pts, R1, R2, L, C, 512, False), Pt.element_field(pts, R1, R2, L, C, 512, True))
    m0 = Pt.element_moments(R1, R2, L, C, 1 << 14, RHO_CU_20, 0.95, False)
    m1 = Pt.element_moments(R1, R2, L, C, 1 << 14, RHO_CU_20, 0.95, True)
    assert m0 == m1
    assert Pt.green_kubo_rho(20.0, 3000, 5, False) == Pt.green_kubo_rho(20.0, 3000, 5, True)
    assert Pt.walk_dT(1e6, 0.003, 0.009, K_CU, 500, 6, False) == Pt.walk_dT(1e6, 0.003, 0.009, K_CU, 500, 6, True)


def test_isocentre_and_axis_converge_to_continuum():
    E4 = MU_0 * C * (math.asinh(L / (2 * R1)) - math.asinh(L / (2 * R2)))
    z = np.array([0.0, 0.1, 0.2])
    pts = np.column_stack([0 * z, 0 * z, z])
    b = Pt.element_field(pts, R1, R2, L, C, 1 << 16, True)
    assert abs(b[0, 2] / E4 - 1) < 1e-10          # zero-variance importance sampling at the centre
    assert np.max(np.abs(b[:, 2] / B_bitter_axis(R1, R2, L, C, z) - 1)) < 1e-3
    assert np.max(np.abs(b[:, :2])) < 1e-9


def test_moments_match_NI_and_power():
    lam, rho = 0.95, RHO_CU_20
    m = Pt.element_moments(R1, R2, L, C, 1 << 18, rho, lam, True)
    lnr = math.log(R2 / R1)
    assert abs(m["NI"] / (C * L * lnr) - 1) < 1e-5
    assert abs(m["P"] / (2 * math.pi * rho * C ** 2 * L * lnr / lam) - 1) < 1e-5


def test_bloch_gruneisen_calibration_and_slope():
    assert abs(Pt.bloch_gruneisen_rho(20.0) / RHO_CU_20 - 1) < 1e-9
    slope = (Pt.bloch_gruneisen_rho(20.5) - Pt.bloch_gruneisen_rho(19.5)) / RHO_CU_20
    assert abs(slope / ALPHA_CU - 1) < 0.03       # model difference, ~1 %


def test_carrier_resistivity_is_unbiased():
    r, rt, se = Pt.green_kubo_rho(40.0, 40000, 3, True)
    assert abs(r / rt - 1) < 4 * se


def test_walkers_match_E18():
    a, b, q = 0.003, 0.009, 5e6
    w, se = Pt.walk_dT(q, a, b, K_CU, 8000, 7, True, dt_rel=1e-3)
    ref = annulus_conduction_dT(q, a, b, K_CU)
    assert abs(w - ref) < 4 * se + 0.01 * ref     # 1 % allowance for time-step bias at dt_rel=1e-3


def test_particles_json_is_consistent():
    d = json.load(open(os.path.join(ROOT, "results", "particles.json"), encoding="utf-8"))
    assert not d["quick"]
    for k, v in d["benchmark"].items():
        assert v["bit_identical"], k
    rows = dict((r["key"], r) for r in d["comparison"])
    for k in ("B0", "NI", "R", "P", "V"):
        assert abs(rows[k]["rel_err"]) < 1e-5, k
    assert abs(rows["dT_cond"]["rel_err"]) < 4 * rows["dT_cond"]["stderr_rel"] + 2e-3
    assert abs(rows["ppm"]["rel_err"]) < 0.02


def test_tolerance_mc_identical_across_workers():
    from bittersim import BitterDesign
    from bittersim.emulation import emulate
    a = emulate(BitterDesign(), realizations=6, seed=3, workers=1, homogeneity=False)
    b = emulate(BitterDesign(), realizations=6, seed=3, workers=2, homogeneity=False)
    assert json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True, default=str)
