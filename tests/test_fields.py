import math
import numpy as np
import pytest
from bittersim import fields, biot_savart as bs, validation as v
from bittersim.constants import MU_0
from bittersim.elliptic import ellipke_agm, loop_field_numba
from scipy.special import ellipk, ellipe


def test_agm_matches_scipy():
    for m in np.linspace(0, 0.9999, 200):
        K, E = ellipke_agm(m)
        assert abs(K - ellipk(m)) < 1e-13 * ellipk(m)
        assert abs(E - ellipe(m)) < 1e-13 * ellipe(m)


def test_numba_loop_matches_scipy_loop():
    r = v.v_elliptic()
    assert r["numba_vs_scipy_loop_max_rel"] < 1e-10


def test_single_loop_on_axis():
    a, I = 0.1, 1000.0
    for z in (-0.2, 0.0, 0.05, 0.3):
        exact = MU_0 * I * a ** 2 / (2 * (a ** 2 + z ** 2) ** 1.5)
        assert loop_field_numba(a, I, 0.0, z)[1] == pytest.approx(exact, rel=1e-14)
        assert loop_field_numba(a, I, 1e-7, z)[1] == pytest.approx(exact, rel=1e-10)


def test_polygon_loop_converges_second_order():
    rows = v.v_single_loop()["segments_center_error"]
    for (n1, e1), (n2, e2) in zip(rows[:-1], rows[1:]):
        assert 3.8 < e1 / e2 < 4.2


@pytest.mark.parametrize("ratio,tol", [(100, 3e-4), (1000, 3e-6)])
def test_long_solenoid_limit(ratio, tol):
    R, n, I = 0.05, 1000.0, 10.0
    coil = fields.CoilLoops(R, R, ratio * R, n * I, "thin", nz=16, z_panels=ratio // 2)
    B = coil.field(np.array([0.0]), np.array([0.0]))[1][0]
    assert abs(B - MU_0 * n * I) / (MU_0 * n * I) < tol


def test_thin_finite_formula_vs_loops():
    R, L, n, I = 0.05, 0.3, 500.0, 20.0
    coil = fields.CoilLoops(R, R, L, n * I, "thin", nz=16, z_panels=4)
    z = np.linspace(-0.4, 0.4, 41)
    assert np.allclose(coil.field(np.zeros_like(z), z)[1], fields.B_thin_axis(R, L, n, I, z), rtol=1e-11)


def test_thick_closed_form_and_bitter_closed_form():
    r = v.v_thick_closed_form()
    assert r["uniform_loops_vs_E3_max_rel"] < 1e-10
    assert r["E3_center_vs_PDF_form_rel"] < 1e-14
    assert r["bitter_loops_vs_E4_max_rel"] < 1e-10


def test_bitter_center_equals_axis_formula():
    assert fields.B_bitter_center(0.05, 0.2, 0.7, 4e5) == pytest.approx(
        float(fields.B_bitter_axis(0.05, 0.2, 0.7, 4e5, 0.0)), rel=1e-14)


def test_segment_engine_pdf_test1_and_richardson():
    r = v.v_segment_convergence()
    assert abs(r["richardson_p"] - 2.0) < 0.01
    assert all(e < 0.00195 for _, e, _ in r["rows_vs_E3"])          # PDF 0.195 % criterion
    assert abs(r["B_richardson"] - r["B_filament_limit"]) / r["B_filament_limit"] < 1e-8


def test_segment_engine_bitter_pdf_test2():
    R1, R2, L, C = 0.05, 0.15, 0.4, 3e5
    s0, s1, c = bs.generate_coil_segments(R1, R2, L, C, "bitter", 40, 16, 180, helical=False)
    B = bs.biot_savart(np.zeros((1, 3)), s0, s1, c)[0, 2]
    assert abs(B - fields.B_bitter_center(R1, R2, L, C)) / fields.B_bitter_center(R1, R2, L, C) < 2e-3


def test_chunking_invariant():
    s0, s1, c = bs.generate_coil_segments(0.05, 0.1, 0.2, 1e5, "bitter", 5, 3, 24)
    obs = np.random.RandomState(1).uniform(-0.03, 0.03, (1234, 3))
    assert np.allclose(bs.biot_savart(obs, s0, s1, c, chunk=100), bs.biot_savart(obs, s0, s1, c, chunk=5000),
                       rtol=0, atol=1e-15)


def test_off_axis_segments_vs_loops():
    R1, R2, L, j = 0.05, 0.1, 0.2, 1e7
    s0, s1, c = bs.generate_coil_segments(R1, R2, L, j, "uniform", 20, 5, 720, helical=False)
    pts = np.array([[0.02, 0.0, 0.01], [0.0, 0.03, -0.05], [0.12, 0.0, 0.0]])
    Bs = bs.biot_savart(pts, s0, s1, c)
    dz = L / 20; dr = (R2 - R1) / 5
    rc = R1 + (np.arange(5) + 0.5) * dr; zc = -L / 2 + (np.arange(20) + 0.5) * dz
    RR, ZZ = np.meshgrid(rc, zc, indexing="ij")
    from bittersim.elliptic import field_from_loops
    rho = np.hypot(pts[:, 0], pts[:, 1])
    Br, Bz = field_from_loops(rho, pts[:, 2].copy(), RR.ravel().copy(), ZZ.ravel().copy(), np.full(RR.size, j * dr * dz))
    assert np.allclose(Bs[:, 2], Bz, rtol=2e-5)
    assert np.allclose(np.hypot(Bs[:, 0], Bs[:, 1]), np.abs(Br), rtol=2e-4, atol=1e-9)


def test_homogeneity_converged():
    from bittersim.design import evaluate_design
    res = evaluate_design()
    rows = v.v_homogeneity_convergence(res)
    assert abs(rows[1]["ppm"] - rows[-1]["ppm"]) < 1e-3
    assert abs(rows[-1]["B0"] - 0.5) < 1e-8
