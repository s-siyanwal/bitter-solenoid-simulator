"""Checks taken from the arXiv papers in the user's ``arxiv _papers`` folder.

Two kinds of test, kept apart:

* code_check -- the simulator must reproduce a closed form printed in a paper
  (Sabulsky et al. arXiv:1309.5330 Eq. 6; Kobelev arXiv:1610.06607 Eq. 2.6/3.1).
* measured   -- the bare solver against a published measurement
  (Hataway et al. arXiv:2607.02813, Claw-ZS Fig. 4a). Tolerances lock the
  current bare-solver residual; they are not a fit and nothing is tuned.
"""
import math
import os

import numpy as np
import pytest

from bittersim import fields, inductance, mechanics
from bittersim.constants import MU_0
from bittersim.elliptic import field_from_loops

DATA = os.path.join(os.path.dirname(__file__), "data")


# ------------------------------------------------- Sabulsky Eq. 6 (code_check)
def _sabulsky_eq6(a1, a2, Ct, z):
    """Single thin Bitter arc, J = C/r: Bz = mu0 C t/2 [1/sqrt(a1^2+z^2) - 1/sqrt(a2^2+z^2)]."""
    return 0.5 * MU_0 * Ct * (1.0 / np.sqrt(a1 ** 2 + z ** 2) - 1.0 / np.sqrt(a2 ** 2 + z ** 2))


def test_e4_thin_plate_limit_matches_sabulsky_eq6():
    # Sabulsky arc radii a1 = 31.75 mm, a2 = 50.80 mm; plate current 400 A.
    a1, a2, I = 31.75e-3, 50.80e-3, 400.0
    Ct = I / math.log(a2 / a1)
    z = np.linspace(-0.05, 0.05, 41)
    ref = _sabulsky_eq6(a1, a2, Ct, z)
    errs = []
    for t in (1e-3, 1e-4, 1e-5):
        B = fields.B_bitter_axis(a1, a2, t, Ct / t, z)
        errs.append(np.max(np.abs(B - ref)) / np.max(ref))
    assert errs[-1] < 1e-8
    # finite-thickness error falls as t^2
    assert errs[0] / errs[1] == pytest.approx(100.0, rel=0.05)


def test_sabulsky_eq6_equals_radial_loop_integral():
    a1, a2, I = 31.75e-3, 50.80e-3, 400.0
    Ct = I / math.log(a2 / a1)
    x, w = np.polynomial.legendre.leggauss(64)
    r = 0.5 * (a2 - a1) * x + 0.5 * (a2 + a1)
    Ir = Ct / r * w * 0.5 * (a2 - a1)          # current of each loop
    z = np.array([0.0, 0.014, 0.03])
    _, Bz = field_from_loops(np.zeros_like(z), z, r, np.zeros_like(r), Ir)
    assert np.allclose(Bz, _sabulsky_eq6(a1, a2, Ct, z), rtol=1e-10)


# ---------------------------------------------- Kobelev long stack (code_check)
def test_bitter_center_tends_to_long_stack_limit():
    # Kobelev Eq. 2.6 with J = C/r: inside the bore B = mu0 C ln(R2/R1) for L -> inf.
    R1, R2, C = 0.05, 0.15, 1e5
    lim = MU_0 * C * math.log(R2 / R1)
    errs = [abs(fields.B_bitter_center(R1, R2, L, C) / lim - 1.0) for L in (4.0, 8.0, 16.0)]
    assert errs[-1] < 1e-4
    # end correction is O((R/L)^2)
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.05)
    assert errs[1] / errs[2] == pytest.approx(4.0, rel=0.05)


def test_hoop_stress_field_in_winding_matches_long_stack_profile():
    # In the winding of a long stack Bz(r) = mu0 C ln(R2/r) (Kobelev Eq. 3.1 form),
    # so the thin-ring hoop stress is sigma = C Bz = mu0 C^2 ln(R2/r).
    R1, R2, L, C = 0.05, 0.15, 4.0, 1e5
    r, Bz, sigma = mechanics.hoop_stress_profile(R1, R2, L, C, n=21)
    B0 = MU_0 * C * math.log(R2 / R1)
    ref = MU_0 * C * np.log(R2 / r)
    assert np.max(np.abs(Bz - ref)) / B0 < 5e-3
    assert np.allclose(sigma, C * Bz)
    assert np.all(np.diff(sigma) < 0)         # peak stress at the bore


# ------------------------------------------- Claw-ZS Fig. 4a (measured, 10 A)
R_IN, R_OUT, I_CLAW = 13.7e-3, 35.0e-3, 10.0


def _claw_layers():
    a = np.genfromtxt(os.path.join(DATA, "claw_zs_layers.csv"), delimiter=",",
                      comments="#", skip_header=5, filling_values=0.0)
    t = a[:, 1] * 25.4e-3
    gap = a[:, 2] * 25.4e-3
    z_bottom = np.concatenate([[0.0], np.cumsum(t + gap)[:-1]])
    return t, z_bottom + t / 2.0, z_bottom[-1] + t[-1]


def _claw_axis(z):
    t, zc, _ = _claw_layers()
    B = np.zeros_like(z)
    for ti, zi in zip(t, zc):
        C = I_CLAW / (ti * math.log(R_OUT / R_IN))
        B += 0.5 * fields.B_bitter_axis(R_IN, R_OUT, ti, C, z - zi)   # 180 deg half-layer
    return B


def _claw_measured():
    m = np.genfromtxt(os.path.join(DATA, "claw_zs_fig4a_Bz_10A.csv"), delimiter=",", names=True)
    return m["z_mm"] * 1e-3, m["Bz_T"]


def test_claw_layer_table_matches_paper():
    t, _, length = _claw_layers()
    assert t.size == 106                       # paper: 106 "C"-shaped layers
    assert length == pytest.approx(0.320, rel=0.01)   # paper: l = 320 mm
    assert t.sum() == pytest.approx(197.485e-3, rel=1e-9)


def test_claw_peak_field_against_measurement():
    z, Bm = _claw_measured()
    zz = np.linspace(0.0, 0.1, 2001)
    Bs = _claw_axis(zz)
    k = np.argmax(Bs)
    assert Bm.max() == pytest.approx(2.760e-3, abs=1e-5)
    rel = (Bm.max() - Bs[k]) / Bs[k]
    assert abs(rel) < 0.02                     # bare solver: -1.1 %
    assert abs(zz[k] - z[np.argmax(Bm)]) < 5e-3


def test_claw_profile_rms_against_measurement():
    z, Bm = _claw_measured()
    res = Bm - _claw_axis(z)
    assert z.size == 46
    assert np.sqrt(np.mean(res ** 2)) < 0.15e-3   # bare solver: 0.100 mT
    assert np.max(np.abs(res)) < 0.25e-3          # bare solver: 0.204 mT, at the z = 0 edge


# ------------------- EPFL bulk-machined spiral coil (measured + cross-code)
# Haeusler et al., "Compact bulk-machined electromagnets for quantum gas
# experiments", SciPost Phys. 6, 048 (2019), arXiv:1901.08791. Flat spiral,
# 31 turns at 1.3 mm radial pitch, 0.92 mm Cu width, 22 mm tall: the current
# density is uniform over the winding (E3), not Bitter 1/r.
EPFL = dict(R1=32e-3, R2=72e-3, H=22e-3, N=31)


def test_epfl_coil_field_at_atoms_against_hall_measurement():
    # Sec. 4.1: Hall probe (Lakeshore 425), 52.2 mm below the mean coil position.
    j = EPFL["N"] * 1.0 / ((EPFL["R2"] - EPFL["R1"]) * EPFL["H"])
    B = fields.B_thick_axis(EPFL["R1"], EPFL["R2"], EPFL["H"], j, 52.2e-3)
    measured = 1.28e-4                         # T/A
    assert abs((measured - B) / B) < 0.02      # bare solver: -0.3 %


def test_epfl_coil_inductance_against_radia():
    # Cross-code, not a measurement: authors' Radia model gives 94 uH (measured 116(2) uH).
    L = inductance.coil_inductance(EPFL["R1"], EPFL["R2"], EPFL["H"], EPFL["N"],
                                   profile="uniform", nr=16)
    assert L == pytest.approx(94e-6, rel=0.05)  # solver: 91.6 uH
