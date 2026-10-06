"""Identities for the explicit layer-stack module (bittersim.stack)."""
import math

import numpy as np
import pytest

from bittersim import fields, stack
from bittersim.constants import MU_0, RHO_CU_20

R1, R2, T = 0.02, 0.05, 1e-3


def test_full_layer_equals_e4():
    z = np.linspace(-0.05, 0.05, 21)
    ly = stack.Layer(R1, R2, 0.01, T)
    C = 1.0 / (T * math.log(R2 / R1))
    assert np.allclose(stack.axis_field([ly], z), fields.B_bitter_axis(R1, R2, T, C, z - 0.01), rtol=1e-14)


def test_two_half_layers_equal_one_full_layer():
    z = np.linspace(-0.05, 0.05, 21)
    half = [stack.Layer(R1, R2, 0.0, T, math.pi), stack.Layer(R1, R2, 0.0, T, math.pi)]
    assert np.allclose(stack.axis_field(half, z), stack.axis_field([stack.Layer(R1, R2, 0.0, T)], z), rtol=1e-14)


def test_mirrored_pairs_helmholtz_and_anti():
    top = [stack.Layer(R1, R2, 0.02, T)]
    anti = stack.mirrored(top, sign=-1.0)
    assert abs(stack.axis_field(anti, 0.0)[0]) < 1e-18
    co = stack.mirrored(top)
    assert stack.axis_field(co, 0.0)[0] == pytest.approx(2 * stack.axis_field(top, 0.0)[0])


def test_uniform_profile_equals_e3():
    ly = stack.Layer(R1, R2, 0.0, 0.02, profile="uniform")
    j = 1.0 / (0.02 * (R2 - R1))
    assert stack.axis_field([ly], 0.03)[0] == pytest.approx(fields.B_thick_axis(R1, R2, 0.02, j, 0.03), rel=1e-14)


def test_bitter_layer_resistance_matches_plate_formula():
    # N full plates: R = 2 pi rho N / (t ln(R2/R1)) (same law as the continuum E30/E31 path)
    layers = [stack.Layer(R1, R2, k * 1.25e-3, T) for k in range(31)]
    R = stack.stack_resistance(layers)
    assert R["R_layers"] == pytest.approx(2 * math.pi * RHO_CU_20 * 31 / (T * math.log(R2 / R1)), rel=1e-12)
    assert R["R_joints"] == 0.0 and R["R_leads"] == 0.0


def test_half_layer_resistance_is_half():
    full = stack.layer_resistance(stack.Layer(R1, R2, 0.0, T))
    assert stack.layer_resistance(stack.Layer(R1, R2, 0.0, T, math.pi)) == pytest.approx(full / 2)


def test_joint_and_lead_terms_are_explicit_and_nonnegative():
    layers = [stack.Layer(R1, R2, 0.0, T)]
    base = stack.stack_resistance(layers)["R_total"]
    R = stack.stack_resistance(layers, n_joints=10, r_joint=5e-6, r_leads=1e-4)
    assert R["R_total"] == pytest.approx(base + 5e-5 + 1e-4)
    with pytest.raises(ValueError):
        stack.stack_resistance(layers, r_joint=-1e-6)


def test_implied_extra_resistance_is_diagnostic_only():
    layers = [stack.Layer(R1, R2, 0.0, T)]
    Rl = stack.stack_resistance(layers)["R_layers"]
    d = stack.implied_extra_resistance(Rl + 2e-4, layers, n_joints=4, contact_area=1e-4)
    assert d["R_gap"] == pytest.approx(2e-4)
    assert d["per_joint"] == pytest.approx(5e-5)
    assert d["specific_contact_ohm_m2"] == pytest.approx(5e-9)
    assert stack.stack_resistance(layers)["R_layers"] == pytest.approx(Rl)   # model unchanged
