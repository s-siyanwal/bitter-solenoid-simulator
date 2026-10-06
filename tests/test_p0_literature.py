"""P0 literature-audit fixes: contact R_c, stacked friction, soft Re, graded holes, BETA."""
import math
import json
import os
import numpy as np
import pytest
from bittersim import BitterDesign, evaluate_design
from bittersim.constants import (
    FRICTION_MULTIPLIER_STACK_MON, R_C_HAO_EXAMPLE_OHM, RE_MIN_CORRELATION,
)
from bittersim.optimize import (
    constraints, correlation_penalty, LIMITS, x_to_design, optimise,
)
from bittersim.validation import v_beta_benchmark
from bittersim import thermal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_published_path_bit_stable_with_p0_defaults():
    R = json.load(open(os.path.join(ROOT, "results", "results.json"), encoding="utf-8"))
    res = evaluate_design(x_to_design(R["optimiser"]["x"]), homogeneity=False)
    assert res["P_elec_W"] == pytest.approx(R["optimal"]["P_elec_W"], rel=1e-12)
    assert res["dp_Pa"] == pytest.approx(R["optimal"]["dp_Pa"], rel=1e-12)
    assert res["Re"] == pytest.approx(R["optimal"]["Re"], rel=1e-12)
    assert res["R_c_ohm"] == 0.0
    assert res["friction_multiplier"] == 1.0
    assert res["hole_layout_mode"] == "uniform"


def test_contact_resistance_in_R_V_P_and_headroom():
    base = evaluate_design(homogeneity=False)
    rc = 5e-6
    r = evaluate_design(BitterDesign(R_c_ohm=rc), homogeneity=False)
    n_if = r["n_interfaces"]
    assert n_if >= 1
    assert r["R_contact_ohm"] == pytest.approx(n_if * rc)
    assert r["R_total_ohm"] == pytest.approx(r["R_copper_ohm"] + r["R_contact_ohm"])
    assert r["V_total_V"] == pytest.approx(r["I_A"] * r["R_total_ohm"])
    assert r["P_elec_W"] == pytest.approx(r["P_copper_W"] + r["P_contact_W"])
    assert r["P_contact_W"] == pytest.approx(r["I_A"] ** 2 * r["R_contact_ohm"])
    assert r["dT_contact_K"] > 0.0
    assert r["T_hot_C"] > base["T_hot_C"]
    # Optimum headroom ~7.6 µΩ/interface (AUDIT M1)
    pub = json.load(open(os.path.join(ROOT, "results", "results.json"), encoding="utf-8"))
    opt = evaluate_design(x_to_design(pub["optimiser"]["x"]), homogeneity=False)
    assert opt["R_c_max_per_interface_ohm"] == pytest.approx(7.6e-6, rel=0.05)


def test_hao_scale_contact_blows_voltage_budget():
    """HAO ~65 µΩ/interface is far above the 8 V headroom of the published optimum."""
    pub = json.load(open(os.path.join(ROOT, "results", "results.json"), encoding="utf-8"))
    r = evaluate_design(x_to_design(pub["optimiser"]["x"], R_c_ohm=R_C_HAO_EXAMPLE_OHM),
                        homogeneity=False)
    assert r["V_total_V"] > 8.0


def test_stacked_friction_scales_dp_not_h():
    smooth = evaluate_design(homogeneity=False)
    stack = evaluate_design(BitterDesign(friction_multiplier=FRICTION_MULTIPLIER_STACK_MON),
                            homogeneity=False)
    assert stack["h_W_m2K"] == pytest.approx(smooth["h_W_m2K"], rel=1e-12)
    # Minor loss kept; frictional part scales by ~15 so total dp ratio is between 1 and 15
    ratio = stack["dp_Pa"] / smooth["dp_Pa"]
    assert 5.0 < ratio < FRICTION_MULTIPLIER_STACK_MON + 0.1
    assert stack["f_darcy"] == pytest.approx(
        smooth["f_darcy_smooth"] * FRICTION_MULTIPLIER_STACK_MON, rel=1e-12)


def test_soft_re_penalty_not_hard_constraint():
    assert "Re_min" in LIMITS
    assert LIMITS["Re_min"] == RE_MIN_CORRELATION
    # constraints array no longer includes Re
    res = evaluate_design(homogeneity=True)
    g = constraints(res)
    assert g.shape == (5,)
    # Below floor → positive penalty; above → zero
    low = dict(res)
    low["Re"] = 0.5 * RE_MIN_CORRELATION
    assert correlation_penalty(low) > 0.0
    high = dict(res)
    high["Re"] = 2.0 * RE_MIN_CORRELATION
    assert correlation_penalty(high) == 0.0


def test_re_min_exposed_as_parameter():
    # Custom floor via limits=
    o = optimise(seed=1, maxiter=2, popsize=4, polish=False,
                 limits={"Re_min": 1.0e4, "Re_penalty_weight": 10.0})
    assert o["limits"]["Re_min"] == 1.0e4
    assert "correlation_penalty" in o


def test_graded_layouts_montgomery_and_vinokur():
    uni = evaluate_design(BitterDesign(hole_layout_mode="uniform"), homogeneity=False)
    mon = evaluate_design(BitterDesign(hole_layout_mode="montgomery"), homogeneity=False)
    vin = evaluate_design(BitterDesign(hole_layout_mode="vinokur"), homogeneity=False)
    assert uni["hole_layout_mode"] == "uniform"
    assert mon["hole_layout_mode"] == "montgomery"
    assert vin["hole_layout_mode"] == "vinokur"
    # Montgomery: equal holes per ring
    lay = BitterDesign(hole_layout_mode="montgomery").hole_layout()
    assert len(set(lay["n_per_row"].tolist())) == 1
    # Geometric edges (spacing grows with r)
    dr = np.diff(lay["edges"])
    assert dr[-1] > dr[0]
    # Graded layouts flatten the radial hot-spot contrast vs uniform
    span_u = uni["T_hot_rows_C"].max() - uni["T_hot_rows_C"].min()
    span_m = mon["T_hot_rows_C"].max() - mon["T_hot_rows_C"].min()
    assert span_m < 0.5 * span_u
    # Docs no longer claim default is graded
    assert "NOT a graded" in " ".join(uni["warnings"])


def test_elongated_aspect_preserves_area_changes_Dh():
    round_ = BitterDesign(elongated_aspect=1.0).hole_layout()
    elong = BitterDesign(elongated_aspect=4.0).hole_layout()
    assert elong["area_hole"] == pytest.approx(round_["area_hole"], rel=1e-12)
    assert elong["D_h"] < round_["D_h"]
    assert elong["perim_hole"] > round_["perim_hole"]


def test_beta_benchmark_table_I_resistance():
    b = v_beta_benchmark()
    # BETA Eq 3 must reproduce Table I R, V, P closely (ρ(T_avg) from our linear law)
    assert b["rel_R"] < 0.03
    assert b["rel_V"] < 0.03
    assert b["rel_P"] < 0.03
    # Inductance within 25 % (geometry idealisation; Table I L = 196.86 µH)
    assert b["rel_L"] < 0.25
    # h at Table II velocities lands near the published band
    hs = [row["h"] for row in b["h_at_velocities"]]
    assert min(hs) < b["table_II"]["h_max"] * 1.5
    assert max(hs) > b["table_II"]["h_min"] * 0.5
    # Measured ΔP anchor present; stack friction raises ΔP toward/above it
    assert b["dp_measured_Pa"] == pytest.approx(7.72e3)
    assert b["dp_beta_analytic_Pa"] == pytest.approx(9.69e3)
    assert b["dp_stack15_Pa"] > b["dp_smooth_Pa"]


def test_friction_multiplier_rejects_nonpositive():
    with pytest.raises(ValueError):
        thermal.channel_flow(2.0, 5e-3, 0.8, 25.0, friction_multiplier=0.0)
