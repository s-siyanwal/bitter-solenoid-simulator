"""Helical current path (segment engine) and the head-preset optimiser outputs."""
import json, math, os
import numpy as np
import pytest
from bittersim import helical as Hx
from bittersim.biot_savart import biot_savart

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = (0.05, 0.30, 1.10, 2.3e5)


def _chains_closed(s0, s1, I):
    # every segment end is the start of the next segment of the same filament, except at breaks;
    # a closed circuit has no dangling ends: the multiset of starts equals the multiset of ends
    a = np.round(np.vstack([s0, s1]) * 1e9).astype(np.int64)
    starts = set(map(tuple, a[:len(s0)]))
    ends = set(map(tuple, a[len(s0):]))
    return starts == ends


@pytest.mark.parametrize("kind", ["uniform", "aligned", "rotating"])
def test_paths_with_bus_are_closed_circuits(kind):
    s0, s1, I = Hx.variant_segments(kind, *G, n_turns=12, radial_steps=2, seg_deg=10.0)
    assert _chains_closed(s0, s1, I)


def test_staircase_turn_count_and_scale():
    (s0, s1, I), start, end = Hx.staircase_path(*G, n_turns=12, overlap_deg=30.0, slit_advance_deg=30.0, radial_steps=1)
    assert end[0] == pytest.approx(11 * math.radians(30) + 2 * math.pi - math.radians(30))
    pts = np.array([[0.0, 0.0, 0.0], [0.005, 0.0, 0.0]])
    dB, sc = Hx.delta_field("rotating", *G, n_turns=24, pts=pts, radial_steps=2, seg_deg=10.0)
    assert sc == pytest.approx(1.0 / (1.0 + 30.0 / 360.0), rel=0.05)
    assert abs(dB[0, 2]) < 1e-9                       # isocentre Bz held fixed by construction


def test_ideal_coil_has_no_tesseral_terms():
    a = Hx.analyse(None, *G, n_turns=100, dsv=0.03, nmax=4)
    assert all(abs(v) < 0.01 for v in a["tesseral_ppm"].values())
    assert a["B_perp_helical_max_uT"] == 0.0


def test_helical_json_consistency():
    d = json.load(open(os.path.join(ROOT, "results", "helical.json")))
    assert not d["quick"]
    for key in ("30mm", "40mm"):
        ideal = d["dsv"][key]["ideal"]
        assert ideal["ladder_ppm"][1][1] < ideal["ladder_ppm"][0][1]
        for kind in ("uniform", "aligned", "rotating"):
            a = d["dsv"][key][kind]
            # Z2/Z4 pairs cannot remove tesseral errors: the shimmed helical coil is never better than ideal
            assert a["ladder_ppm"][1][1] >= ideal["ladder_ppm"][1][1] - 1e-6
    assert "X" in d["dsv"]["30mm"]["aligned"]["needed_shims"]


def test_head_json_consistency():
    d = json.load(open(os.path.join(ROOT, "results", "head.json")))
    assert not d["quick"]
    for r in d["runs"]["power_8V"]:
        assert r["feasible"] and r["V_total_V"] <= 8.0 + 1e-9 and r["ppm_shimmed"] <= 10.0 + 1e-6
    for r in d["runs"]["min_voltage"]:
        assert r["feasible"]
