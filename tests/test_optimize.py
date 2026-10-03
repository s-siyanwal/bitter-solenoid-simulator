import json, os
import numpy as np
from bittersim.optimize import optimise, constraints, x_to_design, BOUNDS
from bittersim import evaluate_design

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_optimiser_returns_feasible_design():
    o = optimise(seed=3, maxiter=8, popsize=8, polish=False)
    assert np.all(o["constraints"] >= 0.0)          # strictly feasible, no tolerance
    assert abs(o["result"]["B0_numeric_T"] - 0.5) < 1e-6
    lo = np.array([b[0] for b in BOUNDS]); hi = np.array([b[1] for b in BOUNDS])
    assert np.all(o["x"] >= lo - 1e-12) and np.all(o["x"] <= hi + 1e-12)


def test_published_optimum_is_feasible_and_reproducible():
    R = json.load(open(os.path.join(ROOT, "results", "results.json")))
    x = R["optimiser"]["x"]
    res = evaluate_design(x_to_design(x))
    g = constraints(res)
    assert np.all(g >= -1e-9)
    assert res["P_elec_W"] == R["optimal"]["P_elec_W"] or abs(res["P_elec_W"] / R["optimal"]["P_elec_W"] - 1) < 1e-9
