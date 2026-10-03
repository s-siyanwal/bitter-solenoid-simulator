import numpy as np
from bittersim.optimize import optimise, constraints, x_to_design, BOUNDS
from bittersim import evaluate_design


def test_optimiser_returns_feasible_design():
    o = optimise(seed=3, maxiter=8, popsize=8, polish=False)
    assert np.all(o["constraints"] > -0.05)
    assert abs(o["result"]["B0_numeric_T"] - 0.5) < 1e-6
    lo = np.array([b[0] for b in BOUNDS]); hi = np.array([b[1] for b in BOUNDS])
    assert np.all(o["x"] >= lo - 1e-12) and np.all(o["x"] <= hi + 1e-12)
