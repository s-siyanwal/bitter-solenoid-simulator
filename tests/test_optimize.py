import json, os
import numpy as np
from bittersim.optimize import optimise, constraints, x_to_design, BOUNDS
from bittersim import evaluate_design

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_optimiser_returns_feasible_design():
    o = optimise(seed=3, maxiter=25, popsize=10, polish=False)
    # tightened from -0.05; a short DE run is not exactly feasible on every SciPy version
    assert np.all(o["constraints"] >= -0.01), o["constraints"]
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


def test_polish_is_judged_on_the_penalised_objective(monkeypatch):
    """Regression test: the SLSQP polish minimises P_total + soft Re penalty, so a polished point
    must be accepted or rejected on that same objective. A point that saves raw power by dropping
    below the soft Re floor, but is worse once the penalty counts, must be rejected."""
    import types
    from bittersim import optimize as O
    x_de = np.array([0.05, 0.24, 1.10, 6.0e-3, 5.2e-3, 1.15, 8.0])
    x_sl = x_de.copy()
    x_sl[5] = 1.0                      # slower water: less pump power, Re below the soft floor
    r_de, r_sl = O._eval(x_de, {}), O._eval(x_sl, {})
    obj = lambda r: r["P_total_W"] / 1e3 + O.correlation_penalty(r)
    # preconditions of the scenario are checked, not assumed
    assert np.all(O.constraints(r_de) > 0) and np.all(O.constraints(r_sl) > 0)
    assert O.correlation_penalty(r_de) == 0.0 < O.correlation_penalty(r_sl)
    assert r_sl["P_total_W"] < r_de["P_total_W"]          # cheaper in raw power ...
    assert obj(r_sl) > obj(r_de)                           # ... but worse on the minimised objective
    monkeypatch.setattr(O, "differential_evolution", lambda *a, **k: types.SimpleNamespace(
        x=x_de, fun=O.penalised(x_de, {}), nfev=1, nit=1, success=False, message="stub"))
    monkeypatch.setattr(O, "minimize", lambda *a, **k: types.SimpleNamespace(
        x=x_sl, success=True, status=0, nit=1, message="stub"))
    o = O.optimise(seed=1, maxiter=1, popsize=5)
    assert o["method"] == "differential_evolution"
    assert np.array_equal(o["x"], x_de)
