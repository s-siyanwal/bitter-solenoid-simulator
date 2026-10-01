"""Multivariable constrained design optimisation (PDF section
"Multivariable Calculus and Constrained Optimization").

Design vector  x = [R1, R2-R1, L, d_plate, D_hole, v_flow, pitch_factor]
(The PDF uses x = [R1, R2, L, V0, v_flow]; here the current/voltage is not a
free variable because it is fixed analytically by the B0 = 0.5 T equality
(E4), and plate thickness, hole diameter and hole pitch factor are added as cooling/geometry
variables.  Supply voltage V0 then follows from the design.)

Objective (E32): minimise P_elec + P_pump.  (PDF: P_elec only; adding the pump
power keeps the optimiser from choosing unlimited flow - ASSUMPTION.)
Constraints:
  g1  Bz(0) = 0.5 T                 (enforced exactly through the current)
  g2  T_hot <= 85 C                 (PDF T_limit)
  g3  V_total <= 8 V                (PDF supply limit)
  g4  R1 >= 0.05 m                  (PDF; lower bound)
  g5  homogeneity <= 100 ppm over a 30 mm DSV         (ASSUMPTION)
  g6  dp <= 5 bar                                     (ASSUMPTION)
  g7  Re >= 1e4 (Dittus-Boelter validity)             (ASSUMPTION)
Bounds follow the PDF code where given (R2 in [0.1, 0.3] m, L in [0.4, 1.5] m,
v in [1, 5] m/s); others are ASSUMPTIONS.

Global search: scipy.optimize.differential_evolution (SciPy >= 0.15) with
quadratic penalties, then SLSQP polish with explicit inequality constraints
(both available in SciPy 1.1).
"""
import numpy as np
from scipy.optimize import differential_evolution, minimize
from .design import BitterDesign, evaluate_design

BOUNDS = [(0.05, 0.10),     # R1 [m]
          (0.05, 0.25),     # R2 - R1 [m] (keeps R2 within the PDF's 0.1-0.3 m)
          (0.40, 1.50),     # L [m]
          (0.5e-3, 6.0e-3), # plate thickness [m]
          (1.0e-3, 6.0e-3), # hole diameter [m]
          (1.0, 5.0),       # water velocity [m/s]
          (1.5, 8.0)]       # hole pitch / hole diameter

LIMITS = dict(T_limit=85.0, V_max=8.0, ppm_max=100.0, dp_max=5e5, Re_min=1e4, R2_max=0.30)


def x_to_design(x, **fixed):
    R1, t, L, dpl, Dh, v, pf = x
    return BitterDesign(R1=R1, R2=R1 + t, L=L, d_plate=dpl, D_hole=Dh, v_flow=v,
                        pitch_factor=pf, **fixed)


def constraints(res, lim=LIMITS):
    """All g_i >= 0 for feasibility (normalised)."""
    return np.array([
        (lim["T_limit"] - res["T_hot_C"]) / lim["T_limit"],
        (lim["V_max"] - res["V_total_V"]) / lim["V_max"],
        (lim["ppm_max"] - res["homogeneity_ppm"]) / lim["ppm_max"],
        (lim["dp_max"] - res["dp_Pa"]) / lim["dp_max"],
        (res["Re"] - lim["Re_min"]) / lim["Re_min"],
        (lim["R2_max"] - res["R2"]) / lim["R2_max"],
    ])


def _eval(x, fixed):
    return evaluate_design(x_to_design(x, **fixed))


def penalised(x, fixed=None):
    res = _eval(x, fixed or {})
    g = constraints(res)
    viol = np.minimum(g, 0.0)
    return res["P_total_W"] / 1e3 + 1e4 * float(np.sum(viol ** 2)) + 1e2 * float(np.sum(-viol))


def optimise(seed=1, maxiter=60, popsize=15, polish=True, fixed=None, verbose=False):
    fixed = fixed or {}
    de = differential_evolution(penalised, BOUNDS, args=(fixed,), seed=seed, maxiter=maxiter,
                                popsize=popsize, tol=1e-8, polish=False, disp=verbose)
    x_best = de.x
    method = "differential_evolution"
    if polish:
        cons = [{"type": "ineq", "fun": (lambda x, i=i: constraints(_eval(x, fixed))[i])} for i in range(6)]
        sl = minimize(lambda x: _eval(x, fixed)["P_total_W"] / 1e3, de.x, method="SLSQP",
                      bounds=BOUNDS, constraints=cons, options={"maxiter": 200, "ftol": 1e-10})
        if sl.success and np.all(constraints(_eval(sl.x, fixed)) >= -1e-9) and \
                _eval(sl.x, fixed)["P_total_W"] < _eval(de.x, fixed)["P_total_W"]:
            x_best = sl.x
            method = "differential_evolution + SLSQP polish"
    res = _eval(x_best, fixed)
    return {"x": x_best, "result": res, "constraints": constraints(res), "method": method,
            "de_nfev": int(de.nfev), "de_fun": float(de.fun)}
