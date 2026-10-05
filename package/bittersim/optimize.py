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
  g7  R2 <= 0.30 m
Soft (not a hard cut — BIR p.4–5 practice):
  Re floor: correlation-validity penalty below Re_min (default 5500).
  The old hard Re >= 1e4 constraint is removed; Re_min is a LIMITS parameter.

Pass friction_multiplier, R_c_ohm, hole_layout_mode, etc. via `fixed=`.
Pass Re_min / Re_penalty_weight via `limits=`.

Global search: scipy.optimize.differential_evolution (SciPy >= 0.15) with
quadratic penalties, then SLSQP polish with explicit inequality constraints
(both available in SciPy 1.1).
"""
import numpy as np
from scipy.optimize import differential_evolution, minimize
from .design import BitterDesign, evaluate_design
from .constants import RE_MIN_CORRELATION

BOUNDS = [(0.05, 0.10),     # R1 [m]
          (0.05, 0.25),     # R2 - R1 [m] (keeps R2 within the PDF's 0.1-0.3 m)
          (0.40, 1.50),     # L [m]
          (0.5e-3, 6.0e-3), # plate thickness [m]
          (1.0e-3, 6.0e-3), # hole diameter [m]
          (1.0, 5.0),       # water velocity [m/s]
          (1.5, 8.0)]       # hole pitch / hole diameter

LIMITS = dict(T_limit=85.0, V_max=8.0, ppm_max=100.0, dp_max=5e5,
              Re_min=RE_MIN_CORRELATION, Re_penalty_weight=50.0, R2_max=0.30)


def x_to_design(x, **fixed):
    R1, t, L, dpl, Dh, v, pf = x
    return BitterDesign(R1=R1, R2=R1 + t, L=L, d_plate=dpl, D_hole=Dh, v_flow=v,
                        pitch_factor=pf, **fixed)


def constraints(res, lim=None):
    """Hard inequality constraints g_i >= 0 for feasibility (normalised).

    Re is intentionally absent: correlation validity is a soft penalty
    (see correlation_penalty), matching BIR practice rather than a hard floor.
    """
    lim = lim or LIMITS
    return np.array([
        (lim["T_limit"] - res["T_hot_C"]) / lim["T_limit"],
        (lim["V_max"] - res["V_total_V"]) / lim["V_max"],
        (lim["ppm_max"] - res["homogeneity_ppm"]) / lim["ppm_max"],
        (lim["dp_max"] - res["dp_Pa"]) / lim["dp_max"],
        (lim["R2_max"] - res["R2"]) / lim["R2_max"],
    ])


def correlation_penalty(res, lim=None):
    """Soft Re floor: quadratic penalty when Re < Re_min (BIR correlation FOS)."""
    lim = lim or LIMITS
    re_min = float(lim["Re_min"])
    if re_min <= 0.0 or res["Re"] >= re_min:
        return 0.0
    w = float(lim.get("Re_penalty_weight", 50.0))
    return w * ((re_min - res["Re"]) / re_min) ** 2


def _eval(x, fixed):
    return evaluate_design(x_to_design(x, **fixed))


def penalised(x, fixed=None, lim=None):
    lim = lim or LIMITS
    res = _eval(x, fixed or {})
    g = constraints(res, lim)
    viol = np.minimum(g, 0.0)
    return (res["P_total_W"] / 1e3
            + 1e4 * float(np.sum(viol ** 2))
            + 1e2 * float(np.sum(-viol))
            + correlation_penalty(res, lim))


def optimise(seed=1, maxiter=60, popsize=15, polish=True, fixed=None, verbose=False,
             limits=None):
    """Run DE (+ optional SLSQP polish).

    fixed : dict passed to BitterDesign (e.g. R_c_ohm, friction_multiplier,
            hole_layout_mode). Use R_c_ohm > 0 to carry a contact-resistance budget.
    limits : override LIMITS entries (T_limit, V_max, Re_min, ...).
    """
    fixed = fixed or {}
    lim = dict(LIMITS)
    if limits:
        lim.update(limits)
    de = differential_evolution(
        lambda x: penalised(x, fixed, lim), BOUNDS,
        seed=seed, maxiter=maxiter, popsize=popsize, tol=1e-8, polish=False, disp=verbose)
    x_best = de.x
    method = "differential_evolution"
    n_hard = 5
    if polish:
        cons = [{"type": "ineq",
                 "fun": (lambda x, i=i: constraints(_eval(x, fixed), lim)[i])}
                for i in range(n_hard)]
        sl = minimize(lambda x: _eval(x, fixed)["P_total_W"] / 1e3
                      + correlation_penalty(_eval(x, fixed), lim),
                      de.x, method="SLSQP", bounds=BOUNDS, constraints=cons,
                      options={"maxiter": 200, "ftol": 1e-10})
        if sl.success and np.all(constraints(_eval(sl.x, fixed), lim) >= -1e-9) and \
                _eval(sl.x, fixed)["P_total_W"] < _eval(de.x, fixed)["P_total_W"]:
            x_best = sl.x
            method = "differential_evolution + SLSQP polish"
    res = _eval(x_best, fixed)
    return {"x": x_best, "result": res, "constraints": constraints(res, lim),
            "correlation_penalty": correlation_penalty(res, lim),
            "limits": lim, "method": method,
            "de_nfev": int(de.nfev), "de_fun": float(de.fun)}
