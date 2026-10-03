"""Time-domain B0 stability: power supply, copper thermal expansion, cooling-water drift.

S1  Field at fixed current scales as 1/length for a uniform isotropic expansion of every
    coil dimension (B = mu0 * N * I * g(R1/L, R2/L) / L):  dB/B = -alpha_L dT_cu.
    `expansion_coefficient_numeric` checks this by re-evaluating E4 on a scaled geometry.
S2  Copper responds to the inlet water temperature through the lumped E19 lag:
    tau dT_cu/dt = T_in(t) + dT_rise - T_cu  (tau from thermal.lumped_time_constant).
S3  Power supply:
    current-regulated: dI/I = ripple(t) + drift(t) (spec values are inputs);
    voltage-regulated: I = V/R(T_cu), so dI/I = -alpha_rho dT_cu / (1 + alpha_rho (T_cu - 20))
    plus voltage ripple filtered by the coil L/R: |H| = 1/sqrt(1 + (w L/R)^2).
S4  Output: total dB/B in ppm and in Hz at the proton Larmor frequency (1 ppm = f_L/1e6 Hz),
    compared with a stability target (default 1 ppm, ASSUMPTION: a common EPI rule of thumb).
Every spec default below is an ASSUMPTION to be replaced with vendor data.
"""
import math
import numpy as np
from .constants import MU_0, ALPHA_CU, GAMMA_PROTON

ALPHA_L_CU = 16.5e-6     # copper linear expansion near 20 C [1/K] (ASSUMPTION: handbook value)

DEFAULT_SPEC = dict(
    mode="current",            # "current" or "voltage" regulated supply
    ripple_ppm=1.0,            # PSU ripple amplitude [ppm of I or V]          ASSUMPTION
    ripple_hz=300.0,           # ripple frequency (6-pulse rectifier on 50 Hz)  ASSUMPTION
    drift_ppm_per_h=2.0,       # PSU drift [ppm/h]                             ASSUMPTION
    water_amp_K=0.1,           # inlet water oscillation amplitude [K]          ASSUMPTION
    water_period_s=300.0,      # chiller cycling period [s]                     ASSUMPTION
    water_drift_K_per_h=0.5,   # slow inlet drift [K/h]                         ASSUMPTION
    duration_s=600.0,          # one fMRI run                                   ASSUMPTION
    target_ppm=1.0,            # EPI stability target, peak-to-peak             ASSUMPTION
)


def bitter_center_field(R1, R2, L, NI):
    """E4 at z = 0 written with the ampere-turns: C = NI / (L ln(R2/R1))."""
    C = NI / (L * math.log(R2 / R1))
    return MU_0 * C * (math.asinh(L / (2 * R1)) - math.asinh(L / (2 * R2)))


def expansion_coefficient_numeric(R1, R2, L, NI, dT=1.0, alpha_L=ALPHA_L_CU):
    """S1 check: (1/B) dB/dT at fixed current for an isotropic expansion alpha_L dT."""
    s = 1.0 + alpha_L * dT
    b0 = bitter_center_field(R1, R2, L, NI)
    b1 = bitter_center_field(R1 * s, R2 * s, L * s, NI)
    return (b1 - b0) / b0 / dT


def first_order_lag(t, u, tau, y0):
    """Exact zero-order-hold discretisation of tau y' = u - y."""
    y = np.empty_like(u)
    y[0] = y0
    a = math.exp(-(t[1] - t[0]) / tau) if np.isfinite(tau) else 1.0
    for i in range(1, len(t)):
        y[i] = a * y[i - 1] + (1 - a) * u[i - 1]
    return y


def simulate(res, tau, L_over_R, spec=None, dt=None):
    """S1-S4 time series for one design result `res` (from evaluate_design).

    The slow terms are sampled every dt; the ripple (hundreds of Hz) is added as its
    worst-case peak-to-peak 2*amplitude. Returns dict with t, ppm components, peak-to-peak values (ppm and Hz) and pass/fail."""
    sp = dict(DEFAULT_SPEC); sp.update(spec or {})
    T = sp["duration_s"]
    dt = dt or 0.1                                  # slow terms; ripple is handled analytically
    t = np.arange(0.0, T + 0.5 * dt, dt)
    w = 2 * math.pi * sp["ripple_hz"]
    f_L = GAMMA_PROTON * res["B0"] / (2 * math.pi)
    T_in = (res["T_in"] + sp["water_amp_K"] * np.sin(2 * math.pi * t / sp["water_period_s"])
            + sp["water_drift_K_per_h"] * t / 3600.0)
    rise = res["T_cu_mean_C"] - res["T_in"]
    T_cu = first_order_lag(t, T_in + rise, tau, res["T_cu_mean_C"])
    dTc = T_cu - res["T_cu_mean_C"]
    expansion = -ALPHA_L_CU * dTc * 1e6
    drift = sp["drift_ppm_per_h"] * t / 3600.0
    if sp["mode"] == "current":
        ripple_amp = sp["ripple_ppm"]
        resist = np.zeros_like(t)
    else:
        ripple_amp = sp["ripple_ppm"] / math.sqrt(1.0 + (w * L_over_R) ** 2)
        T0 = res["T_cu_mean_C"]
        resist = ((1 + ALPHA_CU * (T0 - 20.0)) / (1 + ALPHA_CU * (T_cu - 20.0)) - 1.0) * 1e6
    slow = drift + expansion + resist              # sampled at dt; ripple adds +-ripple_amp on top
    pp = lambda x: float(np.max(x) - np.min(x))
    budget = {"psu_ripple": 2.0 * ripple_amp, "psu_drift": pp(drift), "thermal_expansion": pp(expansion),
              "resistance_drift": pp(resist), "total": pp(slow) + 2.0 * ripple_amp}
    return {"t": t, "slow_ppm": slow, "ripple_amp_ppm": ripple_amp, "T_cu": T_cu, "spec": sp, "f_larmor_Hz": f_L,
            "budget_ppm": budget, "budget_Hz": {k: v * f_L / 1e6 for k, v in budget.items()},
            "meets_target": budget["total"] <= sp["target_ppm"], "tau_s": tau}


def requirements(res, tau, spec=None):
    """Largest disturbances that alone use the whole target (peak-to-peak).

    Water oscillation of amplitude A at period P reaches copper attenuated by 1/sqrt(1+(w tau)^2);
    its p-p field effect is 2 A |H| alpha_L. Voltage mode adds 2 A |H| alpha_rho/(1+alpha_rho dT)."""
    sp = dict(DEFAULT_SPEC); sp.update(spec or {})
    tgt = sp["target_ppm"] * 1e-6
    wP = 2 * math.pi / sp["water_period_s"]
    H = 1.0 / math.sqrt(1.0 + (wP * tau) ** 2)
    k_res = ALPHA_CU / (1 + ALPHA_CU * (res["T_cu_mean_C"] - 20.0))
    return {
        "max_copper_dT_K_current_mode": tgt / ALPHA_L_CU,
        "max_copper_dT_K_voltage_mode": tgt / (ALPHA_L_CU + k_res),
        "max_water_amp_K_current_mode": tgt / (2 * H * ALPHA_L_CU),
        "max_water_amp_K_voltage_mode": tgt / (2 * H * (ALPHA_L_CU + k_res)),
        "max_psu_pp_ppm": sp["target_ppm"],
        "water_filter_gain": H,
    }
