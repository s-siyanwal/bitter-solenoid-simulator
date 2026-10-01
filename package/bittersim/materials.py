"""Temperature-dependent material properties.

Copper: linear resistivity rho(T) = rho20 (1 + alpha (T - 20)).
Water: standard engineering correlations (ASSUMPTION - the PDF only gives 20 C
constants). Valid roughly 0-100 C.  T in degrees Celsius.
"""
import numpy as np
from .constants import RHO_CU_20, ALPHA_CU, CP_WATER


def rho_cu(T):
    """Copper resistivity [Ohm m] at temperature T [C]."""
    return RHO_CU_20 * (1.0 + ALPHA_CU * (np.asarray(T, dtype=float) - 20.0))


def water_density(T):
    """Thiesen-type fit for liquid water density [kg/m^3]."""
    T = np.asarray(T, dtype=float)
    return 1000.0 * (1.0 - (T + 288.9414) / (508929.2 * (T + 68.12963)) * (T - 3.9863) ** 2)


def water_viscosity(T):
    """Dynamic viscosity [Pa s], Vogel equation."""
    T = np.asarray(T, dtype=float)
    return 2.414e-5 * 10.0 ** (247.8 / (T + 273.15 - 140.0))


def water_conductivity(T):
    """Thermal conductivity [W/m/K] (quadratic fit)."""
    T = np.asarray(T, dtype=float)
    return 0.5706 + 1.756e-3 * T - 6.46e-6 * T ** 2


def water_cp(T):
    """Specific heat [J/kg/K]; taken constant (PDF value)."""
    return CP_WATER + 0.0 * np.asarray(T, dtype=float)


def water_props(T):
    """Return dict of (rho, mu, k, cp, Pr, nu) at T [C]."""
    rho = float(water_density(T))
    mu = float(water_viscosity(T))
    k = float(water_conductivity(T))
    cp = float(water_cp(T))
    return {"rho": rho, "mu": mu, "k": k, "cp": cp, "Pr": cp * mu / k, "nu": mu / rho}
