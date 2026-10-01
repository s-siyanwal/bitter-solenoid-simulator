"""bittersim - real-time simulator for a 0.5 T water-cooled copper Bitter solenoid.

Implements the "Computational Blueprint for a Real-Time 0.5 Tesla Water-Cooled
Copper Solenoid Simulator" using only NumPy, SciPy, Numba and Matplotlib APIs
that already existed before July 2018 (Python 3.6 syntax: no dataclasses,
no walrus operator).
"""
from .constants import MU_0
from .design import BitterDesign, evaluate_design
from . import fields, biot_savart, thermal, mechanics, swissroll, inductance

__version__ = "1.0.0"
__all__ = ["MU_0", "BitterDesign", "evaluate_design", "fields", "biot_savart",
           "thermal", "mechanics", "swissroll", "inductance"]
