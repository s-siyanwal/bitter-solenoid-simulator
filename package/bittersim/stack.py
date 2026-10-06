"""Explicit layer stacks: one annular conductor per row of a table.

The continuum path in ``design.py`` smears a Bitter winding into one coil.
Built magnets are tables of plates, half-plates, brass end pieces and
spacers, and their DC resistance has terms the smeared model does not carry.
This module keeps each layer and each extra term visible:

  E4 per layer   Bz(z) = (phi/2pi) mu0 C/2 sum_ends [asinh(zeta/R1) - asinh(zeta/R2)],
                 C = I / (t ln(R2/R1))                     (J = C/r, a slit plate)
  E3 per layer   uniform J = I / (t (R2-R1))               (wound or spiral turn)
  R per layer    J = C/r : R = rho phi / (t ln(R2/R1))
                 uniform : R = rho phi r_mid / (t (R2-R1))
  R_total        sum(R_layer) + n_joints r_joint + R_leads

Joint and lead terms default to 0 and are never fitted here. A caller adds
them only with a measured or separately bounded value.
"""
import math
from collections import namedtuple

import numpy as np

from . import fields
from .constants import RHO_CU_20

Layer = namedtuple("Layer", "r_in r_out z t phi rho profile sign")
Layer.__new__.__defaults__ = (2 * math.pi, RHO_CU_20, "bitter", 1.0)
Layer.__doc__ = ("Annular conductor: radii r_in, r_out, axial centre z and thickness t [m]; "
                 "current path angle phi [rad]; resistivity rho [ohm m]; profile 'bitter' "
                 "(J = C/r) or 'uniform'; sign = current direction (+1/-1).")


def mirrored(layers, sign=1.0):
    """Reflect a stack through z = 0; sign=-1 reverses the mirrored current (anti-Helmholtz)."""
    return list(layers) + [ly._replace(z=-ly.z, sign=ly.sign * sign) for ly in layers]


def axis_field(layers, z, current=1.0):
    """On-axis Bz [T] of the stack carrying ``current`` [A] in series."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    B = np.zeros_like(z)
    for ly in layers:
        frac = ly.phi / (2 * math.pi)
        if ly.profile == "bitter":
            C = current / (ly.t * math.log(ly.r_out / ly.r_in))
            b = fields.B_bitter_axis(ly.r_in, ly.r_out, ly.t, C, z - ly.z)
        elif ly.profile == "uniform":
            j = current / (ly.t * (ly.r_out - ly.r_in))
            b = fields.B_thick_axis(ly.r_in, ly.r_out, ly.t, j, z - ly.z)
        else:
            raise ValueError("unknown profile %r" % (ly.profile,))
        B += ly.sign * frac * b
    return B


def layer_resistance(ly):
    if ly.profile == "bitter":
        return ly.rho * ly.phi / (ly.t * math.log(ly.r_out / ly.r_in))
    if ly.profile == "uniform":
        return ly.rho * ly.phi * 0.5 * (ly.r_in + ly.r_out) / (ly.t * (ly.r_out - ly.r_in))
    raise ValueError("unknown profile %r" % (ly.profile,))


def stack_resistance(layers, n_joints=0, r_joint=0.0, r_leads=0.0):
    """DC resistance with each term reported separately [ohm]."""
    if r_joint < 0 or r_leads < 0 or n_joints < 0:
        raise ValueError("joint and lead terms must be non-negative")
    r_layers = float(sum(layer_resistance(ly) for ly in layers))
    r_j = n_joints * float(r_joint)
    return {"R_layers": r_layers, "R_joints": r_j, "R_leads": float(r_leads),
            "R_total": r_layers + r_j + float(r_leads)}


def implied_extra_resistance(R_measured, layers, n_joints, contact_area=None):
    """What a measured R leaves after the layers: total and per joint.

    Diagnostic only. If ``contact_area`` [m^2] is given, also returns the
    specific contact resistance [ohm m^2] the gap would need if joints alone
    carried it. Nothing is written back into the model.
    """
    gap = float(R_measured) - stack_resistance(layers)["R_layers"]
    out = {"R_gap": gap, "per_joint": gap / n_joints if n_joints else float("nan")}
    if contact_area and n_joints:
        out["specific_contact_ohm_m2"] = out["per_joint"] * float(contact_area)
    return out
