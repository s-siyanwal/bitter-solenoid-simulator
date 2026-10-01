"""Lorentz forces and hoop stress (extension: the PDF mentions Lorentz forces
and hoop stresses qualitatively but gives no equations).

E20 force density f = J x B ; radial f_r = J_theta B_z
E21 thin-ring hoop stress (no radial load sharing, conservative):
    sigma_theta(r) = r J_theta(r) B_z(r)  ->  Bitter: sigma = C B_z(r) / lambda_fill
E22 axial (compressive) force on one half: F_z = sum_{z_i>0} 2 pi r_i I_i B_r(r_i, z_i)
"""
import numpy as np
from .fields import CoilLoops
from .elliptic import field_from_loops


def hoop_stress_profile(R1, R2, L, C, fill=1.0, n=41):
    coil = CoilLoops(R1, R2, L, C, "bitter", nr=48, nz=24, z_panels=4)
    # evaluate between radial quadrature nodes to avoid the filament singularities
    r = np.linspace(R1, R2, n)
    r = 0.5 * (r[1:] + r[:-1])
    _, Bz = coil.field(r, np.zeros_like(r))
    sigma = C * Bz / fill
    return r, Bz, sigma


def axial_force(R1, R2, L, C):
    coil = CoilLoops(R1, R2, L, C, "bitter", nr=16, nz=16, z_panels=4)
    Br = np.zeros_like(coil.a)
    for i in range(coil.a.size):
        mask = np.ones(coil.a.size, dtype=bool)
        mask[i] = False
        br, _ = field_from_loops(coil.a[i:i + 1], coil.z[i:i + 1], coil.a[mask], coil.z[mask], coil.I[mask])
        Br[i] = br[0]
    up = coil.z > 0
    # force on loop: dF = I dl x B ; azimuthal I with radial B gives axial force -I B_r 2 pi r
    Fz = -np.sum(2.0 * np.pi * coil.a[up] * coil.I[up] * Br[up])
    return float(Fz)
