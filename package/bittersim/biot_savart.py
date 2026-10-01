"""Discretised Biot-Savart engine (PDF section "The Discretized Biot-Savart Law").

E6  B_seg(r) = mu0 I/(4 pi) (v x u)/|v x u|^2 (v.w/|v| - u.w/|u|),
    u = r1 - r, v = r2 - r, w = r2 - r1.

Note on sign: with u, v as defined, B for current flowing r1 -> r2 is
mu0 I/(4pi) (u x v)/|u x v|^2 (...) ; v x u = -(u x v).  The PDF prints
"v x u" in the equation but "u x v" in its code.  We use u x v, which is the
form that reproduces mu0 I/(2R) at the centre of a counter-clockwise loop
(+z), verified in tests.

The kernel is Numba-parallel over observation points (prange), and the
driver processes observation points in chunks (default 5000, as in the PDF)
so peak memory is O(chunk + N) instead of O(N*M).
"""
import math
import numpy as np
from numba import njit, prange

MU0 = 4.0e-7 * math.pi


@njit(parallel=True, cache=False)
def biot_savart_kernel(obs, s0, s1, cur):
    M = obs.shape[0]
    N = s0.shape[0]
    B = np.zeros((M, 3))
    for i in prange(M):
        rx = obs[i, 0]
        ry = obs[i, 1]
        rz = obs[i, 2]
        bx = 0.0
        by = 0.0
        bz = 0.0
        for j in range(N):
            ux = s0[j, 0] - rx
            uy = s0[j, 1] - ry
            uz = s0[j, 2] - rz
            vx = s1[j, 0] - rx
            vy = s1[j, 1] - ry
            vz = s1[j, 2] - rz
            wx = vx - ux
            wy = vy - uy
            wz = vz - uz
            cx = uy * vz - uz * vy
            cy = uz * vx - ux * vz
            cz = ux * vy - uy * vx
            c2 = cx * cx + cy * cy + cz * cz
            un = math.sqrt(ux * ux + uy * uy + uz * uz)
            vn = math.sqrt(vx * vx + vy * vy + vz * vz)
            if c2 < 1e-30 or un == 0.0 or vn == 0.0:
                continue  # observation point on the segment line
            f = MU0 * cur[j] / (4.0 * math.pi * c2) * ((vx * wx + vy * wy + vz * wz) / vn
                                                       - (ux * wx + uy * wy + uz * wz) / un)
            bx += f * cx
            by += f * cy
            bz += f * cz
        B[i, 0] = bx
        B[i, 1] = by
        B[i, 2] = bz
    return B


def biot_savart(obs, s0, s1, cur, chunk=5000):
    """Chunked generator driver: yields nothing large, returns (M,3) field."""
    obs = np.ascontiguousarray(np.atleast_2d(obs), dtype=np.float64)
    s0 = np.ascontiguousarray(s0, dtype=np.float64)
    s1 = np.ascontiguousarray(s1, dtype=np.float64)
    cur = np.ascontiguousarray(cur, dtype=np.float64)
    out = np.empty((obs.shape[0], 3))
    for k in range(0, obs.shape[0], chunk):
        out[k:k + chunk] = biot_savart_kernel(obs[k:k + chunk], s0, s1, cur)
    return out


def polygon_loop(R, z, n_ang, I):
    th = np.linspace(0.0, 2 * np.pi, n_ang + 1)
    p = np.column_stack([R * np.cos(th), R * np.sin(th), np.full_like(th, z)])
    return p[:-1], p[1:], np.full(n_ang, I)


def generate_coil_segments(R1, R2, L, amp, profile="bitter", n_turns=60, radial_steps=12,
                           angular_steps=72, helical=True):
    """Filament mesh of a Bitter (J=C/r) or uniform (J=j) coil.

    Corrected version of the PDF ``generate_bitter_geometry``: each of
    ``n_turns`` plates is one turn; each radial filament carries
    I_k = J(r_k) dr dz.  ``helical=True`` advances z by one plate pitch per
    turn (the PDF's helical transition); ``False`` gives planar rings, used
    for axisymmetric validation.
    """
    dz = L / n_turns
    dr = (R2 - R1) / radial_steps
    rc = R1 + (np.arange(radial_steps) + 0.5) * dr
    Jr = amp / rc if profile == "bitter" else amp * np.ones_like(rc)
    Ik = Jr * dr * dz
    a = np.arange(angular_steps)
    th1 = 2 * np.pi * a / angular_steps
    th2 = 2 * np.pi * (a + 1) / angular_steps
    S0, S1, CUR = [], [], []
    for p in range(n_turns):
        z0 = -L / 2.0 + (p + (0.0 if helical else 0.5)) * dz
        if helical:
            z1 = z0 + a / angular_steps * dz
            z2 = z0 + (a + 1.0) / angular_steps * dz
        else:
            z1 = z2 = np.full(angular_steps, z0)
        for k in range(radial_steps):
            r = rc[k]
            S0.append(np.column_stack([r * np.cos(th1), r * np.sin(th1), z1]))
            S1.append(np.column_stack([r * np.cos(th2), r * np.sin(th2), z2]))
            CUR.append(np.full(angular_steps, Ik[k]))
    return np.vstack(S0), np.vstack(S1), np.concatenate(CUR)


def richardson_order(Bh, Bh2, Bh4):
    """E7: p = ln[(B(h)-B(h/2))/(B(h/2)-B(h/4))]/ln 2 and extrapolated value."""
    p = math.log(abs((Bh - Bh2) / (Bh2 - Bh4))) / math.log(2.0)
    Bext = Bh4 + (Bh4 - Bh2) / (2.0 ** p - 1.0)
    return p, Bext
