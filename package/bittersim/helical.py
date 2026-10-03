"""Helical Bitter current path: non-axisymmetric field errors from the plate slit / overlap.

Geometry (ASSUMPTIONS, all overridable):
  - n_turns flat annular plates (one turn each), pitch dz = L / n_turns, J = C / r across the plate,
    discretised into `radial_steps` filaments carrying I_k = C/r_k dr dz.
  - Each plate carries current azimuthally over 2 pi - overlap at constant z. In the overlap
    sector (angle `overlap_deg`) the current crosses to the next plate: modelled as a linear
    ramp in z over that sector (the slit / overlap transition).
  - The slit of plate p sits at phi_p = p * slit_advance_deg. advance = 0 puts every slit in one
    column ("aligned"); advance > 0 rotates the stack so the transitions spiral around.
  - The circuit is closed by a local return bus (ASSUMPTION): radial leads at both ends, an arc
    at radius R_bus in the top end plane, and a straight axial bus at R_bus.
  - A uniform-pitch helix (generate_coil_segments(helical=True)) is the third variant.
The field error is dB = B(path) - B(planar rings of the same mesh at the plate mid-planes), so mesh
discretisation error cancels in the comparison; dB is added to the converged loop model.
"""
import math
import numpy as np
from .biot_savart import biot_savart, generate_coil_segments

OVERLAP_DEG = 30.0           # ASSUMPTION: overlap / transition sector of each plate
BUS_GAP = 0.05               # ASSUMPTION: return bus at R2 + 50 mm


def _filaments(R1, R2, L, C, n_turns, radial_steps):
    dz = L / n_turns
    dr = (R2 - R1) / radial_steps
    rc = R1 + (np.arange(radial_steps) + 0.5) * dr
    return dz, rc, C / rc * dr * dz


def staircase_path(R1, R2, L, C, n_turns, overlap_deg=OVERLAP_DEG, slit_advance_deg=0.0,
                   radial_steps=12, seg_deg=2.5, ramp_segs=8):
    """Segments (s0, s1, I) of the slit/overlap staircase winding, plus end angles."""
    dz, rc, Ik = _filaments(R1, R2, L, C, n_turns, radial_steps)
    ov = math.radians(overlap_deg)
    adv = math.radians(slit_advance_deg)
    n_flat = max(4, int(round((2 * math.pi - ov) / math.radians(seg_deg))))
    phis, zs = [], []
    for p in range(n_turns):
        ph0 = p * adv
        z0 = -L / 2.0 + (p + 0.5) * dz
        f = ph0 + np.linspace(0.0, 2 * math.pi - ov, n_flat + 1)
        phis.append(f[:-1]); zs.append(np.full(n_flat, z0))
        span = ov + adv
        u = np.linspace(0.0, 1.0, ramp_segs + 1)[:-1]
        phis.append(ph0 + 2 * math.pi - ov + span * u); zs.append(z0 + dz * u)
    ph_end = (n_turns - 1) * adv + 2 * math.pi - ov            # last plate stops at its slit
    phis.append(np.array([ph_end])); zs.append(np.array([-L / 2.0 + (n_turns - 0.5) * dz]))
    ph = np.concatenate(phis); z = np.concatenate(zs)
    return _sweep(ph, z, rc, Ik), (0.0, -L / 2.0 + 0.5 * dz), (ph_end, -L / 2.0 + (n_turns - 0.5) * dz)


def _sweep(ph, z, rc, Ik):
    S0, S1, I = [], [], []
    for r, ik in zip(rc, Ik):
        p = np.column_stack([r * np.cos(ph), r * np.sin(ph), z])
        S0.append(p[:-1]); S1.append(p[1:]); I.append(np.full(len(p) - 1, ik))
    return np.vstack(S0), np.vstack(S1), np.concatenate(I)


def return_bus(rc, Ik, start, end, R_bus, arc_deg=2.5):
    """Close every filament: radial lead out at the end, arc at R_bus (top plane), axial bus down
    at the start angle, radial lead in at the start."""
    (ph_s, z_s), (ph_e, z_e) = start, end
    S0, S1, I = [], [], []
    d = (ph_s - ph_e + math.pi) % (2 * math.pi) - math.pi        # shortest arc end -> start
    n_arc = max(1, int(abs(math.degrees(d)) / arc_deg))
    arc = ph_e + d * np.linspace(0.0, 1.0, n_arc + 1)
    for r, ik in zip(rc, Ik):
        pts = [(r * math.cos(ph_e), r * math.sin(ph_e), z_e)]
        pts += [(R_bus * math.cos(a), R_bus * math.sin(a), z_e) for a in arc]
        pts += [(R_bus * math.cos(ph_s), R_bus * math.sin(ph_s), z_s), (r * math.cos(ph_s), r * math.sin(ph_s), z_s)]
        P = np.array(pts)
        S0.append(P[:-1]); S1.append(P[1:]); I.append(np.full(len(P) - 1, ik))
    return np.vstack(S0), np.vstack(S1), np.concatenate(I)


def planar_reference(R1, R2, L, C, n_turns, radial_steps=12, seg_deg=2.5):
    n_ang = int(round(360.0 / seg_deg))
    return generate_coil_segments(R1, R2, L, C, "bitter", n_turns, radial_steps, n_ang, helical=False)


def variant_segments(kind, R1, R2, L, C, n_turns, overlap_deg=OVERLAP_DEG, radial_steps=12, seg_deg=2.5,
                     bus=True, bus_gap=BUS_GAP):
    """kind: 'uniform' (uniform-pitch helix), 'aligned' (staircase, slits in one column),
    'rotating' (staircase, slit advances by the overlap angle each plate)."""
    dz, rc, Ik = _filaments(R1, R2, L, C, n_turns, radial_steps)
    if kind == "uniform":
        n_ang = int(round(360.0 / seg_deg))
        s0, s1, I = generate_coil_segments(R1, R2, L, C, "bitter", n_turns, radial_steps, n_ang, helical=True)
        start, end = (0.0, -L / 2.0), (0.0, L / 2.0)
    else:
        adv = 0.0 if kind == "aligned" else overlap_deg
        (s0, s1, I), start, end = staircase_path(R1, R2, L, C, n_turns, overlap_deg, adv, radial_steps, seg_deg)
    if not bus:
        return s0, s1, I
    b0, b1, bI = return_bus(rc, Ik, start, end, R2 + bus_gap)
    return np.vstack([s0, b0]), np.vstack([s1, b1]), np.concatenate([I, bI])


def delta_field(kind, R1, R2, L, C, n_turns, pts, **kw):
    """dB(pts) = s * B(helical path incl. bus) - B(planar reference mesh), with the current scale
    s chosen so that Bz at the isocentre is unchanged (the supply current is set to hit B0; a
    rotating stack carries 1 + advance/360 turns per plate). Returns (dB, s)."""
    s0, s1, I = variant_segments(kind, R1, R2, L, C, n_turns, **kw)
    p0, p1, pI = planar_reference(R1, R2, L, C, n_turns, kw.get("radial_steps", 12), kw.get("seg_deg", 2.5))
    allp = np.vstack([np.asarray(pts, float), [[0.0, 0.0, 0.0]]])
    Bh = biot_savart(allp, s0, s1, I)
    Bp = biot_savart(allp, p0, p1, pI)
    sc = Bp[-1, 2] / Bh[-1, 2]
    return sc * Bh[:-1] - Bp[:-1], float(sc)


# shim ladder: cumulative sets of ideally nulled Bz harmonics (labels as in harmonics.sh_fit)
LADDER = [
    ("Z2/Z4 shim pairs (currents solved on the ideal coil)", []),
    ("+ retune Z1, Z2, Z4", ["A10", "A20", "A40"]),
    ("+ X, Y", ["A11", "B11"]),
    ("+ ZX, ZY", ["A21", "B21"]),
    ("+ X2-Y2, XY", ["A22", "B22"]),
    ("+ all n = 3 tesseral", ["A31", "B31", "A32", "B32", "A33", "B33"]),
    ("+ all n = 4 tesseral", ["A41", "B41", "A42", "B42", "A43", "B43", "A44", "B44"]),
]
SHIM_NAMES = {"A10": "Z", "A11": "X", "B11": "Y", "A21": "ZX", "B21": "ZY", "A22": "X2-Y2", "B22": "XY",
              "A31": "Z2X", "B31": "Z2Y", "A32": "Z(X2-Y2)", "B32": "ZXY", "A33": "X3", "B33": "Y3",
              "A41": "Z3X", "B41": "Z3Y"}


def analyse(kind, R1, R2, L, C, n_turns, dsv, nmax=8, threshold_ppm=1.0, **kw):
    """Spherical harmonics of |B| (what the spins see; to first order equal to Bz + B_perp^2/(2 B0)),
    transverse field and the shim ladder for one path variant over one DSV. Ladder rows null the
    listed harmonics of |B| ideally (perfect shim coils of those orders), cumulatively.
    kind=None gives the ideal axisymmetric coil (loop model only)."""
    from .harmonics import coil_loopset, solve_shims, pair_set, _real_sh_basis
    from .fields import sphere_points
    r0 = dsv / 2.0
    coil = coil_loopset(R1, R2, L, C)
    a_s, zs, Ish = solve_shims(R1, R2, L, C, r0)
    shims = pair_set(a_s, zs[0], Ish[0]) + pair_set(a_s, zs[1], Ish[1])
    fit_pts = sphere_points(r0, 1200)[:-1]
    fit_pts = np.vstack([fit_pts, 0.5 * fit_pts[::3]])
    hom_pts = sphere_points(r0, 400)
    allp = np.vstack([fit_pts, hom_pts])
    if kind is None:
        dB, sc = np.zeros((len(allp), 3)), 1.0
    else:
        dB, sc = delta_field(kind, R1, R2, L, C, n_turns, allp, **kw)
    B_un = coil.field_xyz(allp) + dB
    B_sh = B_un + shims.field_xyz(allp)
    nf = len(fit_pts)

    def fit(B):
        M, labels = _real_sh_basis(fit_pts, r0, nmax)
        Bm = np.linalg.norm(B[:nf], axis=1)          # |B| is what the spins see
        c, *_ = np.linalg.lstsq(M, Bm, rcond=None)
        res = Bm - M @ c
        return c, labels, float(np.sqrt(np.mean(res ** 2)) / c[0] * 1e6)

    def pp(B):
        Bm = np.linalg.norm(B, axis=1)
        return float((Bm.max() - Bm.min()) / Bm[-1] * 1e6)

    c_un, labels, rms_un = fit(B_un)
    c_sh, _, rms_sh = fit(B_sh)
    lab = ["%s%d%d" % l for l in labels]
    ppm = lambda c: dict((k, float(v / c[0] * 1e6)) for k, v in zip(lab, c))
    coef_un, coef_sh = ppm(c_un), ppm(c_sh)
    Mh, _ = _real_sh_basis(hom_pts, r0, nmax)
    Bh = np.linalg.norm(B_sh[nf:], axis=1)
    ladder = [("unshimmed", pp(B_un[nf:]))]
    removed = []
    for name, keys in LADDER:
        removed += keys
        Bc = Bh.copy()
        for k in removed:
            i = lab.index(k)
            Bc -= c_sh[i] * Mh[:, i]                  # ideal shim of that order (to first order in |B|)
        ladder.append((name, float((Bc.max() - Bc.min()) / Bc[-1] * 1e6)))
    tess = dict((k, v) for k, v in coef_sh.items() if not (k[0] == "A" and k.endswith("0")) and k[0] in "AB" and k != "A00")
    needed = sorted([k for k, v in tess.items() if abs(v) >= threshold_ppm], key=lambda k: -abs(tess[k]))
    perp = np.hypot(B_un[nf:, 0], B_un[nf:, 1])
    dperp = np.hypot(dB[nf:, 0], dB[nf:, 1])
    return {"kind": kind or "ideal", "dsv_m": dsv, "current_scale": sc, "fit_rms_ppm": rms_sh,
            "coef_ppm_unshimmed": dict((k, v) for k, v in coef_un.items() if abs(v) >= 0.01),
            "coef_ppm_shimmed": dict((k, v) for k, v in coef_sh.items() if abs(v) >= 0.01),
            "tesseral_ppm": dict((k, v) for k, v in tess.items() if abs(v) >= 0.01),
            "B_perp_max_uT": float(perp.max() * 1e6), "B_perp_helical_max_uT": float(dperp.max() * 1e6),
            "B_perp_helical_range_uT": float((dperp.max() - dperp.min()) * 1e6),
            "ladder_ppm": ladder, "needed_tesseral": needed,
            "needed_shims": [SHIM_NAMES.get(k, k) for k in needed], "threshold_ppm": threshold_ppm}
