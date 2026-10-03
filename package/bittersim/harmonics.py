"""Spherical-harmonic field description, Z2/Z4 shim loop pairs, tolerance model, head preset.

H1  Bz(r, theta, phi) = sum_n sum_m r^n P_n^m(cos theta) [A_nm cos(m phi) + B_nm sin(m phi)]
    (Bz obeys Laplace's equation inside the bore, so this expansion is exact for n -> inf).
    Coefficients are reported in ppm of B0 at the DSV radius r0: A_nm r0^n / B0 * 1e6.
    P_n^m has no Condon-Shortley phase, so r P_1^1 cos(phi) = x (A11 is the x gradient).
H2  On the axis P_n(1) = 1, so the zonal terms A_n0 are the Taylor coefficients of Bz(0, z).
    `zonal_axis` fits them from the closed-form E4 axis field; `sh_fit` fits all (n, m) from
    Bz sampled on the DSV surface. For an axisymmetric coil both give the same zonal set.
H3  Shim pairs: two symmetric thin-loop pairs (radius a_s, at z = +-z_k) with currents solved
    so that the combined A_20 and A_40 vanish (linear 2x2 solve on the axis Taylor coefficients).
H4  Tolerance model: the coil is rigidly displaced (dx, dy, dz) and tilted about x by `tilt`;
    the field is evaluated in the coil frame and rotated back. This produces tesseral terms.
Assumptions are labelled ASSUMPTION; nothing here changes the published continuum numbers.
"""
import math
import numpy as np
from scipy.special import lpmv
from .constants import MU_0, RHO_CU_20
from .fields import CoilLoops, B_bitter_axis, B_loop_axis, sphere_points
from .elliptic import field_from_loops

N_AXIS = 12            # axis polynomial degree for H2
J_SHIM = 2.0e6         # ASSUMPTION: shim conductor current density [A/m^2] for the power estimate


def axis_nodes(r0, n=33):
    """Chebyshev nodes on [-r0, r0] (same nodes are used in docs/bittersim.js)."""
    k = np.arange(n)
    return r0 * np.cos(np.pi * (k + 0.5) / n)


def poly_fit_axis(z, f, r0, deg=N_AXIS):
    """Least-squares monomial fit in x = z/r0; returns c_n so that f = sum c_n (z/r0)^n."""
    V = np.vander(z / r0, deg + 1, increasing=True)
    c, *_ = np.linalg.lstsq(V, f, rcond=None)
    return c


def zonal_axis(R1, R2, L, C, r0, extra=None):
    """H2 zonal coefficients in ppm of Bz(0) at radius r0 (index n = 0..N_AXIS).
    extra: optional callable z -> additional on-axis Bz (e.g. shim loops)."""
    z = axis_nodes(r0)
    f = B_bitter_axis(R1, R2, L, C, z)
    if extra is not None:
        f = f + extra(z)
    c = poly_fit_axis(z, f, r0)
    return c / c[0] * 1e6, c[0]


def _real_sh_basis(pts, r0, nmax):
    r = np.linalg.norm(pts, axis=1)
    ct = np.where(r > 0, pts[:, 2] / np.where(r > 0, r, 1.0), 1.0)
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    cols, labels = [], []
    for n in range(nmax + 1):
        for m in range(n + 1):
            base = (r / r0) ** n * lpmv(m, n, ct) * (-1) ** m     # no Condon-Shortley phase
            cols.append(base * np.cos(m * phi)); labels.append(("A", n, m))
            if m > 0:
                cols.append(base * np.sin(m * phi)); labels.append(("B", n, m))
    return np.column_stack(cols), labels


def sh_fit(field_xyz, r0, nmax=8, n_pts=1200):
    """H1 least-squares fit of Bz over the DSV surface. Returns ({label: ppm}, B0, rms_resid_ppm)."""
    pts = sphere_points(r0, n_pts)[:-1]       # surface points (drop the centre)
    pts = np.vstack([pts, 0.5 * pts[::3]])     # interior shell pins the radial dependence
    Bz = field_xyz(pts)[:, 2]
    M, labels = _real_sh_basis(pts, r0, nmax)
    c, *_ = np.linalg.lstsq(M, Bz, rcond=None)
    B0 = c[0]
    resid = Bz - M @ c
    coef = {"%s%d%d" % lab: float(v / B0 * 1e6) for lab, v in zip(labels, c)}
    return coef, float(B0), float(np.sqrt(np.mean(resid ** 2)) / B0 * 1e6)


class LoopSet(object):
    """Generic set of coaxial filament loops (radius a, axial position z, current I)."""

    def __init__(self, a, z, I):
        self.a = np.ascontiguousarray(np.asarray(a, dtype=float))
        self.z = np.ascontiguousarray(np.asarray(z, dtype=float))
        self.I = np.ascontiguousarray(np.asarray(I, dtype=float))

    def __add__(self, other):
        return LoopSet(np.r_[self.a, other.a], np.r_[self.z, other.z], np.r_[self.I, other.I])

    def field_xyz(self, pts):
        pts = np.asarray(pts, dtype=float)
        rho = np.ascontiguousarray(np.hypot(pts[:, 0], pts[:, 1]))
        Br, Bz = field_from_loops(rho, np.ascontiguousarray(pts[:, 2]), self.a, self.z, self.I)
        with np.errstate(invalid="ignore", divide="ignore"):
            c = np.where(rho > 0, pts[:, 0] / np.where(rho > 0, rho, 1), 0.0)
            s = np.where(rho > 0, pts[:, 1] / np.where(rho > 0, rho, 1), 0.0)
        return np.column_stack([Br * c, Br * s, Bz])

    def axis(self, z):
        z = np.asarray(z, dtype=float)
        return sum(B_loop_axis(a, I, z - zz) for a, zz, I in zip(self.a, self.z, self.I))


def coil_loopset(R1, R2, L, C, nr=12, nz=32, z_panels=8):
    c = CoilLoops(R1, R2, L, C, "bitter", nr=nr, nz=nz, z_panels=z_panels)
    return LoopSet(c.a, c.z, c.I)


def shim_geometry(R1, gap=0.01, z_factors=(0.5, 1.5)):
    """ASSUMPTION: shim former radius a_s = R1 - gap (it uses `gap` of clear bore);
    pair k sits at z = +-z_factors[k] * a_s."""
    a_s = R1 - gap
    return a_s, [f * a_s for f in z_factors]


def pair_set(a_s, zk, I):
    return LoopSet([a_s, a_s], [zk, -zk], [I, I])


def solve_shims(R1, R2, L, C, r0, gap=0.01, z_factors=(0.5, 1.5)):
    """H3: currents (ampere-turns) of two loop pairs that null A_20 and A_40 at radius r0."""
    a_s, zs = shim_geometry(R1, gap, z_factors)
    z = axis_nodes(r0)
    main = poly_fit_axis(z, B_bitter_axis(R1, R2, L, C, z), r0)
    cols = [poly_fit_axis(z, pair_set(a_s, zk, 1.0).axis(z), r0) for zk in zs]
    A = np.array([[cols[0][2], cols[1][2]], [cols[0][4], cols[1][4]]])
    I = np.linalg.solve(A, -np.array([main[2], main[4]]))
    return a_s, zs, I


def shim_power_W(a_s, I_pairs, J=J_SHIM, rho=RHO_CU_20):
    """Resistive power of the shim pairs: P = rho J |NI| (2 pi a_s) per loop, two loops per pair."""
    return float(sum(2 * rho * J * abs(I) * 2 * math.pi * a_s for I in I_pairs))


def homogeneity_of(field_xyz, r0, n=400):
    pts = sphere_points(r0, n)
    B = np.linalg.norm(field_xyz(pts), axis=1)
    return float((B.max() - B.min()) / B[-1] * 1e6), float(B[-1])


def shimmed_homogeneity(R1, R2, L, C, dsv, gap=0.01, z_factors=(0.5, 1.5), nr=12, nz=32, z_panels=8):
    """Unshimmed and Z2/Z4-shimmed peak-to-peak |B| over the DSV, plus shim currents and power."""
    r0 = dsv / 2.0
    coil = coil_loopset(R1, R2, L, C, nr, nz, z_panels)
    a_s, zs, I = solve_shims(R1, R2, L, C, r0, gap, z_factors)
    shims = pair_set(a_s, zs[0], I[0]) + pair_set(a_s, zs[1], I[1])
    total = coil + shims
    ppm0, B0 = homogeneity_of(coil.field_xyz, r0)
    ppm1, B1 = homogeneity_of(total.field_xyz, r0)
    z0, _ = zonal_axis(R1, R2, L, C, r0)
    z1, _ = zonal_axis(R1, R2, L, C, r0, extra=shims.axis)
    return {"dsv_m": dsv, "ppm_unshimmed": ppm0, "ppm_shimmed": ppm1, "B0_unshimmed_T": B0,
            "B0_shimmed_T": B1, "shim_radius_m": a_s, "shim_z_m": zs, "shim_NI_A": [float(x) for x in I],
            "shim_power_W": shim_power_W(a_s, I), "zonal_ppm_unshimmed": [float(x) for x in z0[:9]],
            "zonal_ppm_shimmed": [float(x) for x in z1[:9]],
            "note": "Shim power assumes copper at J = %.1f A/mm^2 (ASSUMPTION)." % (J_SHIM / 1e6)}


def misaligned(field_xyz, dx=0.0, dy=0.0, dz=0.0, tilt=0.0):
    """H4: return a field function for the coil displaced by (dx,dy,dz) and tilted about x."""
    c, s = math.cos(tilt), math.sin(tilt)
    Rm = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    off = np.array([dx, dy, dz])

    def f(pts):
        local = (np.asarray(pts) - off) @ Rm          # R^T (p - off) for row vectors
        return field_xyz(local) @ Rm.T                # rotate B back to the lab frame
    return f


def tolerance_harmonics(R1, R2, L, C, dsv, dx=0.5e-3, tilt=1e-3, nmax=4):
    """Tesseral terms from an assumed assembly tolerance (ASSUMPTION: 0.5 mm offset, 1 mrad tilt)."""
    coil = coil_loopset(R1, R2, L, C)
    coef, B0, rms = sh_fit(misaligned(coil.field_xyz, dx=dx, tilt=tilt), dsv / 2.0, nmax=nmax)
    tess = {k: v for k, v in coef.items() if not k.endswith("0") or k[0] == "B"}
    zon = {k: v for k, v in coef.items() if k[0] == "A" and k.endswith("0")}
    return {"dx_m": dx, "tilt_rad": tilt, "zonal_ppm": zon, "tesseral_ppm": tess, "fit_rms_ppm": rms}


# ---------------------------------------------------------------- head-bore preset
HEAD = dict(R1=0.19, dsv=0.20, ppm_target=10.0)   # ASSUMPTION: 380 mm clear bore, 200 mm DSV,
                                                   # 10 ppm after Z2/Z4 shims (before active shimming)


def head_preset(R2_grid=None, L_grid=None, cooling=None, copper_usd_per_kg=None, usd_per_kWh=None):
    """Grid search (R2, L) at R1 = 0.19 m for the lowest electrical power whose Z2/Z4-shimmed
    homogeneity over a 200 mm DSV meets HEAD['ppm_target'] and whose hot spot stays <= 85 C.
    Cooling and plate parameters default to the seed-1 optimum. Costs are only reported when a
    price is passed in (no built-in prices)."""
    from .design import BitterDesign, evaluate_design
    R2_grid = np.linspace(0.30, 0.70, 9) if R2_grid is None else R2_grid
    L_grid = np.linspace(1.0, 3.0, 11) if L_grid is None else L_grid
    cool = dict(d_plate=5.99e-3, D_hole=5.915e-3, v_flow=1.697, pitch_factor=7.962)
    cool.update(cooling or {})
    best, rows = None, []
    for R2 in R2_grid:
        for L in L_grid:
            d = BitterDesign(R1=HEAD["R1"], R2=float(R2), L=float(L), dsv=HEAD["dsv"], **cool)
            res = evaluate_design(d, homogeneity=False)
            sh = shimmed_homogeneity(d.R1, d.R2, d.L, res["C_A_per_m"], HEAD["dsv"])
            ok = bool(sh["ppm_shimmed"] <= HEAD["ppm_target"] and res["T_hot_C"] <= 85.0)
            rows.append((float(R2), float(L), res["P_elec_W"], sh["ppm_unshimmed"], sh["ppm_shimmed"], ok))
            if ok and (best is None or res["P_elec_W"] < best[0]["P_elec_W"]):
                best = (res, sh)
    if best is None:
        return {"feasible": False, "grid": rows}
    res, sh = best
    out = {"feasible": True, "R1": res["R1"], "R2": res["R2"], "L": res["L"], "dsv_m": HEAD["dsv"],
           "P_elec_W": res["P_elec_W"], "P_pump_W": res["P_pump_W"], "V_total_V": res["V_total_V"],
           "I_A": res["I_A"], "T_hot_C": res["T_hot_C"], "mass_cu_kg": res["mass_cu_kg"],
           "flow_L_min": res["flow_L_min"], "ppm_unshimmed": sh["ppm_unshimmed"],
           "ppm_shimmed": sh["ppm_shimmed"], "shim_NI_A": sh["shim_NI_A"], "shim_power_W": sh["shim_power_W"],
           "violates_8V_supply": res["V_total_V"] > 8.0, "grid": rows}
    if copper_usd_per_kg is not None:
        out["copper_cost_usd"] = res["mass_cu_kg"] * copper_usd_per_kg
    if usd_per_kWh is not None:
        out["running_cost_usd_per_h"] = (res["P_elec_W"] + res["P_pump_W"] + sh["shim_power_W"]) / 1e3 * usd_per_kWh
    return out
