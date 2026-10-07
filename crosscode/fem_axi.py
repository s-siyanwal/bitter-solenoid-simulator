"""Axisymmetric FEM cross-checks with external open-source codes.

Builds a Gmsh geometry (r-z half plane) from a list of rectangular regions, then
runs either GetDP (ONELAB template Lib_Magnetodynamics2D_av_Cir.pro) or Elmer
(MagnetoDynamics2D) and returns field values at probe points and, for GetDP
harmonic runs, the region voltages of massive conductors.

Regions are dicts:
  {"r1","r2","z1","z2", "kind": "js",  "J": A/m^2}                uniform source
  {"r1","r2","z1","z2", "kind": "bitter", "C": A/m}              J = C / r source
  {"r1","r2","z1","z2", "kind": "iron", "mur": float}            linear magnetic
  {"r1","r2","z1","z2", "kind": "massive", "sigma": S/m, "I": A}  imposed current (AC)

Tool locations default to the portable installs in ../../../tools (override with
the env vars GMSH_EXE, GETDP_EXE, ELMER_BIN). Nothing here is fitted to data.
"""
import math
import os
import re
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.abspath(os.path.join(HERE, "..", "..", "..", "tools"))
GMSH = os.environ.get("GMSH_EXE", os.path.join(TOOLS, "onelab-Windows64", "gmsh.exe"))
GETDP = os.environ.get("GETDP_EXE", os.path.join(TOOLS, "onelab-Windows64", "getdp.exe"))
GETDP_TEMPLATE = os.path.join(os.path.dirname(GETDP), "templates", "Lib_Magnetodynamics2D_av_Cir.pro")
ELMER_BIN = os.environ.get("ELMER_BIN", os.path.join(TOOLS, "ElmerFEM-gui-nompi-Windows-AMD64", "bin"))

AIR, BND = 1, 2
REG0 = 1001               # first physical id of user regions


def _fwd(p):
    return p.replace("\\", "/")


# ------------------------------------------------------------------ geometry
def write_geo(path, regions, Rbox, Zbox, lc_reg, lc_air, lc_axis=None):
    """r-z half plane [0,Rbox]x[-Zbox,Zbox] with embedded rectangles (OpenCASCADE)."""
    L = ['SetFactory("OpenCASCADE");', 'Mesh.MshFileVersion = 2.2;', "Mesh.Algorithm = 6;"]
    for k, g in enumerate(regions, start=1):
        L.append("Rectangle(%d) = {%.9g, %.9g, 0, %.9g, %.9g};" % (k, g["r1"], g["z1"], g["r2"] - g["r1"], g["z2"] - g["z1"]))
    nb = len(regions) + 1
    L.append("Rectangle(%d) = {0, %.9g, 0, %.9g, %.9g};" % (nb, -Zbox, Rbox, 2 * Zbox))
    L.append("v() = BooleanFragments{ Surface{%d}; Delete; }{ Surface{1:%d}; Delete; };" % (nb, len(regions)))
    # after fragmenting, surfaces 1..n keep their tags (tools listed first); the rest is air
    L.append("air() = Surface{:};")
    L.append("air() -= {1:%d};" % len(regions))
    for k in range(1, len(regions) + 1):
        L.append("Physical Surface(%d) = {%d};" % (REG0 + k - 1, k))
    L.append("Physical Surface(%d) = {air()};" % AIR)
    L.append("bb() = Boundary{ Surface{:}; };")
    eps = 1e-5 * max(Rbox, Zbox)
    L.append("outer() = Curve In BoundingBox{%.9g,%.9g,-1, %.9g,%.9g,1};" % (-eps, -Zbox - eps, Rbox + eps, -Zbox + eps))
    L.append("outer() += Curve In BoundingBox{%.9g,%.9g,-1, %.9g,%.9g,1};" % (-eps, Zbox - eps, Rbox + eps, Zbox + eps))
    L.append("outer() += Curve In BoundingBox{%.9g,%.9g,-1, %.9g,%.9g,1};" % (Rbox - eps, -Zbox - eps, Rbox + eps, Zbox + eps))
    L.append("axis() = Curve In BoundingBox{%.9g,%.9g,-1, %.9g,%.9g,1};" % (-eps, -Zbox - eps, eps, Zbox + eps))
    L.append("Physical Curve(%d) = {outer()};" % BND)
    L.append("Physical Curve(%d) = {axis()};" % (BND + 1))
    L.append("MeshSize{ PointsOf{ Surface{air()}; } } = %.9g;" % lc_air)
    L.append("MeshSize{ PointsOf{ Surface{1:%d}; } } = %.9g;" % (len(regions), lc_reg))
    if lc_axis:
        L.append("MeshSize{ PointsOf{ Curve{axis()}; } } = %.9g;" % lc_axis)
    # refinement box around the device (axis to beyond the outer radius), graded to lc_air
    rmax = max(g["r2"] for g in regions)
    zlo, zhi = min(g["z1"] for g in regions), max(g["z2"] for g in regions)
    pad = 0.6 * max(rmax, zhi - zlo)
    L += ["Field[1] = Box;",
          "Field[1].VIn = %.9g; Field[1].VOut = %.9g;" % (lc_axis or lc_reg * 2, lc_air),
          "Field[1].XMin = 0; Field[1].XMax = %.9g;" % (rmax + pad),
          "Field[1].YMin = %.9g; Field[1].YMax = %.9g;" % (zlo - pad, zhi + pad),
          "Field[1].Thickness = %.9g;" % (2 * pad),
          "Field[2] = Constant; Field[2].VIn = %.9g; Field[2].SurfacesList = {1:%d};" % (lc_reg, len(regions)),
          "Field[3] = Min; Field[3].FieldsList = {1, 2};", "Background Field = 3;",
          "Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0;"]
    open(path, "w").write("\n".join(L) + "\n")


def mesh(geo, msh):
    r = subprocess.run([GMSH, geo, "-2", "-format", "msh22", "-o", msh], capture_output=True, text=True)
    if r.returncode or not os.path.isfile(msh):
        raise RuntimeError("gmsh failed:\n" + r.stdout[-2000:] + r.stderr[-2000:])
    n = re.findall(r"(\d+) nodes (\d+) elements", r.stdout)
    return n[-1] if n else None


# ------------------------------------------------------------------ GetDP
def write_pro(path, regions, probes, freq=None):
    G, F, CI = [], [], []
    names_js, names_iron, names_mass = [], [], []
    for k, g in enumerate(regions):
        nm = "R%d" % k
        G.append("  %s = Region[%d];" % (nm, REG0 + k))
        if g["kind"] == "js":
            names_js.append(nm)
            F.append("  js0[%s] = Vector[0, 0, %.12g];" % (nm, g["J"]))
        elif g["kind"] == "bitter":
            names_js.append(nm)
            F.append("  js0[%s] = Vector[0, 0, %.12g / X[]];" % (nm, g["C"]))
        elif g["kind"] == "iron":
            names_iron.append(nm)
            F.append("  nu[%s] = 1 / (%.12g * mu0);" % (nm, g["mur"]))
        elif g["kind"] == "massive":
            names_mass.append(nm)
            F.append("  sigma[%s] = %.12g;" % (nm, g["sigma"]))
            CI.append("      { Region %s; Value %.12g; }" % (nm, g["I"]))
    lin = names_js + names_mass
    pro = ["Flag_Axi = 1;", "Flag_FrequencyDomain = %d;" % (1 if freq else 0),
           "Freq = %.12g;" % (freq or 50), "resPath = \"res/\";",
           "Group {", "  Air = Region[%d]; Bnd = Region[%d];" % (AIR, BND)] + G + [
           "  Vol_Mag = Region[{Air%s}];" % "".join(", " + n for n in names_js + names_iron + names_mass),
           "  Vol_S0_Mag = Region[{%s}];" % ", ".join(names_js),
           "  Vol_C_Mag = Region[{%s}];" % ", ".join(names_mass),
           "}", "Function {", "  mu0 = 4e-7 * Pi;",
           "  nu[Region[{Air%s}]] = 1 / mu0;" % "".join(", " + n for n in lin),
           "  CoefGeos[] = 2 * Pi;"] + F + ["}",
           "Constraint {",
           # a_phi = 0 on the outer boundary AND on the symmetry axis (the axis condition is
           # not implied by the template; without it the solution grows like 1/r near r = 0)
           "  { Name MagneticVectorPotential_2D; Case { { Region Bnd; Value 0; } { Region Region[%d]; Value 0; } } }" % (BND + 1),
           "  { Name Current_2D; Case {"] + CI + ["  } }",
           "  { Name Voltage_2D; Case { } }", "}",
           'Include "%s";' % _fwd(GETDP_TEMPLATE)]
    form = "Magnetodynamics2D_av" if freq else "Magnetostatics2D_a"
    post = form
    # With the axis constraint the template's b is the physical field (verified against E3 to
    # 0.15 % on axis); out-of-plane +z of the (r, z) mesh plane is -phi, so the sign flips.
    ops = ['    Print[ b, OnPoint {%.12g, %.12g, 0}, Format SimpleTable, File > "res/b_pts.txt" ];' % (r, z)
           for r, z in probes]
    if freq:
        ops += ['    Print[ U, OnRegion %s, Format SimpleTable, File > "res/U.txt" ];' % n for n in names_mass]
    pro += ["PostOperation {", "  { Name probes; NameOfPostProcessing %s;" % post, "    Operation {",
            '    CreateDir["res/"];', '    DeleteFile["res/b_pts.txt"];'] + (['    DeleteFile["res/U.txt"];'] if freq else []) + ops + ["    }", "  }", "}"]
    open(path, "w").write("\n".join(pro) + "\n")
    return post


def run_getdp(workdir, regions, probes, Rbox, Zbox, lc_reg, lc_air, freq=None, lc_axis=None):
    os.makedirs(workdir, exist_ok=True)
    geo, msh, pro = (os.path.join(workdir, f) for f in ("case.geo", "case.msh", "case.pro"))
    write_geo(geo, regions, Rbox, Zbox, lc_reg, lc_air, lc_axis)
    nodes = mesh(geo, msh)
    post = write_pro(pro, regions, probes, freq)
    res = "Magnetodynamics2D_av" if freq else "Magnetostatics2D_a"
    r = subprocess.run([GETDP, "case.pro", "-msh", "case.msh", "-solve", res, "-pos", "probes", "-v", "2"],
                       cwd=workdir, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError("getdp failed:\n" + r.stdout[-3000:] + r.stderr[-3000:])
    out = {"nodes_elements": nodes, "b": _read_points(os.path.join(workdir, "res", "b_pts.txt"), complex_=bool(freq))}
    if freq:
        out["U"] = _read_region(os.path.join(workdir, "res", "U.txt"))
    return out


def _read_points(path, complex_=False):
    """Return (r, z, B_r, B_z) per probe; sign flipped to the (r, phi, z) convention."""
    rows = []
    for line in open(path):
        v = [float(x) for x in line.split()]
        if not v:
            continue
        if complex_:      # real/imag pairs per component
            br, bz = complex(v[3], v[4]), complex(v[5], v[6])
        else:
            br, bz = v[3], v[4]
        rows.append((v[0], v[1], -br, -bz))
    return rows


def _read_region(path):
    """OnRegion SimpleTable lines: '<region> <re> <im>'; '#' header lines are skipped."""
    vals = []
    for line in open(path):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        v = [float(x) for x in line.split()]
        vals.append(complex(v[-2], v[-1]) if len(v) >= 3 else complex(v[-1], 0))
    return vals


# ------------------------------------------------------------------ Elmer
def run_elmer(workdir, regions, probes, Rbox, Zbox, lc_reg, lc_air, lc_axis=None):
    """Axisymmetric magnetostatics with MagnetoDynamics2D (+ BSolver), sources and linear iron."""
    os.makedirs(workdir, exist_ok=True)
    geo, msh = os.path.join(workdir, "case.geo"), os.path.join(workdir, "case.msh")
    write_geo(geo, regions, Rbox, Zbox, lc_reg, lc_air, lc_axis)
    nodes = mesh(geo, msh)
    eg = subprocess.run([os.path.join(ELMER_BIN, "ElmerGrid.exe"), "14", "2", "case.msh", "-out", "mesh", "-autoclean"],
                        cwd=workdir, capture_output=True, text=True)
    if eg.returncode:
        raise RuntimeError("ElmerGrid failed:\n" + eg.stdout[-2000:])
    # ElmerGrid renumbers gmsh physical tags compactly in ascending order (bodies and
    # boundaries separately): air (1) -> 1, regions 1001.. -> 2.., outer boundary (2) -> 1.
    body_ids = sorted([AIR] + [REG0 + k for k in range(len(regions))])

    def body(pid):
        return body_ids.index(pid) + 1

    bnd_of = {BND: 1, BND + 1: 2}

    S = ['Header\n  CHECK KEYWORDS Warn\n  Mesh DB "." "mesh"\nEnd',
         "Simulation\n  Coordinate System = Axi Symmetric\n  Simulation Type = Steady\n"
         "  Steady State Max Iterations = 1\n  Output Intervals = 0\nEnd",
         "Constants\n  Permeability of Vacuum = %.12g\nEnd" % (4e-7 * math.pi)]
    mats = [("Relative Permeability = 1", None)]
    blocks = []
    bf = []
    S.append("Body 1\n  Target Bodies(1) = %d\n  Equation = 1\n  Material = 1\nEnd" % body(AIR))
    nb, nm, nf = 1, 1, 0
    for k, g in enumerate(regions):
        nb += 1
        mat = 1
        force = None
        if g["kind"] == "iron":
            nm += 1
            S.append("Material %d\n  Relative Permeability = %.12g\nEnd" % (nm, g["mur"]))
            mat = nm
        elif g["kind"] == "js":
            nf += 1
            force = nf
            S.append("Body Force %d\n  Current Density = %.12g\nEnd" % (nf, g["J"]))
        elif g["kind"] == "bitter":
            nf += 1
            force = nf
            S.append('Body Force %d\n  Current Density = Variable Coordinate 1\n    Real MATC "%.12g/tx"\nEnd' % (nf, g["C"]))
        S.append("Body %d\n  Target Bodies(1) = %d\n  Equation = 1\n  Material = %d\n%sEnd"
                 % (nb, body(REG0 + k), mat, ("  Body Force = %d\n" % force) if force else ""))
    S.append("Material 1\n  Relative Permeability = 1\nEnd")
    S.append("Equation 1\n  Active Solvers(2) = 1 2\nEnd")
    S.append('Solver 1\n  Equation = MgDyn2D\n  Procedure = "MagnetoDynamics2D" "MagnetoDynamics2D"\n'
             "  Variable = Potential\n  Linear System Solver = Direct\n  Linear System Direct Method = UMFPack\nEnd")
    S.append('Solver 2\n  Equation = MgDyn2DPost\n  Procedure = "MagnetoDynamics2D" "BSolver"\n'
             "  Discontinuous Galerkin = False\n  Linear System Solver = Iterative\n"
             "  Linear System Iterative Method = BiCGStab\n  Linear System Preconditioning = ILU0\n"
             "  Linear System Max Iterations = 2000\n  Linear System Convergence Tolerance = 1e-12\nEnd")
    # SaveLine reads the coordinates as pairs (one line per pair): each probe is a degenerate pair
    pts = "\n".join("    %.9g %.9g\n    %.9g %.9g" % (r, z, r, z) for r, z in probes)
    S.append('Solver 3\n  Exec Solver = After All\n  Equation = SaveProbes\n  Procedure = "SaveData" "SaveLine"\n'
             '  Filename = "probes.dat"\n  Polyline Coordinates(%d,2) = Real\n%s\n  Polyline Divisions(%d) = %s\nEnd'
             % (2 * len(probes), pts, len(probes), " ".join(["1"] * len(probes))))
    S.append("Boundary Condition 1\n  Target Boundaries(1) = %d\n  Potential = 0\nEnd" % bnd_of.get(BND, BND))
    open(os.path.join(workdir, "case.sif"), "w").write("\n\n".join(S) + "\n")
    open(os.path.join(workdir, "ELMERSOLVER_STARTINFO"), "w").write("case.sif\n1\n")
    r = subprocess.run([os.path.join(ELMER_BIN, "ElmerSolver.exe"), "case.sif"], cwd=workdir, capture_output=True, text=True)
    if r.returncode or "ERROR" in r.stdout:
        raise RuntimeError("ElmerSolver failed:\n" + r.stdout[-3000:])
    return {"nodes_elements": nodes, "b": read_elmer_probes(os.path.join(workdir, "probes.dat"))}


def read_elmer_probes(path):
    """SaveLine output -> [(r, z, B_r, B_z)] using the column names in <path>.names."""
    names = open(path + ".names").read()
    cols = {}
    for line in names.splitlines():
        line = line.strip()
        if ":" in line and line.split(":")[0].strip().isdigit():
            k, v = line.split(":", 1)
            cols[v.strip().lower()] = int(k) - 1
    rows = []
    for line in open(path):
        v = [float(x) for x in line.split()]
        if v:
            rows.append((v[cols["coordinate 1"]], v[cols["coordinate 2"]], v[cols["b 1"]], v[cols["b 2"]]))
    return rows[0::2]                    # two identical rows per degenerate probe line
