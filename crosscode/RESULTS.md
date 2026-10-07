# Cross-code comparison: bittersim vs GetDP vs Elmer (2026-10-08)

Tools (portable installs, not in the repo): GetDP 4.0.0 + Gmsh 5.0.0 (ONELAB bundle, complex build), Elmer 26.2 (`ElmerFEM-gui-nompi`). `crosscode/fem_axi.py` generates the Gmsh geometry and the GetDP and Elmer input from a layer table, runs them, and parses the field. `crosscode/run_cases.py` holds the cases, with output in `crosscode/results/*.json`. Every FEM model uses the same layer table as the bittersim tests, and nothing is tuned to a measurement.

## Setup lessons (both now verified)
- **GetDP axisymmetric template** (`Lib_Magnetodynamics2D_av_Cir.pro`, `Flag_Axi = 1`): a = 0 must be imposed on the symmetry axis as well as on the outer boundary. Without it the solution behaves like C/r near the axis and the field is wrong by factors. With it, the template's `b` is the physical field; the sign flips because out-of-plane +z of the (r, z) mesh is −φ. A Grok CLI answer claiming no axis condition is needed was contradicted by the closed-form test and was not used.
- **Elmer**: ElmerGrid renumbers Gmsh physical tags compactly, so body and boundary indices are ranks. `SaveLine` reads coordinates as pairs, so each probe is a degenerate pair. `BSolver` works but is flagged obsolete.

## Case A: EPFL spiral, uniform J, verification against E3
Relative error versus the closed form E3 on the axis (z = 0, 26.1 mm, 52.2 mm):

| code | mesh (finest) | error |
|---|---|---|
| Elmer | 62 k nodes | −0.12 %, −0.25 %, −0.12 % |
| GetDP | 62 k nodes | −0.15 %, −0.60 %, −0.70 % |

Field per amp at 52.2 mm: measured 1.28 G/A; bittersim 1.2835; Elmer 1.2820; GetDP 1.2746. All three are within the measurement agreement.

## Case B: Claw-ZS axial profile, 106 half-layers smeared to ½ current, 10 A
| | peak [mT] | RMS vs measured (46 pts) |
|---|---|---|
| measured | 2.760 | – |
| bittersim (E4) | 2.786 | 0.0996 mT |
| Elmer (217 k nodes) | 2.786 | 0.0998 mT |

Elmer vs bittersim RMS: **0.0005 mT** (0.02 %).

## Case C: NASA MDF coil 1, 3 kW, linear-μ steel vessel
| rung | GetDP | Elmer | earlier in-house FEM | (meas − sol)/sol, verified |
|---|---|---|---|---|
| air core, t = 1.0 mm, N = 113 | 50.8 | 50.8 | 50.8 | +169.6 % |
| steel μr = 200, far shell | 63.8 | 63.8 | 85.5 | +114.5 % |
| steel μr = 2000, far shell | 66.2 | 66.2 | 93.6 | +106.9 % |
| steel μr = 2000, as-built shell | 65.9 | 65.9 | 100.1 | +107.7 % |
| thin disks t = 0.20, N = 140, as-built shell | 95.0 | 95.0 | 132.4 | +44.1 % |
| thin disks t = 0.20, N = 150, far shell | 100.9 | 100.9 | 137.2 | +35.7 % |

Convergence at μr = 2000: the result is 66.18 / 66.21 / 66.18 mT for coarse, fine (0.7 mm) and a 4 m box, identical in both codes.

**Finding.** The air core agrees exactly across all three codes, but the earlier in-house FEM overstated the steel's on-axis boost by about 2.8×. It had +42.8 mT at μr = 2000 against +15.4 mT in both established codes. The in-house solver was a tensor-product mesh with materials assigned by element centroid. Its "< 5 %" thin-disk closure therefore does not survive. With verified FEM the NASA field stays at +36 % to +44 % even under the thin-disk winding hypothesis, and +107 % with the published winding. The NASA discrepancy is open again.

## Still to do
- **Case D:** AC impedance L(f), R(f) of the Claw-ZS and Olsen stacks with massive conductors carrying an imposed series current (GetDP `Current_2D`), compared with the bittersim PEEC model and the measured data.
- **Optional:** a 3-D Elmer run with the Olsen 36° spacer sector.
