# bitter-solenoid-simulator

This is a real-time simulator for a **0.5 T water-cooled copper Bitter solenoid** with Swiss-roll RF metamaterials for early-stage fMRI. It implements *"Computational Blueprint for a Real-Time 0.5 Tesla Water-Cooled Copper Solenoid Simulator"*.

The constraints come from the blueprint:
- It runs on a 32 GB Ryzen 5 laptop, Kaggle and Colab.
- It uses only Python-3.6-compatible code and NumPy / SciPy / Numba / Matplotlib (plus ipywidgets) APIs that existed before July 2018. It uses no Magpylib, PyCharge or other later field libraries.

Every number below was produced by `examples/run_all.py` on python 3.13.5, numpy 2.5.3, scipy 1.18.1, numba 0.68.0, matplotlib 3.11.2. Nothing is hand-typed: `examples/make_docs.py` renders this file from `results/results.json`.

## Layout

```
package/bittersim/   simulation package
  elliptic.py        exact loop fields (scipy.special.ellipk/ellipe + Numba AGM kernel, prange-parallel)
  fields.py          ideal/thin/thick/Bitter closed forms, Gauss-loop winding model, DSV homogeneity
  biot_savart.py     blueprint straight-segment Biot-Savart engine (Numba, chunked), helical Bitter mesh, Richardson
  inductance.py      Maxwell mutual inductance, winding self-inductance, Nagaoka check
  design.py          coupled electrical / thermal / hydraulic evaluation of a design
  catalog.py         material allow-list (conductors, coolants, insulators, housings)
  emulation.py       mesoscopic emulation (Drude RVE, contact scatter, ONB flags); not molecular dynamics
  thermal.py         Re, Pr, Dittus-Boelter, Gnielinski, friction, pressure drop, hot spot, lumped transient
  mechanics.py       Lorentz force, hoop stress, axial compression
  swissroll.py       Pendry/Lorentzian mu_eff, skin-effect losses, heuristic SNR gain
  optimize.py        differential_evolution + SLSQP constrained optimisation
  validation.py      V&V suite used by VALIDATION.md
  cli.py             command line interface (python -m bittersim)
tests/               pytest: analytic limits, convergence, energy balance, web-demo parity, emulation
notebooks/           Colab/Jupyter notebook with ipywidgets sliders
docs/                static web demo (index.html + bittersim.js), no install needed
examples/            run_all.py (reproduces all results/figures), make_docs.py, run_emulation.py
figures/ results/    generated outputs
VALIDATION.md        validation tables and plots    DESIGN_SUMMARY.md  equation-labelled summary for review
EMULATION.md         mesoscopic emulation grades, formulas, and non-goals
```

## Quickstart

```bash
pip install -r requirements.txt          # or: pip install -r requirements-2018.txt on Python 3.6.6
pip install -e .                         # optional; tests/examples also work without installing
python -m pytest -q tests
python -m bittersim evaluate             # blueprint initial-guess design
python -m bittersim evaluate --R2 0.3 --L 1.102 --d-plate 0.00599 --D-hole 0.005915 --v-flow 1.697 --pitch-factor 7.962
python -m bittersim optimise --maxiter 100
python -m bittersim field --rho 0 0.01 0.02 --z 0 0 0.01
python -m bittersim swissroll
python -m bittersim catalog
python -m bittersim emulate --realizations 50 --seed 1
python examples/run_emulation.py
python examples/run_all.py && python examples/make_docs.py   # regenerate everything
```

(Without `pip install -e .`, prefix the commands with `PYTHONPATH=package`.)

### Open in Colab

The repository is **private**, so the usual "Open in Colab" badge link will not work for other people. To use it in Colab:
1. Upload `notebooks/bitter_solenoid_realtime.ipynb` (File → Upload notebook).
2. Make the repo available in the session, either by uploading the repo folder or zip, or with `!git clone https://<TOKEN>@github.com/s-siyanwal/bitter-solenoid-simulator.git`.
3. Run all cells.

Kaggle works the same way: add the repo as a dataset.

### Web demo (no install)

Open `docs/index.html` in any browser, including from a file URL. It is a plain-JS port of the analytic core: E4 field, exact elliptic loop fields for homogeneity, the thermal and hydraulic model, and Swiss-roll μ_eff. `tests/test_webdemo.py` checks that it agrees with the Python model to 1e-9, and that the PDF-guess and seed-1 presets match the printed precision. The page draws the plate stack, a coarse Bz map, constraint chips, and a button-triggered emulation (quantile bars, a Langevin cloud in one representative volume, a pump-failure sketch). Every emulation panel is labelled mesoscopic emulation, not molecular dynamics. An id that is not in the catalog is refused. The living write-up is [PROJECT_REPORT.md](PROJECT_REPORT.md) and [report/bitter_solenoid_report.tex](report/bitter_solenoid_report.tex). `examples/update_report.py` refreshes the numerical TeX tables after `run_all.py` or `run_emulation.py`.

### Simulation vs emulation

`evaluate` and `optimise` are the continuum simulator. `emulate` is a separate layer: a Drude representative volume, contact-resistance scatter, saturation and onset-of-nucleate-boiling flags, and percentile bands. It does not replace the continuum model, and it is not a particle model of the ~2575 kg magnet. Heuristics are labelled in the UI and in JSON as `model_grade`. The continuum optimum in the tables below is unchanged. See [EMULATION.md](EMULATION.md).

The owner authorized making this repository public so GitHub Pages can serve `docs/` from `main`. The site is public at `https://s-siyanwal.github.io/bitter-solenoid-simulator/`. Private Pages would need Enterprise Cloud. Details are in [docs/HOSTING.md](docs/HOSTING.md).

## Physics summary

| block | model |
|---|---|
| Field, closed forms | ideal μ0nI; finite thin solenoid (cos θ1 + cos θ2); thick uniform solenoid (log form); Bitter J = C/r on axis: Bz(0) = μ0C[asinh(L/2R1) − asinh(L/2R2)] |
| Field, off-axis | winding as Gauss–Legendre coaxial loops with exact K(m), E(m) loop fields (Numba, prange); the blueprint's straight-segment Biot–Savart engine is used as an independent check |
| Bitter electrics | Laplace → J = V0/(2πρr), I_plate = V0 d ln(R2/R1)/(2πρ), P = V0² L ln(R2/R1)/(2πρ); fill factor for insulation and holes; current fixed so Bz(0) = 0.5 T |
| Thermal | ρ(T) iteration; axial holes, Dittus–Boelter (blueprint) plus Gnielinski check; ṁc_p dT/dz = q′; Newton cooling; in-plate conduction to the hole; hot spot at R1 / outlet; lumped transient (normal operation and pump failure) |
| Hydraulics | Re, Petukhov/laminar friction, Δp with minor losses, pump power |
| Mechanics | J×B, thin-ring hoop stress, axial compressive force |
| Swiss roll | μ_eff = 1 − Fω²/(ω²−ω0²+jωΓ) with Pendry's geometric ω0 and Γ, skin-effect sheet resistance, heuristic SNR gain. **RF only:** μ_eff(0) = 1, so it has no effect on static B0 or its homogeneity in this model |
| Optimisation | differential_evolution then SLSQP; minimise P_elec + P_pump subject to 0.5 T, T_hot ≤ 85 °C, V ≤ 8 V, R1 ≥ 50 mm, ≤ 100 ppm, Δp ≤ 5 bar, Re ≥ 1e4 |

The full equation list (E1–E32) is in [DESIGN_SUMMARY.md](DESIGN_SUMMARY.md).

## Results: optimal design

The differential-evolution run (seed 1, 10605 evaluations) took 112 s. The SLSQP polish did not improve it.

| ID | quantity | value |
|---|---|---|
| R1 | Optimal inner radius R1 | 50.01 mm |
| R2 | Optimal outer radius R2 | 300.00 mm |
| R3 | Optimal stack length L | 1102.0 mm |
| R4 | Plate thickness / insulator | 5.990 mm / 0.25 mm (assumed) |
| R5 | Number of plates (turns) | 176.6 |
| R6 | Cooling holes: diameter, pitch/D, count | 5.915 mm, 7.962, 115 holes in 5 rows |
| R7 | Ampere-turns NI for Bz(0) = 0.5 T | 454948 A |
| R8 | Current I | 2576.1 A |
| R9 | Supply voltage / voltage per plate V0 | 4.543 V / 25.725 mV |
| R10 | Resistance | 1.7635 mOhm |
| R11 | Copper current density at R1 / R2 | 4.856 / 0.810 A/mm^2 |
| R12 | Electrical power (rho at mean Cu temperature) | 11.704 kW |
| R13 | Pump shaft power (eta = 0.7) | 81.0 W |
| R14 | Isocentre field from exact loop model | 0.499999999 T |
| R15 | Homogeneity, peak-to-peak over 30 mm DSV | 99.95 ppm |
| R16 | Water velocity / total flow | 1.697 m/s / 321.8 L/min (5.35 kg/s) |
| R17 | Reynolds / Prandtl number | 10068 / 6.901 |
| R18 | h Dittus-Boelter / Gnielinski | 8098 / 8118 W/m^2K |
| R19 | Pressure drop (Darcy f = 0.0314) | 10.57 kPa |
| R20 | Mixed water outlet temperature (inlet 20 C) | 20.522 C |
| R21 | Hottest channel: water rise / film dT / conduction dT | 2.370 / 2.900 / 0.675 K |
| R22 | Maximum (hot-spot) copper temperature, at R1, outlet | 25.946 C |
| R23 | Mean copper temperature | 20.901 C |
| R24 | Max hoop stress (E21, field profile) / simple bound C B0/lambda | 0.1176 / 0.1214 MPa |
| R25 | Axial compressive force on each half | 3656 N |
| R26 | Inductance / stored energy / L/R time constant | 1.1589 mH / 3846 J / 0.657 s |
| R27 | Lumped thermal time constant (63 %) / pump-failure time to 85 C | 29.8 s / 4922 s |
| R28 | Copper mass | 2575 kg |
| R29 | Larmor frequency at 0.5 T | 21.288739 MHz |
| R30 | Swiss roll tuned exactly to f_L: mu_eff(f_L), Q | 1.000 + 52.81j, Q = 96.8 |
| R31 | Best roll tuning f0/f_L, mu_eff(f_L), heuristic SNR gain | 1.01475, 17.367 + 5.701j, 4.745 |
| R32 | Swiss roll mu_eff at DC | 1.0 + 0.0j (no static-field effect) |
| R33 | PDF initial guess (R2 = 0.15 m, L = 0.8 m, v = 2.5 m/s, 2 mm plates) | P = 17.974 kW, V = 19.537 V (violates 8 V), 154.15 ppm (violates 100 ppm), T_hot = 21.15 C |
| R34 | Energy-balance residual (sum of cell heats vs E31) | -1.6e-16 |

At the optimum:
- **Active limits:** homogeneity (100 ppm), R2 = 0.3 m, R1 = 50 mm, Re ≥ 1e4, and plate thickness and hole pitch near their upper bounds.
- **Large margins:** hot-spot temperature and supply voltage.

At 0.5 T the design is therefore driven by field quality and size, not by cooling. The blueprint's initial guess (R33) violates both the 8 V limit and the homogeneity target.

## Validation (details in [VALIDATION.md](VALIDATION.md))

| ID | quantity | value |
|---|---|---|
| R35 | Richardson order of convergence, segment engine (theory 2) | 1.999976 |
| R36 | Segment engine vs E3 (thick uniform), 360 angular steps | 1.090e-04 relative |
| R37 | Bitter 1/r helical filaments vs E4 (691200 segments) | 2.586e-04 relative |
| R38 | Loop model vs E3 / E4 closed forms (61 axial points) | 8.8e-14 / 1.9e-13 |
| R39 | Long-solenoid limit L/R = 1000 vs mu0 n I | 2.00e-06 relative |
| R40 | Numerical sheet inductance vs Nagaoka (worst of 3) | 2.5e-11 relative |

## Figures

![field](figures/fig1_field_homogeneity.png)
![convergence](figures/fig2_convergence.png)
![thermal](figures/fig3_thermal.png)
![swiss roll](figures/fig4_swissroll.png)
![mechanics](figures/fig5_mechanics.png)

## Assumptions (not specified in the blueprint)

- **Geometry:**
  - Insulator thickness is 0.25 mm.
  - Cooling holes are circular, axial and run the full stack length L. They sit on uniform-pitch rows, are graded so the 1/r current profile is undisturbed, and the hole area reduces the copper cross-section uniformly.
  - The helical lead and slit are neglected in the axisymmetric field (the segment engine quantifies their effect: about 5e-4 T of transverse field at the centre of the test coil, see VALIDATION §5).
- **Field quality:** homogeneity limit 100 ppm peak-to-peak over a 30 mm DSV, with no shimming.
- **Materials:**
  - Copper: α = 3.93e-3 /K, k = 390 W/mK, density 8960 kg/m³, c_p = 385 J/kgK. ρ is uniform in each evaluation (no radial ρ(T) variation in the current distribution).
  - Water properties are temperature-dependent correlations evaluated at the mean bulk temperature. Inlet water is 20 °C (blueprint code).
- **Hydraulics:** minor-loss coefficient 1.5, pump efficiency 0.7, Δp ≤ 5 bar, Re ≥ 1e4 (Dittus–Boelter validity). Laminar flow uses Nu = 4.36.
- **Hot spot:** innermost hole row at the outlet. Water rise in that channel uses only that cell's heat (no mixing). Film ΔT is taken over the copper part of the wall. In-plate conduction is modelled as an annulus with an insulated outer boundary.
- **Transient:** single lumped copper node, quasi-steady water.
- **Optimiser:** the objective adds pump power to electrical power. Plate thickness, hole diameter and pitch are extra design variables; voltage is derived, not free. Bounds are listed in `optimize.py`.
- **Swiss roll:**
  - Geometry: r = 5 mm, lattice 12 mm, N = 28 turns, ε_r = 3, 10 µm Cu foil. The gap is solved for resonance. Loss is 50× the ideal-foil value, giving Q ≈ 97, typical of reported MRI rolls.
  - SNR model: a heuristic (reciprocity coupling plus magnetic-circuit flux guide plus μ″ noise term, E28). It is illustrative only.
  - Copper diamagnetism (χ ≈ −1e-5), which could cause ppm-level B0 distortion in a real build, is not modelled.

## Not implemented / deviations

- **Particle / DFT model of the winding:** not implemented. `emulation.py` is a mesoscopic emulation of a representative volume plus catalog and contact scatter.
- **FEniCS curl-curl FEA cross-validation:** FEniCS 2017.2/2018.1 cannot be pip-installed on current Colab/Kaggle or here. It is replaced by independent cross-checks: closed forms, the loop model and the segment engine.
- **"Complete 3D thermal map" and full conjugate heat transfer:** replaced by a 1-D axial channel model plus a local conduction correction for the hottest cell.
- **Analytic Jacobians/Hessians for SLSQP/trust-constr:** not used. The hole layout makes the objective piecewise, so the run uses DE plus SLSQP with finite differences (trust-constr is available but not needed).
- **Blueprint text vs. this code:**
  - The blueprint code's objective plugs the *total* 8 V into the per-plate power formula E31. Here V0 is the per-plate voltage, so P = V_total·I.
  - The blueprint states the turns-per-metre and power figures of its wire-wound example without derivation; they are not reproduced.
- **Versions:** numba 0.39.0 is listed in the blueprint as June 2018, but it was released 2018-07-11, so `requirements-2018.txt` pins 0.38.1.

## CI

`.github/workflows/tests.yml` runs pytest on Python 3.10 and 3.12 with current packages. A second job runs it in a `python:3.6.6` container with the pinned 2018 stack; that job is allowed to fail.

