# bitter-solenoid-simulator

`bittersim` is a Python simulator for a **0.5 T water-cooled copper Bitter solenoid** with Swiss-roll RF metamaterials for early-stage fMRI. It models the magnetic field and its homogeneity, the coupled electrical / thermal / hydraulic behaviour of the perforated copper plate stack, the mechanics, and an RF layer, and it searches the design space for the lowest-power magnet that meets the field, cooling and supply limits. It implements *"Computational Blueprint for a Real-Time 0.5 Tesla Water-Cooled Copper Solenoid Simulator"*. One full design evaluation takes about 10 ms after Numba warm-up, which is what drives the interactive notebook and the browser demo.

**Published design** (seed-1 optimum, rows R1–R34 below): {{HEADLINE}}.

**[Live demo](https://s-siyanwal.github.io/bitter-solenoid-simulator/)** (no install) · [Equations E1–E32](DESIGN_SUMMARY.md) · [Validation](VALIDATION.md) · [Emulation notes](EMULATION.md) · [Project report](PROJECT_REPORT.md)

Highlights:
- **Physics:** closed-form solenoid fields, exact elliptic-integral loop fields for off-axis homogeneity, and the blueprint's straight-segment Biot–Savart engine as an independent check; the Laplace-derived J ∝ 1/r Bitter current distribution; a ρ(T)-coupled thermal–hydraulic model with a per-ring hot spot; Lorentz stresses, inductance and transients.
- **Design optimisation:** differential evolution plus an SLSQP polish over seven geometry and cooling variables.
- **fMRI layer:** spherical-harmonic field description, Z2/Z4 shim pairs, an assembly-tolerance model, a B0 stability budget and a quasi-static RF SNR model.
- **Cross-checks:** a validation suite against closed forms, and a particle emulation that re-derives the continuum results with independent estimators.
- **Reproducibility:** every result in this README is rendered from `results/*.json` by `examples/make_docs.py` (no result is hand-typed, and a test fails if the README is stale); `examples/reproduce.py` re-derives the published numbers and figures in about half a minute.

Design constraints taken from the blueprint:
- It runs on a 32 GB Ryzen 5 laptop, Kaggle and Colab.
- It uses only Python-3.6-compatible code and NumPy / SciPy / Numba / Matplotlib (plus ipywidgets) APIs that existed before July 2018. It uses no Magpylib, PyCharge or other later field libraries.

Contents: [Quickstart](#quickstart) · [How it works](#how-it-works) · [Results](#results-published-design) · [Reproducing the results](#reproducing-the-results) · [fMRI layer](#fmri-layer-shims-head-preset-b0-stability-rf-snr) · [Particle emulation](#particle-emulation-vs-continuum-simulation) · [Helical current path](#helical-current-path-non-axisymmetric-field-errors) · [Head preset](#head-preset-optimisation) · [Validation](#validation) · [Assumptions](#assumptions-not-specified-in-the-blueprint) · [Deviations](#not-implemented--deviations) · [Repository layout](#repository-layout) · [CI](#ci)

## Quickstart

```bash
git clone https://github.com/s-siyanwal/bitter-solenoid-simulator.git
cd bitter-solenoid-simulator
pip install -r requirements.txt      # NumPy, SciPy, Numba, Matplotlib, pytest
python -m pytest -q tests            # full test suite, about a minute
python examples/reproduce.py         # re-derive the published numbers and figures 1-5
```

`requirements-2018.txt` pins the historical pre-July-2018 stack instead (use it with Python 3.6.6).

Command line (`pip install -e .` installs the package and a `bittersim` command; without it, prefix these commands with `PYTHONPATH=package`):

```bash
python -m bittersim evaluate             # blueprint initial-guess design (R33)
python -m bittersim evaluate --R2 0.3 --L 1.102 --d-plate 0.00599 --D-hole 0.005915 --v-flow 1.697 --pitch-factor 7.962   # published optimum, rounded inputs
python -m bittersim optimise --maxiter 100   # current optimiser; see "Reproducing the results"
python -m bittersim field --rho 0 0.01 0.02 --z 0 0 0.01
python -m bittersim swissroll
python -m bittersim catalog
python -m bittersim emulate --realizations 50 --seed 1
```

### Notebook (Colab / Kaggle / Jupyter)

`notebooks/bitter_solenoid_realtime.ipynb` has ipywidgets sliders for the design explorer, the Swiss roll, a short optimiser run, the fMRI layer and the particle estimators.
- **Locally:** `pip install -r requirements-notebook.txt`, then open the notebook with Jupyter from the repository.
- **Colab:** upload the notebook (File → Upload notebook), make the repo available in the session (upload the repo folder or zip, or `!git clone https://github.com/s-siyanwal/bitter-solenoid-simulator.git`), and run all cells. Kaggle works the same way: add the repo as a dataset.

### Web demo (no install)

The [live demo](https://s-siyanwal.github.io/bitter-solenoid-simulator/) is served from `docs/` by GitHub Pages (`docs/index.html` is a short overview; `docs/demo.html` also works from a file URL). The demo shows the basic sliders and results first. Everything except the analytic core sits in a collapsed 'Advanced (work in progress)' section. The demo is a plain-JS port of the analytic core: E4 field, exact elliptic loop fields for homogeneity, the thermal and hydraulic model, and Swiss-roll μ_eff. `tests/test_webdemo.py` checks that it agrees with the Python model to 1e-9, and that the PDF-guess and seed-1 presets match the printed precision. The page draws the plate stack, a coarse Bz map, constraint chips, and a button-triggered emulation (quantile bars, a Langevin cloud in one representative volume, a pump-failure sketch). Every emulation panel is labelled mesoscopic emulation, not molecular dynamics. A reduced-N particle panel ports P1 (current elements: axis field and its error against E4, identical to Python at equal N; `tests/test_webdemo_particles.py`) and P2 (Green–Kubo carriers, statistically tested). An id that is not in the catalog is refused. Hosting details are in [docs/HOSTING.md](docs/HOSTING.md).

## How it works

`evaluate_design()` in `package/bittersim/design.py` turns a design (radii, length, plate thickness, cooling holes, water velocity) into every derived quantity:

1. **Current from the field target.** For J = C/r the on-axis field has a closed form (E4), so the current amplitude C is solved directly from Bz(0) = 0.5 T. Then NI = C L ln(R2/R1), and the plate count n = L / (d_plate + d_ins) is treated as continuous (I = NI / n).
2. **Electrics.** The copper fill factor is λ = d/(d + d_ins) · (1 − f_h), with the hole area fraction f_h from the hole layout. P = 2πρC²L ln(R2/R1)/λ; R and V follow, plus an optional contact resistance per plate interface.
3. **Thermal–hydraulic fixed point.** ρ(T) is iterated together with the mean copper temperature. The per-hole channel flow gives Re, Pr, h (Dittus–Boelter, with Gnielinski as a check), Petukhov friction and Δp. The mixed water rise and the mean film drop set the copper temperature.
4. **Hot spot.** For every ring of holes, at the outlet: the channel's water rise, plus the film ΔT, plus in-plate conduction from the hottest cell to the hole wall (E18), iterated on ρ(T_hot). The hottest ring is reported.
5. **Field quality.** The winding is replaced by 12 × 64 Gauss–Legendre coaxial loops with exact elliptic-integral fields (Numba, parallel). Homogeneity is the peak-to-peak |B| over 404 points on the DSV sphere.
6. **Post-processing.** Hoop stress and axial force, inductance, lumped transients (normal operation and pump failure), and the Swiss-roll μ_eff.

`optimise()` in `package/bittersim/optimize.py` runs `differential_evolution` over seven variables [R1, R2 − R1, L, d_plate, D_hole, v, pitch/D] with exterior penalties, then an SLSQP polish with explicit inequality constraints. The polished point is kept only if it is feasible and lowers P_elec + P_pump.

| block | model |
|---|---|
| Field, closed forms | ideal μ0nI; finite thin solenoid (cos θ1 + cos θ2); thick uniform solenoid (log form); Bitter J = C/r on axis: Bz(0) = μ0C[asinh(L/2R1) − asinh(L/2R2)] |
| Field, off-axis | winding as Gauss–Legendre coaxial loops with exact K(m), E(m) loop fields (Numba, prange); the blueprint's straight-segment Biot–Savart engine is used as an independent check |
| Bitter electrics | Laplace → J = V0/(2πρr), I_plate = V0 d ln(R2/R1)/(2πρ), P = V0² L ln(R2/R1)/(2πρ); fill factor for insulation and holes; current fixed so Bz(0) = 0.5 T |
| Thermal | ρ(T) iteration; axial holes, Dittus–Boelter (blueprint) plus Gnielinski check; ṁc_p dT/dz = q′; Newton cooling; in-plate conduction to the hole; hot spot at R1 / outlet; lumped transient (normal operation and pump failure) |
| Hydraulics | Re, Petukhov/laminar friction × optional stack multiplier (MON 10–20× on f only), Δp with minor losses, pump power |
| Contact R_c | Optional per-interface resistance in continuum R/V/P + overlap heat; headroom report |
| Hole grading | uniform (default) / montgomery / vinokur; per-ring hot spot |
| Mechanics | J×B, thin-ring hoop stress, axial compressive force |
| Swiss roll | μ_eff = 1 − Fω²/(ω²−ω0²+jωΓ) with Pendry's geometric ω0 and Γ, skin-effect sheet resistance, heuristic SNR gain. **RF only:** μ_eff(0) = 1, so it has no effect on static B0 or its homogeneity in this model |
| Optimisation | differential_evolution then SLSQP; minimise P_elec + P_pump subject to Bz(0) = 0.5 T (exact, through the current), T_hot ≤ 85 °C, V ≤ 8 V, ≤ 100 ppm over the DSV, Δp ≤ 5 bar, R2 ≤ 0.30 m and the variable bounds (R1 ≥ 50 mm). The Reynolds number is **not** a hard limit: below Re_min = {{RE_MIN}} a soft penalty {{RE_W}}·((Re_min − Re)/Re_min)² is added to the objective (in kW) to keep the heat-transfer correlations in their validity range |

The full equation list (E1–E32) is in [DESIGN_SUMMARY.md](DESIGN_SUMMARY.md).

**Simulation vs emulation.** `evaluate` and `optimise` are the continuum simulator. `emulate` is a separate layer: a Drude representative volume, contact-resistance scatter, saturation and onset-of-nucleate-boiling flags, and percentile bands. It does not replace the continuum model, and it is not a particle model of the ~2575 kg magnet. Heuristics are labelled in the UI and in JSON as `model_grade`. The continuum optimum in the tables below is unchanged. See [EMULATION.md](EMULATION.md). `particles.py` (Numba-parallel current elements, carriers and heat walkers) re-derives the continuum numbers from particles; see [Particle emulation](#particle-emulation-vs-continuum-simulation) below.

## Results: published design

The published optimum comes from a differential-evolution run (seed 1, {{NFEV}} evaluations) that took {{OPT_SECONDS}} s. The SLSQP polish did not improve it. The core design results were computed on {{ENV}}.

**Reynolds number.** {{RE_NOTE}}

{{RESULTS_TABLE}}

At the published optimum:
- **Active limits:** homogeneity (100 ppm), R2 = 0.3 m, R1 = 50 mm, and plate thickness and hole pitch near their upper bounds (for Re, see the note above).
- **Large margins:** hot-spot temperature and supply voltage.

At 0.5 T the design is therefore driven by field quality and size, not by cooling. The blueprint's initial guess (R33) violates both the 8 V limit and the homogeneity target.

![field](figures/fig1_field_homogeneity.png)
*Figure 1. Left: on-axis Bz of the Bitter coil (E4) compared with a uniform thick coil (E3), a thin coil (E2) and the ideal solenoid (E1) at the same NI. Right: |B| deviation in ppm near the centre (exact loop fields), with the 30 mm DSV.*

![thermal](figures/fig3_thermal.png)
*Figure 3. Left: water, wall and copper hot-spot temperature along the hottest channel. Right: lumped transient for normal cooling and for pump failure, with the 85 °C limit.*

![swiss roll](figures/fig4_swissroll.png)
*Figure 4. Left: Swiss-roll μ′ and μ″ near the Larmor frequency. Right: the legacy heuristic SNR gain (E28) against the roll tuning (superseded by the RF model, R49).*

![mechanics](figures/fig5_mechanics.png)
*Figure 5. Left: Bitter current density J ∝ 1/r in the copper. Right: midplane hoop stress and Bz(r, z = 0).*

## Reproducing the results

Every result in this README is rendered from `results/*.json` by `examples/make_docs.py`. The test suite runs `make_docs.py --check`, which fails if `README.md` or `DESIGN_SUMMARY.md` is out of date.

**Quick check** (no optimiser run; writes only to `repro_out/`, which is git-ignored):

```bash
python examples/reproduce.py         # exit code 1 if any quantity does not match
```

It re-evaluates the published design vector (`results/results.json` → `optimiser.x`) with the current code, recomputes the quantities behind rows R1–R51 with the same calls as the scripts below, compares them with the stored JSON values (relative tolerance 1e-6; quantities that are round-off-level errors or quantised to a grid get a matching absolute tolerance), and redraws figures 1–5. CI runs it on every push.

**Full regeneration.** These scripts overwrite the committed files under `results/` and `figures/`. `run_all.py` and `run_emulation.py` also rewrite `VALIDATION.md` and the report tables (`report/`, `docs/report/`), so compare with `git diff`. The runtimes are the ones recorded in the published runs, on {{CPUS}} CPU threads.

| README rows | script | writes | recorded runtime |
|---|---|---|---|
| R1–R40, figures 1–5 | `examples/run_all.py` | `results/results.json`, `results/validation.json`, `results/RESULTS.md`, `figures/fig1`–`fig5`, `VALIDATION.md`, report tables | optimiser {{OPT_SECONDS}} s; the whole script takes a few minutes |
| R41–R51, figure 7 | `examples/run_fmri.py` | `results/fmri.json`, `figures/fig7_fmri.png` | {{FMRI_SECONDS}} s |
| R52–R56, figures 8–9 | `examples/run_particles.py` | `results/particles.json`, `figures/fig8_particles.png`, `figures/fig9_particle_convergence.png` | not recorded |
| R57–R59, figure 10 | `examples/run_helical.py` | `results/helical.json`, `figures/fig10_helical.png` | {{HX_MINUTES}} min |
| R60–R61, figure 11 | `examples/run_head.py` | `results/head.json`, `figures/fig11_head.png` | {{HEAD_MIN_SECONDS}} min |
| mesoscopic emulation (around the blueprint initial guess) | `examples/run_emulation.py` | `results/emulation_compare.json`, `figures/fig6_emulation.png` (not committed), report tables | not recorded |
| all rendered tables | `examples/make_docs.py` | `README.md`, `DESIGN_SUMMARY.md` | seconds |

`run_fmri.py`, `run_particles.py` and `run_helical.py` read the published optimum from `results/results.json`; they do not rerun the optimiser. `run_particles.py`, `run_helical.py` and `run_head.py` accept `--quick` for a coarse smoke run whose numbers are not the published ones.

{{RE_REPRO}}

## fMRI layer: shims, head preset, B0 stability, RF SNR

`examples/run_fmri.py` writes `results/fmri.json` and `figures/fig7_fmri.png` from the published optimum (it does not rerun the optimiser). Modules: `harmonics.py` (spherical harmonics H1-H4, Z2/Z4 shim loop pairs, tolerance model, head preset), `stability.py` (S1-S4) and `rfsnr.py` (Q1-Q6). Every default input is an **assumption** and is copied into `fmri.json` under `assumptions`: shim J, head bore and DSV, PSU and chiller specs, the 1 ppm EPI target, coil, slab and tissue geometry, tissue conductivity, and the Swiss-roll loss multiplier.

{{FMRI_TABLE}}

How to read these numbers:
- **Shims.** Two thin-loop pairs inside the bore cancel Z2 and Z4. What is left is mostly Z6.
- **Head preset.** It sits at the upper end of the R2 grid in `head_preset()`, so it is not a converged optimum. The design is heavy and needs more than the 8 V supply. A copper-cost figure appears only when you pass a price: `head_preset(copper_usd_per_kg=...)`.
- **Stability.** At fixed current, B falls with copper temperature by the R46 coefficient (thermal expansion). A voltage-regulated supply adds the copper resistance coefficient on top, which gives the voltage-mode total in R47. So the magnet needs a current-regulated supply and chiller water held within the R48 amplitude.
- **RF.** This model is quasi-static: a surface loop, a laterally infinite uniaxial Swiss-roll slab, and a conducting half-space for the tissue. The roll array improves on the same coil held off over an air gap. It does not beat putting the coil straight on the tissue. The old E28 heuristic (R31) overstated the gain.
- **Helical path.** Tesseral terms from the helical current path are computed in the 'Helical current path' section below.

![fmri](figures/fig7_fmri.png)

## Particle emulation vs continuum simulation

`examples/run_particles.py` writes `results/particles.json`, `figures/fig8_particles.png` and `figures/fig9_particle_convergence.png` for the seed-1 optimum. The code is in `package/bittersim/particles.py`. It is **not molecular dynamics**: the magnet holds about 10^28 atoms. Instead, three particle estimators each solve the same physics as one continuum formula, so each pair can be compared:

- **P1, current elements.** {{P_FIELD_N}} samples × 32 Biot–Savart point elements (rings of 16 plus the z-mirror) give B on the axis and on the 30 mm DSV. The same samples give NI, P, R and V. Continuum counterparts: E4, the exact loop model, and E31.
- **P2, conduction carriers.** {{P_CARRIERS_N}} classical carriers follow an exact Ornstein–Uhlenbeck velocity process with relaxation time τ(T). The Einstein (Green–Kubo) relation σ = n e² D / (kT) turns this into ρ. τ(T) comes from a Bloch–Grüneisen lattice model (Θ_R = {{P_THETA}} K and RRR = {{P_RRR}}, both assumptions) calibrated to ρ(20 °C) = 1.68e-8 Ω m. The continuum uses a linear ρ(T).
- **P3, heat walkers.** {{P_WALKERS_N}} 2-D Bessel random walks run from the hot cell's edge to the cooling-hole wall. Feynman–Kac turns their mean exit time into the conduction rise. Continuum counterpart: E18.

**Parallelism.** All kernels use Numba `prange`. Each particle owns a counter-based random stream (splitmix64 keyed by seed and index), and reductions are summed serially in fixed chunks. As a result, serial and parallel runs are bit-identical, and the run script checks this with `np.array_equal`. `emulate(..., workers=n)` draws every random number first, then evaluates the tolerance realizations in a spawn process pool, so its output is identical for any worker count. Each realization takes {{P_MC_MS}} ms serially, so process start-up dominates at this size (measured speedup in the last row below). The particle kernels are where the parallel speedup lives.

{{PARTICLES_BENCH}}

{{PARTICLES_TABLE}}

{{PARTICLES_ROWS}}

**Where they must agree, and where they should not:**
- **Agree to estimator error** (same equations, different numerical method): B on the axis, NI, I, R, V, P at the same ρ, and the conduction rise at the same heat source.
  - B0 at the isocentre agrees exactly, by construction. The P1 sampling density equals the E4 integrand, so that estimator has zero variance. B0 is therefore an implementation check, not an independent test. The residual {{P_B0_ERR}} is the loop model's own quadrature error against E4. The independent tests are the off-centre axis points and the DSV.
  - The axis error is largest near the coil ends: {{P_AXIS_INNER}} max for |z| ≤ 0.4 L, and {{P_AXIS_ALL}} including the ends. That is where the importance density is worst.
  - V = P / I, so the shared ⟨1/r⟩ error cancels in V.
- **Homogeneity.** The ppm figure needs about 1e-7 relative field accuracy, so it converges slowly with N (fig9). The continuum reference is the loop model evaluated on the same 404 points ({{P_PPM_SAME}} ppm; R15 reports {{P_PPM_R15}} ppm).
- **Should differ (different physics):**
  - ρ(T). The Bloch–Grüneisen slope at 20 °C differs from the linear α = 0.00393 1/K. The power and hot-spot temperature differ accordingly. Between 21 and 26 °C the model gap is small ({{P_RHO_MODEL}} at the mean Cu temperature), smaller than the 1e6-carrier sampling error ({{P_RHO_SE}}). It grows to {{P_RHO_85}} at the 85 °C trip limit.
  - The emulated hot spot still takes the water rise and film drop from the continuum correlations, rescaled by the carrier ρ. Only conduction and resistivity come from particles.
- **Continuum contact resistance:** optional per-interface `R_c_ohm` (default 0) enters R, V, P and overlap-sector heat; V headroom is reported as max R_c/interface. Emulation still samples log-normal contacts.
- **Neither model includes:** a full Florida-Bitter elongated/staggered plate, or a turbulence-resolving CFD channel model. Classical carriers reproduce Drude σ under the relaxation-time approximation. Quantum (Fermi–Dirac) statistics change v_F and the mean free path, not σ.

![particles](figures/fig8_particles.png)
![particle convergence](figures/fig9_particle_convergence.png)

## Helical current path: non-axisymmetric field errors

`examples/run_helical.py` writes `results/helical.json` and `figures/fig10_helical.png`. It uses the straight-segment Biot–Savart engine (`biot_savart.py`) and the path builder in `helical.py` to model the real Bitter current path for the seed-1 optimum. All geometry inputs below are assumptions:

- {{HX_NT}} flat plates, each carrying one turn with J = C/r across the plate, discretised into {{HX_NR}} radial filaments.
- A {{HX_OV}}° overlap sector in which the current crosses to the next plate. In the model the crossing is a linear ramp in z over that sector.
- The circuit is closed by a return bus: radial leads at both ends and one axial bus at r = {{HX_BUS}} mm.

Three path variants are compared with the ideal axisymmetric coil:

| variant | description |
|---|---|
| uniform | uniform-pitch helix |
| aligned | flat plates with every slit in one angular column |
| rotating | flat plates with each slit advanced by the overlap angle, plate to plate |

The rotating stack carries more than one turn per plate. Its current is rescaled by {{HX_SCALE}} so that B0 is unchanged, and the same rescaling (≈1 for the other variants) is applied throughout.

**Method.**
- The helical error is dB = B(path) − B(planar rings on the same mesh), added to the converged loop model. Mesh error therefore cancels.
- The table reports the spherical harmonics of |B|, which is what the spins see. To first order |B| = Bz + B_perp²/(2 B0), so the transverse field enters too.
- Radial discretisation was checked at 48 vs 96 filaments (second-order convergence). The largest Richardson error estimate on any tesseral coefficient is {{HX_RICH}} ppm overall and {{HX_RICH_AL}} ppm for the aligned stack.
- The ladder below nulls whole harmonic orders cumulatively with ideal shims. The orders "needed" come from the shortest prefix of the ladder that brings the coil within 1 ppm of the ideal coil's Z2/Z4-shimmed value (an assumed budget). Within that prefix, only orders whose step lowers p-p by at least 0.1 ppm are listed. "Terms ≥ 1 ppm" lists the individual coefficients at or above 1 ppm.

{{HELICAL_TABLE}}

Peak-to-peak |B| [ppm] along the shim ladder:

{{LADDER_TABLE}}

**What this shows:**
- The Z2/Z4 pairs only remove zonal terms. The helical errors are mostly tesseral.
- With aligned slits, the slit column and its return bus act as an off-axis axial current. That produces a transverse field of the order of mT across the bore (R59) and X/Y gradient terms of |B|, so X and Y shims (plus higher orders, per the table) are needed.
- The uniform helix is close to the ideal coil.
- A stack whose slits rotate by the overlap angle each plate forms a short-pitch helix of transition currents. That produces large tesseral terms up to high order, which shims of order ≤ 4 cannot remove. A real design should avoid this stacking, or use a coaxial or bifilar return.
- Not modelled: plate-to-plate contact distribution, hole perforation of the current path, and lead geometry beyond the single bus.

![helical](figures/fig10_helical.png)

## Head preset: optimisation

`examples/run_head.py` writes `results/head.json` and `figures/fig11_head.png` (run time about {{HEAD_MIN_SECONDS}} min). It runs differential evolution with seeds 1–3 over R2, L, plate thickness, hole diameter, water velocity and pitch, at R1 = 190 mm and a 200 mm DSV.

- **Bounds (assumptions):** R2 up to 1.5 m, L up to 5 m, plates up to 20 mm, holes up to 8 mm.
- **Constraints:** ≤ 10 ppm after Z2/Z4 shims, T_hot ≤ 85 °C, Δp ≤ 5 bar, Re ≥ 1e4, with and without V ≤ 8 V. Unlike `optimize.py`, the head-preset optimiser (`harmonics.head_optimise`) still treats Re ≥ 1e4 as a hard constraint.
- **Objective:** electrical + pump + shim power. A third run minimises the supply voltage instead.

{{HEAD_TABLE}}

How to read it:
- **Power alone does not have an interior optimum here.** Power keeps falling slowly as R2 grows, so every seed ends near the widened R2 bound, with a very large copper mass. The hole diameter sits at its upper bound.
- **A practical head design needs a mass or cost term,** or a fixed R2. The table states which bounds bind.
- **The 8 V limit does not bind.** Thick plates mean fewer turns and a higher current at a low voltage.
- **The minimum supply voltage is set by the plate-thickness bound.**
- The old grid result (R45) is kept for reference.

{{NEXT_ROWS}}

![head](figures/fig11_head.png)

## Validation

Details, convergence tables and plots are in [VALIDATION.md](VALIDATION.md).

{{VALIDATION_TABLE}}

![convergence](figures/fig2_convergence.png)
*Figure 2. Relative error of Bz(0) for the straight-segment Biot–Savart engine on a thick uniform solenoid against the number of angular steps per turn, with a slope −2 reference.*

## Assumptions (not specified in the blueprint)

- **Geometry:**
  - Insulator thickness is 0.25 mm.
  - Cooling holes are axial and run the full stack length L. Default layout is **uniform-density** round holes (equal pitch, n_per_row ∝ r) — this is *not* graded. Optional `hole_layout_mode=montgomery|vinokur` implements MON (a1/r)² or BETA Vinokur grading; `elongated_aspect>1` makes area-preserving elliptical holes. Hole area reduces the copper cross-section via f_h.
  - The helical lead and slit are neglected in the axisymmetric field (the segment engine quantifies their effect: about 5e-4 T of transverse field at the centre of the test coil, see VALIDATION §5).
- **Field quality:** homogeneity limit 100 ppm peak-to-peak over a 30 mm DSV, with no shimming.
- **Materials:**
  - Copper: α = 3.93e-3 /K, k = 390 W/mK, density 8960 kg/m³, c_p = 385 J/kgK. ρ is uniform in each evaluation (no radial ρ(T) variation in the current distribution).
  - Water properties are temperature-dependent correlations evaluated at the mean bulk temperature. Inlet water is 20 °C (blueprint code).
- **Hydraulics:** minor-loss coefficient 1.5, pump efficiency 0.7, Δp ≤ 5 bar. Correlation validity enters the optimiser as a soft penalty below Re = {{RE_MIN}}, not as a hard limit (see [Results](#results-published-design) for the floor used by the published run). Laminar flow (Re < 2300) uses Nu = 4.36.
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

## Repository layout

```
package/bittersim/   simulation package
  elliptic.py        exact loop fields (scipy.special.ellipk/ellipe + Numba AGM kernel, prange-parallel)
  fields.py          ideal/thin/thick/Bitter closed forms, Gauss-loop winding model, DSV homogeneity
  biot_savart.py     blueprint straight-segment Biot-Savart engine (Numba, chunked), helical Bitter mesh, Richardson
  inductance.py      Maxwell mutual inductance, winding self-inductance, Nagaoka check
  design.py          coupled electrical / thermal / hydraulic evaluation of a design
  catalog.py         material allow-list (conductors, coolants, insulators, housings)
  emulation.py       mesoscopic emulation (Drude RVE, contact scatter, ONB flags; tolerance MC with workers=n); not molecular dynamics
  helical.py         helical Bitter current path (slit/overlap staircase, return bus) for the segment engine; harmonics and shim ladder
  particles.py       parallel particle emulation: Biot-Savart current elements, Green-Kubo carriers, Feynman-Kac heat walkers
  thermal.py         Re, Pr, Dittus-Boelter, Gnielinski, friction, pressure drop, hot spot, lumped transient
  mechanics.py       Lorentz force, hoop stress, axial compression
  swissroll.py       Pendry/Lorentzian mu_eff, skin-effect losses, heuristic SNR gain
  optimize.py        differential_evolution + SLSQP constrained optimisation
  harmonics.py       spherical-harmonic fit, Z2/Z4 shim loop pairs, tolerance (tesseral) model, head-bore preset
  stability.py       time-domain B0 stability: PSU ripple/drift, copper expansion, water-temperature drift
  rfsnr.py           quasi-static RF receive model: coil, Swiss-roll slab and tissue losses, SNR ratios
  validation.py      V&V suite used by VALIDATION.md
  cli.py             command line interface (python -m bittersim)
examples/            reproduce.py (quick check), run_all.py, run_fmri.py, run_particles.py, run_helical.py, run_head.py,
                     run_emulation.py (regenerate results/ and figures/), make_docs.py (renders README / DESIGN_SUMMARY),
                     update_report.py (report tables)
tests/               pytest: analytic limits, convergence, energy balance, web-demo parity, emulation, rendered docs
notebooks/           Colab/Jupyter notebook with ipywidgets sliders
docs/                GitHub Pages: index.html (overview), demo.html + bittersim.js (interactive demo, no install)
results/ figures/    generated outputs (JSON and PNG) behind every number and figure in this README
report/              LaTeX project report; PROJECT_REPORT.md is the living write-up
requirements*.txt    current stack, pinned 2018 stack, notebook extras
VALIDATION.md        validation tables and plots    DESIGN_SUMMARY.md  equation-labelled summary for review
EMULATION.md         mesoscopic emulation grades, formulas, and non-goals
```

## CI

`.github/workflows/tests.yml` runs pytest on Python 3.10 and 3.12 with current packages, followed by `examples/reproduce.py`. A second job runs pytest in a `python:3.6.6` container with the pinned 2018 stack (`requirements-2018.txt`); it must pass too. `.github/workflows/pages.yml` deploys `docs/` to GitHub Pages on pushes to `main`.
