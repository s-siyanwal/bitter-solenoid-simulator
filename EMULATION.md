# Simulation vs emulation

The continuum solver in `design.py` is the reference. This file describes the separate layer in `emulation.py` and `catalog.py`. It does not replace that solver and it does not change the published optimum.

Every emulation object is labeled **mesoscopic emulation, not molecular dynamics**. The magnet is about 2575 kg of copper, order 10^28 atoms. A 1 µm³ copper volume is still about 10^11 atoms. A laptop or a GitHub Pages tab cannot emulate the device atomistically. Two-temperature electron–phonon models matter on picosecond scales, not for this DC magnet.

## Model grades

| grade | meaning | where it is used |
|---|---|---|
| established | identity or a standard fluctuation formula | Drude `rho = m/(n e^2 tau)`, `v_d = J/(n e)`, Johnson `S_v = 4 kT R`, shot-vs-Johnson comparison, Antoine `Tsat` |
| approximate | engineering correlation or a stated scatter model | Kohler `Δρ/ρ ≈ a (ω_c τ)^2` with `a` of order 1, Langevin RVE picture, log-normal contact resistance, hole and plate tolerances, flow maldistribution stand-in, Bergles–Rohsenow ONB, housing convection, optional radial `ρ(T)` note |
| speculative | not a measurement on this coil | Hooge-like 1/f (default off), Swiss-roll SNR gain |

`model_grade` is written on the JSON object from `emulate`.

## Drude RVE

At the inner-radius copper, with catalog `n` and `ρ(T)`:

- `τ` from `ρ = m_e / (n e² τ)`
- Fermi speed 1.57×10⁶ m/s for copper (catalog value; Al 1350 uses its own)
- thermal speed `sqrt(8 kT / (π m))`
- drift `v_d = J / (n e)`
- mean free path `λ = v_F τ`
- transit time across one plate `d / v_d`, and across one mean free path `τ`
- Hall angle `ω_c τ`. At 0.5 T this is far below 1, so magnetoresistance is negligible here
- optional Kohler term, default on for the emulation only, off for the identity check

With Kohler off, `J·E` equals the continuum cell heat `ρ J²`. They agree well inside 1%.

The Langevin step is only inside that RVE, for a few thousand notional carriers, fixed seed:

`dv = -(v/τ) dt + (e/m) E dt + sqrt(2 kT / (m τ)) dW`

with `(e/m) E = v_d / τ`. The report is the ensemble mean against `v_d`, plus the standard error. At room temperature the thermal speed dwarfs the drift, so a few thousand carriers do not resolve `v_d`. That is the point of showing the error. It is not a sample of the magnet.

## Noise and scatter

Independent of the Langevin picture, a Monte Carlo (default 200 draws, seed required) perturbs the continuum design and compares each draw with `evaluate_design` on the same nominal inputs.

- Johnson–Nyquist `S_v = 4 k_B T R`, integrated from 1 Hz to 10 kHz by default (supply sense). A separate 1 Hz bin at the Larmor frequency is a spectral density only. It is not the MRI noise floor.
- Shot noise `2 e I` is compared with Johnson current noise and is negligible in this metal.
- Hooge-like 1/f is optional, default off, grade speculative.
- Plate-to-plate contact resistance is log-normal. Default median 1 µΩ per interface. The log-sigma (1) is an assumption, flagged as such. The sum is in series with the bulk resistance.
- Hole diameter and plate thickness ±1% (assumption), inlet sensor σ = 0.2 K, resistivity lot ±2%, inner-row flow ±10%. The continuum model has one velocity, so that velocity is the stand-in for the inner row.
- Radius and length use small normal tolerances (0.2 mm and 0.5 mm) so B0 can move when the supply current is held at the nominal value.
- Outputs: B0, V, P, T_hot, ppm, drift, mean free path, Johnson voltage. Bias is the median minus the deterministic value. Percentiles are 5 / 50 / 95.
- Trip probability: T_hot above 85 °C, ONB margin below 0, hoop stress above 0.6 times the selected temper yield.

`examples/run_emulation.py` writes `results/emulation_compare.json` and `figures/fig6_emulation.png`. It does not overwrite `results/results.json`.

## Thermal add-ons

Still continuum, not a CHF map.

- Saturation temperature from site pressure by the Antoine equation, intended for about 0.07–0.2 MPa.
- Onset of nucleate boiling for water uses Bergles–Rohsenow (1964). The correlation is named in the output. Other coolants do not get that water formula.
- Mixed outlet should stay below Tsat − 10 K.
- Velocity notes: below 0.5 m/s fouling, Re below 10⁴ correlation invalid, above about 8 m/s erosion caution for aerated water, above 15 m/s outside this catalog.
- The adiabatic pump-failure trajectory is flagged and stops being interpreted when the wall temperature reaches Tsat. The published time to 85 °C is not a boiling model.

## Catalog

Only ids in `catalog.py` are accepted. Unknown ids raise `ValueError` and list the allow-list.

Conductors: `ofhc_cu` (default), `cu_ag`, `cu_zr`, `al_1350`. The silver and zirconium entries are catalog approximations (higher strength, slightly higher resistivity), not heat certificates. Aluminium is allowed with a warning and is not recommended for this stack.

Coolants, liquid only: `di_water` (existing correlations; deionized, preferably deaerated), `water_glycol_30` (lower k, higher viscosity, electrically conductive), `galden_ht135` (dielectric, lower h).

Insulators: `polyimide_kapton` (default stand-in for 0.25 mm), `mica`, `ptfe`, `g10`. G-10 is a poor wet insulation choice and warns. Housing (`ss304`, `g10_bore`, `aluminum_6061`) is used for an external-convection estimate and a mass estimate. A conductive bore adds a shorted-turn warning, not a finite-element eddy solution.

## Audit warnings that are not “fixed”

The continuum optimum is unchanged. Evaluations carry a `warnings` list instead of silently editing old numbers.

1. Uniform `ρ(T)` with no radial feedback into `J = C/r`. An optional correction is reported only, default off.
2. Graded round holes are not a Florida-Bitter plate.
3. Helical slit and contact resistance are omitted from the continuum resistance. The emulation samples contact.
4. The thin-ring hoop stress is an upper estimate and is far below yield at 0.5 T. The design is not stress-limited.
5. No CHF map and no water chemistry. The pump-failure time is adiabatic and invalid at boiling.
6. Swiss-roll SNR is a heuristic. The emulation does not claim an SNR measurement.
7. Interactive output lists constraints that sit on their bounds.
8. There is no FEniCS model.

## Parallel particle emulation (`particles.py`)

A second, particle-level layer re-derives continuum results with particle estimators. It reports the difference from the continuum rather than replacing it:

- **P1:** Biot–Savart current elements give B0, the axis field, DSV homogeneity, NI, P, R and V.
- **P2:** Ornstein–Uhlenbeck conduction carriers give ρ(T) through the Green–Kubo / Einstein relation, with τ(T) from Bloch–Grüneisen. Θ_R = 343 K and RRR = 100 are assumptions.
- **P3:** Feynman–Kac heat walkers give the hot-cell conduction rise (E18).

All three use Numba `prange` with per-particle counter-based random streams, so the serial and parallel builds are bit-identical. `emulate(..., workers=n)` runs the tolerance Monte Carlo in a process pool and gives identical output for any `n`. `examples/run_particles.py` writes `results/particles.json` and figures 8 and 9. The README section "Particle emulation vs continuum simulation" holds the generated comparison and benchmark tables, and states which quantities must agree and which differ for physical reasons. Grade: established for P1 and P3 (exact estimators of the same equations). P2 is approximate, because the Bloch–Grüneisen τ(T) is a model; the carrier ensemble itself is an unbiased estimator.

## Non-goals

- Molecular dynamics, DFT, or a particle simulation of the whole magnet.
- An MRI noise-floor or SNR measurement.
- A certified magnet design or a manufacturing drawing.
- A critical-heat-flux map or a boiling two-phase solver.
- A finite-element shorted-turn or eddy-current solution.
- Changing `results/results.json` or the committed optimum.
