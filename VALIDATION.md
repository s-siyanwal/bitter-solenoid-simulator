# VALIDATION

All numbers below are produced by `examples/run_all.py` -> `bittersim.validation` (environment: python 3.13.5, numpy 2.5.3, scipy 1.18.1, numba 0.68.0, matplotlib 3.11.2, cpus 8). The same checks run as pytest assertions in `tests/`.

## 1. Elliptic integrals and loop fields

| check | max relative difference |
|---|---|
| Numba AGM K(m), E(m) vs scipy.special.ellipk/ellipe, m in [0, 0.999] | 1.09e-14 |
| Numba loop field vs SciPy loop field (1600 off-axis points) | 2.04e-12 |
| On-axis loop (E5) vs code, axis branch | 1.93e-16 |
| On-axis loop (E5) vs elliptic branch at rho = 1e-7 m | 7.50e-13 |

Segment engine (E6) at the centre of an N-gon loop vs mu0 I/(2a):

| N | rel. error |
|---|---|
| 36 | 2.546e-03 |
| 72 | 6.351e-04 |
| 144 | 1.587e-04 |
| 288 | 3.967e-05 |

## 2. Long-solenoid limit (thin sheet, R = 50 mm, n = 1000 /m, I = 10 A)

| L/R | B loops [T] | mu0 n I [T] | rel. to mu0 n I | rel. to E2 |
|---|---|---|---|---|
| 10 | 0.012322340 | 0.012566371 | 1.942e-02 | 2.891e-12 |
| 100 | 0.012563858 | 0.012566371 | 1.999e-04 | 1.367e-14 |
| 1000 | 0.012566345 | 0.012566371 | 2.000e-06 | 1.463e-14 |

The deviation from mu0 n I falls as 2(R/L)^2, as expected from E2.

## 3. Thick-solenoid and Bitter closed forms (R1 = 50 mm, R2 = 150 mm, L = 400 mm)

| check | max relative difference |
|---|---|
| Gauss-loop model, uniform J, vs E3 (61 axial points) | 8.84e-14 |
| E3 at z=0 vs PDF centre formula | 0.00e+00 |
| Gauss-loop model, J = C/r, vs E4 (61 axial points) | 1.87e-13 |

## 4. PDF test 1: uniform thick solenoid with straight segments + Richardson extrapolation

40 planar turns x 8 radial filaments; B_z(0) exact (E3) = 2.240235636 T.

| angular steps | segments | rel. error vs E3 (all discretisation) | polygon error only |
|---|---|---|---|
| 45 | 14400 | 4.254e-04 | 3.214e-04 |
| 90 | 28800 | 1.844e-04 | 8.037e-05 |
| 180 | 57600 | 1.241e-04 | 2.009e-05 |
| 360 | 115200 | 1.090e-04 | 5.024e-06 |
| 720 | 230400 | 1.052e-04 | 1.256e-06 |

Richardson order from 180/360/720 steps: **p = 1.999976** (theory 2). Extrapolated B = 2.240468599 T equals the filament limit 2.240468599 T; the residual 1.04e-04 vs E3 is the radial/axial filament (midpoint) error of the 40x8 cross-section mesh, not the segment law. The PDF's < 0.195% criterion is met at every resolution tested.

## 5. PDF test 2: Bitter 1/r filaments vs closed form E4

80 plates x 24 radial x 360 angular (691200 segments). E4: 0.375520951 T.

| geometry | B_z(0) [T] | rel. error | transverse B [T] |
|---|---|---|---|
| planar rings | 0.375427782 | 2.481e-04 | 5.74e-19 |
| helical (PDF generator) | 0.375423845 | 2.586e-04 | 5.36e-04 |

The helical lead adds a small transverse field (the net axial current of the helix); it is reported but not part of the axisymmetric models.

## 6. Inductance and energy

| R [m] | length [m] | N | L numeric E9 [H] | L Nagaoka E10 [H] | rel. | Nagaoka / long-solenoid |
|---|---|---|---|---|---|---|
| 0.050 | 0.050 | 50 | 2.593288024e-04 | 2.593288024e-04 | 2.52e-11 | 0.525510 |
| 0.050 | 0.500 | 200 | 7.264760539e-04 | 7.264760539e-04 | 2.97e-12 | 0.920093 |
| 0.050 | 5.000 | 1000 | 1.957264411e-03 | 1.957264411e-03 | 2.48e-12 | 0.991562 |

Thick winding -> sheet limit (R1 = 50 mm, length 0.5 m, 200 turns):

| thickness [m] | rel. difference |
|---|---|
| 5e-03 | 3.341e-02 |
| 5e-04 | 3.568e-03 |
| 5e-05 | 3.593e-04 |

(The difference falls linearly with thickness, i.e. converges.)

Energy, long solenoid (R = 50 mm, 5 m, N = 1000, 10 A): L I^2/2 = 9.786322e-02 J, B^2/(2 mu0) x volume = 9.869604e-02 J, ratio 0.991562 (end-field deficit, equals the Nagaoka factor).

## 7. Homogeneity convergence (optimal design, 30 mm DSV)

| radial nodes | axial nodes | ppm | B0 [T] |
|---|---|---|---|
| 6 | 16 | 435.924137 | 0.499587537037 |
| 12 | 64 | 99.950895 | 0.499999998977 |
| 16 | 144 | 99.950247 | 0.499999999999 |
| 24 | 256 | 99.950247 | 0.500000000000 |

## 8. Thermal energy conservation

Sum of per-cell Joule heat over all hole rows vs analytic P (E31): relative residual -1.55e-16 (optimal design), 4.05e-16 (PDF initial guess). Mixed-outlet energy balance V0 I = m_dot cp dT holds by construction (dT computed from it); the test suite checks it to 1e-12 for random voltages/velocities.

## 9. Performance / memory (chunked Numba engine)

20000 segments x 20000 observation points in 0.746 s (5.36e+08 segment-point evaluations/s, 8 threads); process peak RSS 347 MB.

## Figures

![field](figures/fig1_field_homogeneity.png)

![convergence](figures/fig2_convergence.png)

![thermal](figures/fig3_thermal.png)

![swiss roll](figures/fig4_swissroll.png)

![mechanics](figures/fig5_mechanics.png)


## 10. BETA experimental benchmark (Bates et al., RSI 2018)

Anchor: Bates, Birmingham & Romero-Talamás, *Rev. Sci. Instrum.* **89**, 054704 (2018) — Tables I–II and measured ΔP on p.8. Produced by `bittersim.validation.v_beta_benchmark()`.

### Table I — electrical (BETA Eq 3 at T_avg = 44.5 °C)

Geometry: r1 = 20 mm, r2 = 69.86 mm, ℓ = 80.5 mm, t = 0.5 mm, λ = 0.8113, N = 77.625, I = 1175 A.

| quantity | Table I | continuum (this repo) | rel. difference |
|---|---|---|---|
| R | 17.679 mΩ | 17.705 mΩ | 1.5e-3 |
| V | 20.773 V | 20.803 V | 1.4e-3 |
| P | 24.41 kW | 24.44 kW | 1.4e-3 |
| L | 196.86 µH | 197.43 µH (E9 bitter) | 2.9e-3 |

### Table II — heat transfer

Elongated-hole Dh proxy ≈ 2.93 mm (2.5 mm radial × 4.5° span). Dittus–Boelter h at the published velocity band:

| v [m/s] | h (this repo) [W/m²K] | Table II band |
|---|---|---|
| 1.18 | ~6170 | 7191–8089 |
| 1.96 | ~9260 | 7191–8089 |

Wall temperatures in Table II (53.7 → 36.1 °C analytic) require the paper's multi-ring thermal-resistance network; the continuum hot-spot model here is a single-channel film+conduction estimate and is not re-solved for the BETA geometry in this check.

### Pressure drop

| source | ΔP |
|---|---|
| BETA measured | 7.72 kPa |
| BETA analytic (2 mm effective roughness) | 9.69 kPa |
| this repo, smooth Petukhov at mean-flow proxy | ~9.25 kPa |
| this repo, friction_multiplier=15 (MON stack rule) | much higher (rule is for smooth-tube baseline on large MRI stacks) |

The smooth-pipe ΔP at the proxy velocity lands close to BETA's analytic 9.69 kPa. The MON 10–20× multiplier is the continuum default recommendation for *smooth-tube* baselines on stacked MRI-scale channels; BETA already folds roughness into its analytic model.

## Not validated

* FEniCS curl-curl cross-validation (PDF): FEniCS 2018.1 is not pip-installable on current Colab/Kaggle or this environment; replaced by the independent closed-form and segment/loop cross-checks above.
* Swiss-roll SNR gain: heuristic model with no experimental calibration here.
