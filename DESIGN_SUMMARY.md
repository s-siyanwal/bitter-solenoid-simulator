# Design summary for independent physics review: 0.5 T water-cooled copper Bitter solenoid simulator

Equations are E1–E32 and results R1–R40. All R values come from `examples/run_all.py` (python 3.13.5, numpy 2.5.3, scipy 1.18.1, numba 0.68.0, matplotlib 3.11.2). **[A]** marks an assumption not in the blueprint.

## 1. Magnetostatics
- E1 Ideal solenoid: B = μ0nI.
- E2 Thin finite solenoid on axis: Bz = (μ0nI/2)(cosθ1 + cosθ2).
- E3 Uniform thick solenoid on axis: Bz(z) = (μ0j/2) Σ_{ζ=L/2∓z} ζ ln[(R2+√(R2²+ζ²))/(R1+√(R1²+ζ²))]. At z = 0 this is the blueprint formula.
- Bitter plate: Laplace's equation gives V = V0(2π−θ)/2π, so:
  - E29 J = V0/(2πρr)
  - E30 I_plate = V0 d ln(R2/R1)/(2πρ)
  - E31 P = V0² L ln(R2/R1)/(2πρ), with V0 the voltage per plate.
- E4 Bitter coil with J = C/r, on axis: Bz(z) = (μ0C/2) Σ_ζ [asinh(ζ/R1) − asinh(ζ/R2)], so Bz(0) = μ0C[asinh(L/2R1) − asinh(L/2R2)]. The current is chosen so that Bz(0) = 0.5 T exactly.
- E5 Loop on axis: μ0Ia²/[2(a²+z²)^{3/2}]. Off axis we use the exact elliptic-integral loop field, Bz = μ0I/(2π√((a+ρ)²+z²))[K + (a²−ρ²−z²)/((a−ρ)²+z²)E] (Bρ is analogous). Windings are represented as Gauss–Legendre loop sets (12×64). Numba computes K and E by AGM, matching scipy.special to 1e-14.
- E6 Segment Biot–Savart (blueprint engine): B = μ0I/4π (u×v)/|u×v|² (v·w/|v| − u·w/|u|). The blueprint equation prints v×u, which has the wrong sign; its own code uses u×v. The engine runs in Numba prange over 5000-point chunks.
- E7 Richardson order: p = ln[(B_h − B_{h/2})/(B_{h/2} − B_{h/4})]/ln2.
- Homogeneity: peak-to-peak |B|/|B0| over 400 points on a 30 mm DSV sphere **[A]**.

## 2. Inductance
- E8 Maxwell coaxial-loop mutual inductance: M = μ0√(ab)[(2/k−k)K − (2/k)E].
- E9 L = (N/NI)² ∬J(r)J(r′) 2∫(len−u)M du dr dr′.
- E10 Nagaoka sheet inductance: L = μ0πR²N²K_N/len, with K_N = 4/(3πk′)[(k′²/k²)(K−E) + E − k].

## 3. Electrical and thermal-hydraulic model
- Fill factor λ = d/(d+d_ins)·(1−f_hole) **[A: d_ins = 0.25 mm; holes graded so the 1/r profile is preserved]**. Then NI = CL ln(R2/R1) and P = 2πρC²L ln(R2/R1)/λ.
- Resistivity ρ(T) = 1.68e-8[1+0.00393(T−20)], iterated to self-consistency with the mean and hot-spot temperatures.
- Hole rows have a uniform pitch (pitch/D is a design variable). All holes carry the same velocity.
- E11 Re = ρvD/μ, Pr = c_pμ/k, h = Nu·k/D. Water properties are temperature-dependent **[A]**.
- E12 Dittus–Boelter: Nu = 0.023Re^0.8Pr^0.4 (used).
- E13 Gnielinski correlation (cross-check).
- E14 Darcy friction: 64/Re laminar, Petukhov turbulent.
- E15 Δp = (fL/D + 1.5)ρv²/2; pump power = ΔpQ/0.7 **[A]**.
- E16 Water heating: ṁc_p dT_f/dz = q′.
- E17 Wall temperature: T_s = T_f + q″/h, with q″ taken over the copper part of the hole wall.
- E18 Conduction in the hottest cell (annulus, uniform source ρJ²(R1)): ΔT = q_v/(4k)[2b²ln(b/a) − (b²−a²)].
- Hot spot = T_in + water rise (inner row, outlet) + film + conduction.
- E19 Lumped transient: C_th dT/dt = P20[1+α(T−20)] − (T−T_w)/R_th, solved for normal cooling and for pump failure.
- Energy check: V·I = ṁc_pΔT holds to 1e-12 in the tests.

## 4. Mechanics (the blueprint covers this only qualitatively)
- E20 Force density f = J×B.
- E21 Thin-ring hoop stress σθ = rJθBz = C·Bz(r)/λ (no load sharing, so conservative).
- E22 Axial force on one half: F_z = −Σ_{z>0} 2πrIB_r.

## 5. Swiss roll (acts at RF only)
- E23 Larmor frequency f_L = γB0/2π.
- E24 EMF = −jωμ0NSH.
- E25 μ_eff = 1 − Fω²/(ω² − ω0² + jωΓ).
- E26 Pendry (1999): F = πr²/a², ω0² = dc0²/(2π²ε_r r³(N−1)), Γ = 2σ_s/(μ0r(N−1)), with Γ multiplied by 50 **[A]**. Geometry r = 5 mm, a = 12 mm, N = 28, ε_r = 3, 10 µm foil **[A]**; the gap d is solved to hit ω0.
- E27 σ_s = ρ/(δ(1−e^{−t/δ})).
- E28 Heuristic SNR **[A; the blueprint gives no SNR model]**:
  - η0 = (a_c²/(a_c²+ℓ²))^{3/2}
  - η = 1/(1 + (1/η0−1)/μ′), set to 0 if μ′ ≤ 0
  - G = (η/η0)/√(1 + 0.05μ″)
- Static field: μ_eff(0) = 1, so the rolls are modelled as having **no effect on static B0 or its homogeneity** and are absent from every B0 calculation. Copper diamagnetism is not modelled.

## 6. Optimisation
- E32: minimise P_elec + P_pump over [R1, R2−R1, L, d, D_hole, v, pitch/D]. The current is fixed by the 0.5 T equality.
- Blueprint constraints: T_hot ≤ 85 °C, V ≤ 8 V, R1 ≥ 50 mm, R2 ≤ 0.3 m.
- Assumed constraints **[A]**: ≤ 100 ppm over the DSV, Δp ≤ 5 bar, Re ≥ 1e4.
- Method: differential_evolution (seed 1) with penalties, then an SLSQP polish. The polish did not improve the result, because the hole count is discrete.

## 7. Results

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
| R27 | Lumped thermal time constant (63 %) / pump-failure time to 85 C: lumped mean / hot-spot adiabatic bound | 54.3 s / 4922 s / 453 s |
| R28 | Copper mass | 2575 kg |
| R29 | Larmor frequency at 0.5 T | 21.288739 MHz |
| R30 | Swiss roll tuned exactly to f_L: mu_eff(f_L), Q | 1.000 + 52.81j, Q = 96.8 |
| R31 | Best roll tuning f0/f_L, mu_eff(f_L), legacy heuristic SNR gain (E28, superseded by R49) | 1.01475, 17.367 + 5.701j, 4.745 |
| R32 | Swiss roll mu_eff at DC | 1.0 + 0.0j (no static-field effect) |
| R33 | PDF initial guess (R2 = 0.15 m, L = 0.8 m, v = 2.5 m/s, 2 mm plates) | P = 17.974 kW, V = 19.537 V (violates 8 V), 154.15 ppm (violates 100 ppm), T_hot = 21.15 C |
| R34 | Energy-balance residual (sum of cell heats vs E31) | -1.6e-16 |
| R35 | Richardson order of convergence, segment engine (theory 2) | 1.999976 |
| R36 | Segment engine vs E3 (thick uniform), 360 angular steps | 1.090e-04 relative |
| R37 | Bitter 1/r helical filaments vs E4 (691200 segments) | 2.586e-04 relative |
| R38 | Loop model vs E3 / E4 closed forms (61 axial points) | 8.8e-14 / 1.9e-13 |
| R39 | Long-solenoid limit L/R = 1000 vs mu0 n I | 2.00e-06 relative |
| R40 | Numerical sheet inductance vs Nagaoka (worst of 3) | 2.5e-11 relative |
| R41 | Zonal Z2 / Z4 over 30 mm DSV (ppm at DSV radius, unshimmed) | -66.609 / -0.0585 |
| R42 | Peak-to-peak over 30 mm / 40 mm DSV: unshimmed -> Z2+Z4 shim pairs | 99.95 -> 0.264 ppm / 177.74 -> 1.364 ppm |
| R43 | Shim pairs: radius, z positions, NI, power (J = 2 A/mm^2 assumed) | 40.0 mm, +-20.0 / +-60.0 mm, 1.91 / 38.90 A, 0.689 W |
| R44 | Tesseral terms from 0.5 mm offset + 1 mrad tilt (assumed tolerance): A11 / B21 | -2.220 / 0.0888 ppm |
| R45 | Head preset (R1 = 190 mm, 200 mm DSV, <= 10 ppm after Z2/Z4; best on grid) | R2 700 mm, L 2600 mm, P 38.7 kW + pump 0.86 kW, 14.85 V, 2603 A, 31.5 t Cu, 1645 L/min, 1024 -> 9.92 ppm, T_hot 22.07 C (violates 8 V) |
| R46 | Field drift vs copper temperature at fixed current (numeric, isotropic expansion) | -16.500 ppm/K |
| R47 | B0 stability over 10 min, current-regulated PSU: ripple / drift / expansion / total | 2.000 / 0.333 / 2.527 / 4.439 ppm (94.5 Hz at f_L); voltage-regulated total 601.9 ppm |
| R48 | Max copper dT / water oscillation amplitude for 1 ppm (current / voltage mode) | 0.0606 / 2.54e-04 K ; 0.0459 / 1.93e-04 K |
| R49 | RF SNR, Swiss-roll slab vs same coil over air / vs coil on tissue (loss x50, best tuning f0/f_L) | 1.503 / 0.255 (f0/f_L = 1.0550, mu = 5.78 + 0.44j) |
| R50 | RF SNR vs air for loss multiplier 1 / 10 / 50 | 3.231 / 2.366 / 1.503 |
| R51 | RF resistances with slab: coil / tissue / slab | 0.0356 / 0.0359 / 0.1354 ohm |
| R52 | Particle emulation: speedup on 8 threads (field / carriers / walkers), serial == parallel bit-identical | 6.57x / 7.57x / 7.93x, yes |
| R53 | Particle vs continuum: B0, P (relative difference) | 2.05e-09 / 2.77e-07 |
| R54 | Particle vs continuum: DSV homogeneity | 99.73 ppm vs 99.95 ppm |
| R55 | Carrier (Bloch-Gruneisen) vs linear resistivity at mean Cu T | +0.085 % (+/- 0.141 % stat.) |
| R56 | Hot-spot copper T: particle / continuum | 25.946 / 25.946 C |

**At the optimum:**
- Homogeneity, R1 and R2 limits are active, and Re is at its floor.
- Temperature (R22) and voltage (R9) have large margins.

**Questions for reviewers:**
- Is E4 correct?
- Is the λ mapping valid for a perforated, insulated stack?
- Is the 1/r profile still valid with holes and with radial ρ(T)?
- Is the hot-spot model adequate?
- Is E26 correct? At exact tuning E25 gives μ′ = 1 and μ″ = FQ (R30), which absorbs rather than guides flux.
- Is E28 defensible?
