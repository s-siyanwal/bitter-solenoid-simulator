# Design summary for independent physics review: 0.5 T water-cooled copper Bitter solenoid simulator

Equations are E1–E32 and results R1–R40. All R values come from `examples/run_all.py` ({{ENV}}). **[A]** marks an assumption not in the blueprint.

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

{{RESULTS_TABLE}}

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
