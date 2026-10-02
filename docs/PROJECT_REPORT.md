# Project report

Living record for `s-siyanwal/bitter-solenoid-simulator`.
The longer derivation is `report/bitter_solenoid_report.tex`.
`examples/update_report.py` rewrites `report/generated_numbers.tex` from `results/results.json` and `results/emulation_compare.json` at the end of `examples/run_all.py` and `examples/run_emulation.py`, and copies this file to `docs/PROJECT_REPORT.md`.
A later change of physics, of the page, or of a reported number appends a changelog entry in the same commit:

```bash
python examples/update_report.py --append "what changed"
```

The continuum solver is the reference. The emulation layer is separate. It is mesoscopic emulation, not molecular dynamics. It is not a particle model of the magnet.

## Published continuum numbers

Seed-1 differential evolution, 10605 evaluations. Re-runs on this tree reproduced the same printed optimum. Do not replace this table unless the algorithm itself changes; a changed run gets its own row.

| quantity | PDF initial guess | seed-1 optimum |
|---|---:|---:|
| R1 | 50.00 mm | 50.01 mm |
| R2 | 150.00 mm | 300.00 mm |
| L | 800.0 mm | 1102.0 mm |
| voltage | 19.537 V | 4.543 V |
| electrical power | 17.974 kW | 11.704 kW |
| hot-spot copper | 21.153 C | 25.946 C |
| homogeneity, 30 mm DSV | 154.15 ppm | 99.95 ppm |
| flow | 1834.9 L/min | 321.8 L/min |
| copper mass | 322 kg | 2575 kg |
| time to 85 C, lumped, cooling off | adiabatic, invalid once boiling | 4922 s |

The PDF guess violates 8 V and 100 ppm. The optimum is on the homogeneity, radius, and Reynolds bounds. Hoop stress on the optimum is 0.1176 MPa. The design is not stress-limited. Swiss-roll SNR gain 4.745 is a heuristic. mu_eff(0) = 1.

## Emulation, 50 draws, seed 1, PDF guess

Label: mesoscopic emulation, not molecular dynamics. Python is the reference. The browser subset uses another seed stream and does not resample ppm.

| quantity | nominal | p05 | p50 | p95 |
|---|---:|---:|---:|---:|
| B0 (T) | 0.500000 | 0.496394 | 0.500184 | 0.503441 |
| V (V) | 19.5366 | 19.3608 | 19.9514 | 20.4862 |
| P (W) | 17973.6 | 17811.8 | 18355.2 | 18847.2 |
| T_hot (C) | 21.1528 | 21.0021 | 21.2889 | 21.5206 |
| ppm | 154.146 | 152.976 | 154.021 | 154.999 |

Trip fractions on temperature, boiling onset, and 0.6 yield were 0. The median voltage sits above the continuum voltage because contact resistance is sampled only in the emulation. The adiabatic PDF-guess path reaches saturation near 479 s and is invalid after that. That is not the 4922 s optimum figure.

## Algorithms

Continuum evaluation locks `C` so the closed-form centre field is `B0`, iterates copper temperature and the mixed water rise, then places the hot spot on the inner hole row. Off-axis field is Gauss–Legendre loops with AGM elliptic integrals. Optimisation is differential evolution, then SLSQP, on `P_elec + P_pump`.

```
procedure OPTIMISE(seed, bounds, limits):
    population <- differential_evolution(objective, bounds, seed)
    best <- feasible member with smallest P_elec + P_pump
    polished <- SLSQP(best) with finite differences
    return polished if it improves the objective else best

procedure EMULATE(design, n, seed):
    nominal <- EVALUATE(design)
    for i in 1..n:
        q <- perturbed geometry, flow, inlet, resistivity lot
        Rc <- sum of log-normal contact draws
        q.B0 <- field of the nominal current in q
        sample <- EVALUATE(q)
        scale <- (R + Rc) / R * lot
        record B0, scaled V, P, T_hot, drift, mean free path, Johnson voltage
    return p05, p50, p95, trip fractions, model_grade
```

Drude identities, shared with the browser: `rho = m_e / (n e^2 tau)`, `v_d = J / (n e)`, `lambda = v_F tau`, Johnson `S_V = 4 kT R`. The Langevin step is only inside that representative volume.

## How to run

```bash
PYTHONPATH=package python -m pytest -q tests
PYTHONPATH=package python -m bittersim evaluate
PYTHONPATH=package python -m bittersim emulate --realizations 50 --seed 1
python examples/run_emulation.py
python examples/run_all.py
```

Open `docs/index.html` with no network. Presets: "PDF initial guess" and "seed-1 optimum".

## Changelog

### 2026-10-02 — emulation layer

Branch work that became `cursor/mesoscopic-emulation-5253`, then this Pages branch.
Added `catalog.py` and `emulation.py` without changing the continuum optimum.
Algorithms added: Drude RVE, Langevin cloud, Johnson and shot comparison, log-normal contact Monte Carlo, Antoine saturation, Bergles–Rohsenow ONB flag.
Numbers: PDF guess still prints P = 17.974 kW, V = 19.537 V, T_hot = 21.153 C. Optimum table unchanged. 50-draw seed-1 table recorded above.
Tests: 47 passed on the full tree after the missing modules were present (`python -m pytest -q tests`).

### 2026-10-02 — Pages demo and living report

Files: `docs/index.html`, `docs/bittersim.js`, `docs/PROJECT_REPORT.md`, `report/bitter_solenoid_report.tex`, `examples/update_report.py`, `PROJECT_REPORT.md`, tests.
What changed: one-page visual demo (scaled stack, Bz map, presets, quantile bars, Langevin cloud, pump-failure sketch with separate 4922 s and 479 s captions) and this report, which regenerates its number tables when results scripts run.
Algorithms: no change to the continuum optimum or to the Python Monte Carlo. Browser field map is a coarse loop grid for drawing only. ppm still comes from `evaluate`.
Numbers: unchanged from the tables above.
Tests: run in the same commit as this entry; see the changelog line appended by the test run if the count is repeated there.

Pages URL, once `main` is public and Pages is set to `/docs`: https://s-siyanwal.github.io/bitter-solenoid-simulator/

### 2026-10-02

Tests after the visual demo and the living report: 51 passed. Published continuum figures unchanged. PDF built locally from report/bitter_solenoid_report.tex (6 pages).
