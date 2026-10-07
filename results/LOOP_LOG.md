# Hourly loop log

One entry per run of the hourly progress loop. Every run appends here and commits on `tests/arxiv-checks`, then pushes to GitHub, even when no code changed. The residual-ml repo is local-only (no remote), so its status is recorded here.

| run (IST) | work | simulator tests | pushed | resources |
|---|---|---|---|---|
| 2026-10-07 ~01:00 | JQI Bitter Ioffe-Pritchard vs Hall data: shape closes, amplitude 9–12 % high (notebook turn count) | 112 passed | `855c485` | CPU 8 %, RAM 19 GB free, C: 13.2 GB |
| 2026-10-07 ~02:00 | `bittersim/stack.py` layer stacks; joint diagnostic: one joint term does not generalize | 121 passed | `fffafba` | CPU 26 %, C: 13.1 GB |
| 2026-10-07 ~03:00 | DC stack inductance; JQI L vs Radia/measurement | 125 passed | `98b2818` | CPU 26 %, C: 13 GB |
| 2026-10-07 ~04:00 | PEEC AC impedance model | 129 passed | `bbb3c78` | CPU 16 %, C: 13 GB |
| 2026-10-07 ~05:00 | Olsen \|Z(f)\| held-out score: AC +22–25 % vs DC-L +63–81 % | 130 passed | `757eb90` | CPU 27 %, C: 13.2 GB |
| 2026-10-07 ~06:00 | Claw-ZS impedance with phase: L(f) within 3 %, nothing fitted | 133 passed | `c8f9e63` | CPU 15 %, C: 12.7 GB |
| 2026-10-07 ~07:00 | cross-campaign scorecard `results/GENERALIZATION.md` | 133 passed | `0df2b9f`, `a6847f2` | CPU 18 %, C: 12.7 GB |
| 2026-10-07 ~08:00 | residual-ml on Windows: 58/62/58/38 pass; bates-wt 69 pass only with numpy imported before torch (Anaconda OpenMP clash, not a repo bug); nothing to push | 133 passed | none (log started next run) | CPU 11 %, C: 12.0 GB (falling; not from this work) |
| 2026-10-07 08:54 | started this log; each run now commits and pushes an entry | 133 passed | this commit | C: 12.0 GB |
| 2026-10-07 ~09:30 (loop stopped; manual run) | 1972 papers read (RRE bound check); figures fig12-fig18 (`examples/make_arxiv_figures.py`); LaTeX report Section 8 "Validation against measured magnets" (`report/measured_validation.tex`), changelog appended | 133 passed | this commit | C: ~12 GB |
| 2026-10-08 ~00:45 | cross-code: GetDP 4.0 + Elmer 26.2 installed (portable); `crosscode/` driver; cases A (EPFL, both FEM within 0.1-0.7 % of E3), B (Claw-ZS profile, Elmer = bittersim to 0.02 %), C (NASA: GetDP = Elmer to 0.1 mT; earlier in-house FEM overstated steel boost ~2.8x, so the <5 % NASA closure is withdrawn) | 133 passed | this commit | C: ~11 GB |
| 2026-10-08 ~02:00 | cross-code case D: GetDP a-v harmonic, massive conductors (6-plate check: DC R = plate law, L within 0.6 % of PEEC). Claw-ZS L: GetDP/PEEC/meas within 2 %. Olsen \|Z\| +18-22 % residual reproduced by FEM, so it is not a PEEC artefact | 133 passed | this commit | CPU ~10 %, C: 6.6 GB |
