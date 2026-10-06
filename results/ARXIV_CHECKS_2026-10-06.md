# arXiv-folder checks and new comparison data (2026-10-06)

Branch `tests/arxiv-checks` (local, uncommitted). Bare solver only, no fitted parameters, no corrector.
Convention: relative residual = (measured − solver)/solver.

## Windows portability fixes
- `bittersim/validation.py` imported the POSIX-only `resource` module, so `tests/test_fields.py` and `tests/test_p0_literature.py` failed to collect on Windows. The import is now guarded and peak RSS reports NaN without it.
- 41 text-mode `open()` calls in `examples/`, `tests/` and `validation.py` had no encoding, so `make_docs.py --check` failed under cp1252. All now pass `encoding="utf-8"`.
- Suite on Windows (Anaconda Python 3.13, numpy 2.3.5, scipy 1.16.3, numba 0.62.1): 101 passed before the new tests were added.

## New tests (`tests/test_arxiv_checks.py`)
| test | kind | source | result |
|---|---|---|---|
| E4 thin-plate limit = Sabulsky Eq. 6 | code_check | arXiv:1309.5330 | 1e-8 at t = 10 µm; error ∝ t² |
| Sabulsky Eq. 6 = radial loop integral | code_check | arXiv:1309.5330 | rtol 1e-10 |
| E4 centre → μ0 C ln(R2/R1) as L → ∞ | code_check | Kobelev arXiv:1610.06607 Eq. 2.6 | O((R/L)²) |
| in-winding Bz used for hoop stress = μ0 C ln(R2/r), long stack | code_check | Kobelev Eq. 3.1 | 0.21 % of B0 (limit 0.5 %) |
| Claw-ZS layer table: 106 layers, 317.5 mm vs l = 320 mm | data check | CZS_radia.py | pass |
| Claw-ZS B0 peak at 10 A | measured | arXiv:2607.02813 Fig. 4a | −1.1 % (limit 2 %) |
| Claw-ZS profile, 46 points | measured | same | RMS 0.100 mT, max 0.204 mT |
| EPFL coil, B at 52.2 mm | measured | arXiv:1901.08791 Sec. 4.1 | 1.2835 vs 1.28 G/A: **−0.3 %** |
| EPFL coil, L vs Radia | cross-code | same | 91.6 vs 94 µH (−2.5 %) |

## New dataset: EPFL bulk-machined coil (arXiv:1901.08791, SciPost Phys. 6, 048)
A flat spiral (31 turns, 1.3 mm radial pitch, 32–72 mm, 22 mm tall), so the current density is uniform (E3), not Bitter 1/r. It exercises a different solver path from the Bitter campaigns.

| channel | measured | bare solver | rel. residual | note |
|---|---|---|---|---|
| B on axis, 52.2 mm, per A | 1.28 G/A (Hall) | 1.2835 G/A | −0.3 % | closes |
| R_dc | 10.4(10) mΩ | 8.41 mΩ (Cu, 20 °C, centre-line length) | +24 % (2.0 σ) | authors' own estimate is 8.7 mΩ; leads/joints are not modelled |
| L | 116(2) µH | 91.4 µH (nr → ∞) | +27 % | authors' Radia gives 94 µH; solver and Radia agree, both miss |
| water rise | 0.45 K/kW at 0.23 L/s | 1.04 K/kW (P/(ρ c Q)) | — | measured rise is less than half the energy-balance value; not a solver target until resolved |

## Pattern across campaigns (for model improvement, not for fitting)
- **Field closes on every campaign where the winding geometry is known:** EPFL −0.3 %, Claw-ZS −1.1 % at the peak, and Bates −2 to −3 %. The Olsen arc-model spread is a geometry question.
- **DC resistance is always under-predicted:** EPFL +24 %, Claw-ZS +42 %, and Olsen ×4.8 (a contact defect). All of these are copper-only models, so the common missing term is the joints, leads and contacts. The next physics step is an explicit lead/joint term with a measured or bounded value. Scaling the resistivity is not the fix.
- **L has a model-form split:** the smeared axisymmetric L is low for the flat spiral (EPFL) and high for the half-layer stacks (Claw-ZS, Olsen). AC vs DC (eddy/proximity currents) and lead inductance are the candidates. L should stay out of any corrector until that is modelled.

## Searched, not usable yet
- NHMFL 41.5 T all-resistive magnet: about 41.4 T at 48 kA (field factor ≈ 0.86 T/kA), first tested 21 Aug 2017. The coil geometry is only in Toth & Bole, IEEE TASC 28(3) 2018, doi:10.1109/TASC.2017.2775578 (paywalled). Blocked on dimensions.
- arXiv:2308.15476 is a magnetometer paper, not a magnet dataset.
- JQI Bitter Ioffe-Pritchard (arXiv:2008.11181, RSI 92, 033201): Hall maps at 150 A, plus R and L. Done: see the JQI section below.

## JQI Bitter Ioffe-Pritchard (added 2026-10-07)
Geometry from the authors' design notebook (`Cloverleaf_Trap_Reversed_v9.nb`): each stack is a 5.207 mm brass ring followed by 1 mm Cu layers; widths 12.7 mm; mid-radii 22.225/36.195 mm (curvature) and 61.722/75.692 mm (anti-bias); first layer at 19.05 mm from the centre. Bitter J = C/r, 1 A per layer.

| channel (pair, per A) | measured | bare solver | rel. residual |
|---|---|---|---|
| curvature B0 | 156.94(1) uT/A | 171.82 | −8.7 % |
| curvature B'' | 49.4(4) uT/(cm² A) | 54.05 | −8.6 % |
| anti-bias B0 (Helmholtz) | 160.46(4) uT/A | 176.87 | −9.3 % |
| anti-bias B' (anti-Helmholtz) | 24.2(1) uT/(cm A) | 27.42 | −11.7 % |
| curvature B''/B0 | 0.3148 /cm² | 0.3146 | −0.1 % |
| R curvature (Cu + brass, ρ_brass 6.6e-8 assumed) | 9.2(6) mΩ | 4.90 | +88 % |
| R anti-bias | 13.0(9) mΩ | 13.11 | −0.8 % |

Reading: the shape closes, which means the radii and axial placement are right, but the amplitude is uniformly about 9–12 % high. The geometric R reproduces the authors' own estimates (5.0 and 13.6 mΩ), so this is the layer table the authors used. The paper nonetheless claims Radia agrees with the measurement to better than 3 %, and the uniform-J (Radia-like) version of this table also comes out 10.7 % high on curvature B0. The as-built layer count or current calibration is therefore the open question. Nothing was tuned. The small curvature coil repeats the cross-campaign pattern of under-predicted R (leads, CuCr rods at 80 % IACS, notch contacts).

## Joint/lead resistance: does one term generalize? (2026-10-07)
New module `bittersim/stack.py` handles explicit layer tables. It computes E4/E3 per layer, uses the arc fraction φ/2π, mirrors pairs, and gives R per layer, with joint and lead terms that default to 0 and are never fitted. `examples/joint_diagnostic.py` (→ `results/joint_diagnostic.json`) divides each measured-R gap by the number of current-carrying interfaces:

| campaign | R measured | R layers | gap | interfaces | implied per interface |
|---|---|---|---|---|---|
| EPFL spiral | 10.4(10) mΩ | 8.41 mΩ | 1.99 mΩ | 0 (monolithic) | n/a |
| Claw-ZS | 5.3(2) mΩ | 3.66 mΩ | 1.64 mΩ | 210 | 7.8 ± 1.0 µΩ |
| JQI curvature | 9.2(6) mΩ | 4.90 mΩ | 4.30 mΩ | 20 | 215 ± 30 µΩ |
| JQI anti-bias | 13.0(9) mΩ | 13.11 mΩ | −0.11 mΩ | 22 | −5 ± 41 µΩ |

**Verdict: no single joint term generalizes.** The two JQI coils share their construction, yet they disagree by about 4σ. Claw-ZS sits 30× below the JQI curvature value, and EPFL has a 2 mΩ gap with no joints at all. The extra resistance is coil- and setup-specific (lead length, terminal blocks, where the voltage taps sit). It is not a material constant, so the solver keeps R_joint = R_leads = 0 by default.

Consequences:
- A measured R on a small coil (below about 10 mΩ) is not a clean model target until the voltage-tap positions are known.
- No corrector on R.
- Field remains the channel that transfers across campaigns.
