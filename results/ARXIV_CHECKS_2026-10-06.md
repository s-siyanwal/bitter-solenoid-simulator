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
- JQI Bitter Ioffe-Pritchard (arXiv:2008.11181, RSI 92, 033201): Hall maps at 150 A, plus R and L. The design repo is already in `experiment data/Bitter-Ioffe-Pritchard-master.zip`. Next candidate.
