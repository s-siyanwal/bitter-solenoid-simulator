# Test data provenance

| file | what | source |
|---|---|---|
| `claw_zs_layers.csv` | 106 Claw-ZS half-layers: Cu thickness and following spacer gap, inches | Parsed verbatim from `CZS_radia.py` (`unscaledZSdisc`, `unscaledgapsZS`) in github.com/olsenlab-science/Claw-ZS; file sha256 in the CSV header. The deposited 106th gap (20 in) is a plot sentinel and is left blank. |
| `claw_zs_fig4a_Bz_10A.csv` | 46 measured on-axis Bz points at 10 A | Vector-digitised from Fig. 4a of Hataway et al., arXiv:2607.02813v1 (bare-solver pass of 2026-10-06, `arxiv _papers/results_2026-10-06/B_axis_fig4a_measured_10A.csv`, sha256 `f691d188…bf6`). The Radia/COMSOL curves in that figure are not targets. |

Values quoted inline in `test_arxiv_checks.py` (no file):

- EPFL bulk-machined coil, Häusler et al., SciPost Phys. 6, 048 (2019), arXiv:1901.08791, Sec. 2 and 4.1: R1 = 32 mm, R2 = 72 mm, 22 mm tall, 31 turns, 1.3 mm pitch, 0.92 mm Cu width. Measured 1.28 G/A on axis 52.2 mm from the coil (Hall probe), R = 10.4(10) mΩ, L = 116(2) µH. Authors' Radia L = 94 µH, intrinsic R estimate 8.7 mΩ. The "832.2 G for a 104.4 mm pair" is a projection, not a measurement.
- Sabulsky et al., arXiv:1309.5330, Eq. 6 (thin Bitter arc on axis) and arc radii a1 = 31.75 mm, a2 = 50.80 mm.
- Kobelev, arXiv:1610.06607, Eq. 2.6 / 3.1 (long-stack field B = μ0 C ln(R2/r) in the winding).
