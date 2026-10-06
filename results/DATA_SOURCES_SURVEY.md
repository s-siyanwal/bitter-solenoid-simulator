# Survey of public comparison data (2026-10-07)

Searched: IEEE TASC/TMAG (abstracts plus open author copies on public.magnet.fsu.edu, OSTI and EMFL), the APS DPP/DAMOP meeting archives, HAL/LNCMI, GitHub, Zenodo, MDSOAR (UMBC), and arXiv.
A source is usable only if it has a **measured** target **and** enough geometry to run the bare solver without guessing.

## Usable now
| source | measured | geometry | status |
|---|---|---|---|
| EPFL bulk-machined coil, arXiv:1901.08791 / SciPost Phys. 6, 048 | 1.28 G/A Hall point, R 10.4(10) mΩ, L 116(2) µH | complete (Sec. 2) | **in tests**: B −0.3 % |
| Claw-ZS, arXiv:2607.02813 + github.com/olsenlab-science/Claw-ZS | Fig. 4a B(z), V, Z(f), switching | `CZS_radia.py` layer table | **in tests**: peak −1.1 % |
| JQI Bitter Ioffe-Pritchard, arXiv:2008.11181 (RSI 92, 033201) + github.com/JQIamo/Bitter-Ioffe-Pritchard (zenodo 10.5281/zenodo.3999624) | Hall probe at 150 A, pair 3.81 cm: curvature B0 -156.94(1) uT/A, B'' -49.4(4) uT/(cm^2 A); anti-bias B0 160.46(4) uT/A, B' -24.2(1) uT/(cm A); R 9.2(6)/13.0(9)/32(2) mOhm; L 5.0(4)/21(1)/68(5) uH | layer table from the authors' design notebook `Cloverleaf_Trap_Reversed_v9.nb` (tests/data/jqi_bitter_ip_design.json); as-built counts not in the paper | **in tests**: shape ratios close (curvature B''/B0 0.3146 vs 0.3148 /cm^2; anti-bias B'/B0 +2.8 %); amplitude 9-12 % solver-high on all four (flagged, not tuned); R anti-bias 13.11 vs 13.0(9) mOhm, curvature 4.90 vs 9.2(6) mOhm (-47 %, leads/contacts) |

## Blocked (missing geometry; nothing guessed)
| source | what is public | missing |
|---|---|---|
| NHMFL 41.5 T all-resistive (Toth & Bole, IEEE TASC 28(3) 2018, doi:10.1109/TASC.2017.2775578) | about 41.4 T at 48 kA, 32 mm bore, 6 coils, 1 m OD (NHMFL highlight PDF) | coil-by-coil radii, grading, plate counts (paywalled paper) |
| HFML 38 T (Radboud repository 2066/161973) | 37.53 T at 40 kA, 20.5 MW, 0.936–0.938 T/kA | as-built A–E coil table (see full report) |
| UMBC ALPHA 7 T Bitter (APS DPP 2017/2023/2024/2025 abstracts) | 15 cm bore, mirror ratio about 1.1, L = 1.888 mH, R reported "0.364 mΩ", within 2 % of theory | variable turn-density table; field map not published yet. Watch for a thesis on MDSOAR. |
| Tohoku HFLSM inner double-Bitter | 19.0 T, 32 mm bore, 7.2 MW | geometry |
| Sabulsky (arXiv:1309.5330) | 1.93(2) G/A at 14 mm, R, L, ΔT | 9 dimensions (see `sabulsky_blocked.md` in the residual-ml results) |

## Not magnet targets
- giaccone/Bfield_measured (GitHub): a 5-turn 1 kHz loop with no axial extent given. It could serve as a loop-model check, but has low value.
- Rtavakol/MagLab (GitHub): shim-coil field data, not Bitter.
- arXiv:2308.15476: a magnetometer paper.
- arXiv:2008.05046: a ⁶Li apparatus paper; the "1.28 G/A" hit actually comes from arXiv:1901.08791.
- EMFL ISABEL D9.4 (2025): pre-designs only.
- Bates/Birmingham/Romero-Talamás IEEE TMAG 53(3) 2017 (MDSOAR): an optimiser paper; its BETA measurements are already in the bundle.

## Code-check material (not data)
- Akhmeteli et al., IEEE TASC 2018 (open copy on public.magnet.fsu.edu): non-circular bore stress (von Mises ×3.18 unsupported, ×1.67 with partial support). A candidate reference if an elliptic-disk solver is added.

## Patents (Google Patents / USPTO / CN / JP / WO)
None of these has a **measured** field. All give design values only, so they can serve at most as cross-code checks.
| patent | content | verdict |
|---|---|---|
| US4660013A (GE, 1987), resistive whole-body MRI magnet | Tables 1–2: six-coil geometry (radii, widths, axial distance, turns, current 685/827 A). Air-core design stated as "0.25 T" | **ambiguous**. If "axial distance" means the coil centre, the middle coils overlap their mirror images (impossible). Reading it as the inner face gives 0.274 T from E3 (+10 % vs the rounded 0.25 T). Not used. |
| US4774487A, US4748429, US4808956, US4736176, US4823101 (Thomson-CGR, 1980s), Bitter-disk MRI magnets | requirements only (0.15–0.5 T, 1–10 ppm on a 40 cm DSV); no dimension table | no data |
| US6876288B2, transverse-field Bitter magnet | concept | no data (already in `non-circular/`) |
| US12191073B2 (KIT, 2020), HTS "Bitter-principle" device | targets only (3 T at 100 A/mm², about 10 cm³) | no data, and not resistive copper |
| JP3706900B2 (Tohoku), YBCO Bitter-type plates | design: 80/150 mm, 40 µm film, 21–23 kA at 20 K | no data, superconducting |
| CN122192447A (CAS Hefei), Bitter water-flow test rig | measurement apparatus | no magnet data |
| CN117236080B (CHMFL) | stress-analysis recipe | method only |

## Industry and facility documents
- HFML: 38 T from five Florida-Bitter coils at 20.5 MW and 40 kA; resistive magnets cooled at 140 L/s, 20 bar. No coil tables.
- LNCMI Grenoble (EMFL pages): 36–37 T in a 34 mm bore, 31 T in 50 mm, and 20 T in 170 mm, each at 24 MW. No coil tables.
- NHMFL cell 6 page and user-committee reports: 41.5 T, 32 mm bore. No coil tables.
- Science Museum record (RRE Malvern / IRD Bitter solenoid): 23 T in a 3 cm bore at 8 MW. Historical.
- FNAL-hosted 1972 Magnet Technology conference papers: Mulhall (IRD), "Bitter magnets for 20 T" (C720919/p360), and Parkinson et al., RRE Malvern (p450). The server returns 403 to automated fetches; worth a manual download, since period papers often print full coil tables.
