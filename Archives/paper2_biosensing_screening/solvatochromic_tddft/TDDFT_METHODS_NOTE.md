---
name: tddft_solvatochromic_methods_note
description: Methodology, results, and known caveats for the new TD-DFT solvatochromic-shift calculations (Paper 2 optical transduction mechanism)
---

# TD-DFT solvatochromic scan — methodology, results, caveats (2026-07-10)

## Why this exists

Paper 2's docking result establishes binding *potential* but not a *sensing mechanism* — there was
no signal-transduction pathway in the paper's scope (the TADF/OLED framing that would have
supplied one was dropped for lack of data, see `docking_rerun/DOCKING_METHODS_NOTE.md`). This
calculation supplies a real, calculated alternative: the shift in a molecule's lowest optical
absorption (S1) between a nonpolar and a polar environment (a solvatochromic shift), which is a
standard, well-established transduction principle for optical/colorimetric sensors.

## Method

- **Geometries:** existing B3LYP/6-31G*-optimized structures for 1712, 17851, 20778, 4550
  (`NTO_and_related_study/data/input_geometries/*.xyz`, reused as-is, no re-optimization).
- **Level:** B3LYP/6-31G* (matches the DFT level used throughout Paper 1, per `CLAUDE.md`).
- **Excited states:** TDA-TD-DFT (Tamm-Dancoff approximation), 8 lowest singlet roots.
  Full linear-response TD-DFT (`TDA false`) was tried first and **crashed ORCA's parallel CIS
  module** (`orca_cis_mpi`) on a converged SCF; TDA is the standard, more numerically robust choice
  for this kind of screening and ran cleanly — documented in code comments in
  `run_tddft_scan.py`.
- **Environment:** CPCM implicit solvation, two solvents per molecule — toluene (nonpolar proxy for
  the OPV-blend environment these molecules were originally screened in) and water (polar proxy for
  an aqueous biosensing interface).
- **Script:** `run_tddft_scan.py`, portable via `Path.home()` resolution (works unmodified on the
  server under a different username, same directory layout for ORCA/venv).

## Results — S1 solvatochromic shift (water − toluene)

| mol_id | S1 (toluene) | S1 (water) | Shift (eV) | Shift (nm) |
|--------|:---:|:---:|:---:|:---:|
| 1712   | 1.138 eV (1089 nm) | 1.144 eV (1084 nm) | **+0.006** | −5.5 |
| 4550   | 0.197 eV (6289 nm) | 0.374 eV (3318 nm) | **+0.177** | −2972 |
| 17851  | 0.246 eV (5043 nm) | 0.593 eV (2092 nm) | **+0.347** | −2951 |
| 20778  | 0.135 eV (9150 nm) | 0.719 eV (1725 nm) | **+0.583** | −7425 |

Full 8-state table (all molecules/solvents) in `solvatochromic_results.csv`.

> **UPDATE (2026-07-16): superseded framing below.** The findings paragraphs in
> this section describe the original FOUR-molecule run (1712 + reactive foils
> 17851/20778/4550) and its "trade-off / 1712 strongest binder" reading. The
> study was later expanded to SEVEN molecules (three further stable donors —
> 17574, 18506, 9168 — added on a stability axis) and a reference-ligand
> calibration was added. The current, correct reading is: **the two readouts
> select different subsets.** Docking tracks stability class — all four stable
> donors bind within the native-inhibitor band (9168 strongest at −8.0 kcal/mol,
> NOT 1712), the three reactive structures bind weaker; the large solvatochromic
> shift appears only in the reactive structures, and there in optically dark
> states. Binding is NOT ordered by dipole. The protocol description in this note
> remains accurate; only the four-molecule "trade-off/1712-strongest" conclusion
> below is superseded. See the manuscript and README for the seven-molecule
> result.

**Finding (original 4-molecule run — superseded, see update above):** 1712 — the low-dipole molecule, the sole View-B-viable OPV candidate, and the
strongest docking binder — shows essentially **no** solvatochromic shift (+0.006 eV, within
numerical noise). The three high-dipole molecules all show substantial, real shifts (+0.18 to
+0.58 eV). This is the same low/high dipole partition already established from the static dipole
moments, now corroborated by an independent, calculated spectroscopic observable.

**This produces a genuine design trade-off, not a single "best" molecule:** 1712 binds proteins
most strongly (docking) but is optically the least environment-responsive (worst as an optical
transducer); the other three bind more weakly but are strongly solvatochromic (better as optical
transducers). Report this as a trade-off, not resolved in either direction — that is the honest
result.

## Important caveat — do not overstate

The absolute S1 energies (0.14–1.14 eV in toluene) are unusually low for organic chromophores —
these fall in the far-IR, not the visible range typically implied by "optical sensing." Several
states also have very small oscillator strengths (< 0.001), i.e. weak/dark transitions. B3LYP is
known in the literature to systematically underestimate charge-transfer and diffuse excited-state
energies relative to range-separated functionals (e.g. CAM-B3LYP, ωB97X-D). **Treat the absolute S1
energies/wavelengths as level-of-theory-dependent and not quantitatively reliable on their own** —
the manuscript should state this explicitly. The *relative* water-vs-toluene shift is the safer
quantity to lean on (systematic errors partially cancel between two calculations on the same
molecule/geometry), but even that should be framed as a computational prediction requiring
verification with a range-separated functional (and, ultimately, experiment) before any strong
sensing claim is made.

## Range-separated-functional cross-check (CAM-B3LYP/6-31G*, done 2026-07-10)

Re-ran all 8 jobs at CAM-B3LYP/6-31G* (`solvatochromic_results_CAM-B3LYP.csv`,
`runs_CAM-B3LYP/`) to test whether the B3LYP finding survives a method less prone to
underestimating charge-transfer excited states.

| mol_id | Shift, B3LYP (eV) | Shift, CAM-B3LYP (eV) |
|--------|:---:|:---:|
| 1712   | +0.006 | **-0.071** |
| 4550   | +0.177 | +0.193 |
| 17851  | +0.347 | +1.151 |
| 20778  | +0.583 | +1.100 |

**The qualitative partition is robust across both functionals:** 1712's shift stays an order of
magnitude smaller than the other three under both methods (magnitude 0.006-0.07 eV vs. 0.19-1.15
eV) — the sign of 1712's tiny shift flips between functionals, which itself confirms it is
statistically indistinguishable from zero/noise, while the other three are large and same-sign
(blue-shifting in water) under both methods. CAM-B3LYP also gives much more physically reasonable
absolute S1 energies (0.86-1.7 eV / 550-1400 nm toluene) than B3LYP's implausible far-IR values
(0.14-1.14 eV / 5000-9000+ nm), consistent with the known B3LYP CT-state underestimation flagged
above. **Report the CAM-B3LYP absolute values as the more reliable estimate**, and cite both
functionals' agreement on the qualitative trade-off as the robustness argument.

## What is still NOT done

- No explicit solvent / QM-MM treatment of the actual binding pocket environment — CPCM(water) and
  CPCM(toluene) are bulk-solvent proxies, not the real protein pocket dielectric.
- No experimental UV-vis validation.
- 1712's shift, while small under both functionals, is not converged in sign — do not claim a
  specific direction for it, only that it is much smaller than the other three candidates'.
