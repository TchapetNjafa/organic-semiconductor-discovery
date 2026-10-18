# Organic Semiconductor Discovery — gap observables in high-throughput OPV screening

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18201812.svg)](https://doi.org/10.5281/zenodo.18201812)
[![Code: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Supporting code and data for the manuscript *Class-Dependent Errors of Cheap Band-Gap Observables Undermine
High-Throughput Screening of Organic Photovoltaic Absorbers* (P. S. Mvoto Kongo, S. C. Teguia Kouam,
J.-P. Tchapet Njafa, S. G. Nana Engo). Concept DOI of the deposit: 10.5281/zenodo.18201812.

Every script runs from the repository root as `python -I Paper1_AdvTheorySimul/<script>.py .`

## What the paper shows

- Kohn–Sham gaps stored in PubChemQC and sTDA-xTB excitation energies select almost disjoint candidate sets from the
  same 1.1–2.2 eV window (33 vs 559 molecules, 4 shared, 15,952 molecules).
- The raw error depends on the chemical class: sTDA-xTB overestimates 34 donor–acceptor optical gaps by 1.0 eV but
  reproduces four acenes to 0.12 eV. A linear calibration brings each observable to about 0.1 eV within the
  donor–acceptor class and worsens the acenes.
- Disagreement between the two observables flags defective database records: recomputation confirms 21 of 29 flagged
  records and none of 29 controls.
- A validity-first, calibrated cascade keeps 37 of 17,458 records, 10 of them inside the chemical domain of the
  calibration.

## What is here

| Path | Content |
|---|---|
| `ZENODO_UPLOAD/dataset_pubchemqc_opv_17458.csv` | 17,458 PubChemQC records (B3LYP/6-31G*//PM6, block `CHNOPSFClNaKMgCa500`, IDs 1–24,999): SMILES, formula, stored frontier energies and gap |
| `DATA_external/HOPV15/` | HOPV15 v4 (Lopez et al., *Sci. Data* 2016, CC BY 4.0), experimental optical gaps used for calibration |
| `NTO_and_related_study/scripts/00*, 01*` | SMILES → RDKit ETKDG (seed 42) + UFF → xtb4stda/stda (sTDA-xTB) for the PubChemQC shard |
| `Paper1_AdvTheorySimul/build_stda_all.py`, `join_ks_stda.py` | collect S1 and the lowest bright state (f ≥ 0.01) from the raw spectra; join with the stored gaps |
| `Paper1_AdvTheorySimul/hopv_stda.py`, `hopv_stda_dftgeom.py` | sTDA-xTB on HOPV15 (UFF and DFT geometries) |
| `Paper1_AdvTheorySimul/ks_gap_check.py`, `b3lyp_reopt.py`, `std2_b3lyp.py` | PySCF B3LYP recomputations, geomeTRIC re-optimisation, std2 sTDA/XsTDA on a B3LYP ground state |
| `Paper1_AdvTheorySimul/flag_recall.py` | recall sample of the record check: B3LYP/6-31G* recomputation of 120 unflagged records (seed 11) |
| `Paper1_AdvTheorySimul/gate_cascade.py`, `applicability.py`, `pipeline/` | validity rules R1–R6, record check, calibrated window, Scharber gates, applicability domain |
| `Paper1_AdvTheorySimul/make_tables.py`, `make_si_tables.py`, `make_figures.py` | Tables 1–2, Tables S1, S4, S5 and Figures 2–5 from the files in `outputs/extract/` |
| `Paper1_AdvTheorySimul/formal/` | Lean 4 (v4.34.1) + Mathlib proofs F1–F3 (`lake build`) |
| `outputs/extract/` | all intermediate results the manuscript uses (CSV/TXT), sTDA-xTB spectra of HOPV15, std2 outputs |
| `outputs/analysis/analysis-ledger.md` (+ `archive/`) | analysis record: method, verification source, result, interpretation and caveats of every number in the manuscript |
| `Paper1_AdvTheorySimul/manuscript/si/` | Supporting Information (LaTeX source and PDF) |

The raw sTDA-xTB output of the 17,448 shard molecules (`xtb_results/<mol_id>/{mol.xyz,tda.dat}`, 320 MB compressed)
is deposited on Zenodo only, as `xtb_results_shard_tda.tar.gz`. Extract it under `NTO_and_related_study/data/` to
rerun `build_stda_all.py`. std2 Molden inputs (≈1.7 GB) are not deposited; `std2_b3lyp.py` regenerates them.

## Reproducing the tables

```bash
python -I Paper1_AdvTheorySimul/gate_cascade.py .     # outputs/extract/gate_cascade.{csv,txt}
python -I Paper1_AdvTheorySimul/applicability.py .    # outputs/extract/applicability.{csv,txt}
python -I Paper1_AdvTheorySimul/make_tables.py .      # Table 1, Table 2
python -I Paper1_AdvTheorySimul/make_si_tables.py .   # SI tables
```

## Software

Analysis and tables: Python 3.12, RDKit 2023.09.3, pandas 2.1.4, NumPy 1.26.4, SciPy 1.11.4, Matplotlib 3.6.3.
Quantum chemistry: PySCF 2.13.0, geomeTRIC 1.0.2, xtb4stda 1.0, stda 1.6.3, std2 2.0.1 (environment with RDKit
2026.03.2, pandas 3.0.2, NumPy 2.3.5). Formal proofs: Lean 4 v4.34.1 with Mathlib.

## Reproducibility note

The Lean project pins Mathlib to revision `d13f23b` (tag v4.34.1): run `lake exe cache get && lake build` in
`Paper1_AdvTheorySimul/formal/`. The proofs were checked locally against that Mathlib through a path dependency; the
pinned git form in this deposit has not been rebuilt from scratch.

The three `NTO_and_related_study/scripts/` files are deposited exactly as run and contain absolute paths of the
original machine; edit those paths before rerunning them.

Rerunning the scripts on this folder reproduces Tables 1–2 and S1, S4, S5 byte for byte. `applicability.csv` was
produced with RDKit 2026.03.2; with RDKit 2023.09.3 one record of the random 2,000-record database sample receives a
different sp² count, which changes no reported number. Other CSVs agree to floating-point precision (≤ 1e-15).

## Archives

`Archives/` holds material from earlier analyses of the same PubChemQC shard that the paper does not use: the
PCE_SAScore screen (`paper1_opv_screening/`), the biosensing study (`paper2_biosensing_screening/`), notebooks,
scripts and figures of the first project phase, and the previous README (`README_previous.md`) and citation file. It
is kept for traceability only; nothing in the paper depends on it.

## Licences

Code: MIT (`LICENSE`). Data produced in this work (`outputs/`): CC BY 4.0. Third-party data keep their own terms:
HOPV15 (CC BY 4.0, cite Lopez et al. 2016, 10.1038/sdata.2016.86); PubChemQC records (cite Nakata & Maeda,
*J. Chem. Inf. Model.* 2023, 10.1021/acs.jcim.3c00899).
