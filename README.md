# Organic Semiconductor Discovery

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18201812.svg)](https://doi.org/10.5281/zenodo.18201812)
[![Code: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Data and code for a two-paper programme built on one high-throughput screen of the
PubChemQC database.

- **Paper 1 (organic photovoltaics):** *Efficiency and Accessibility Are Not
  Enough: Chemical-Validity Filtering in High-Throughput Screening of Organic
  Photovoltaic Materials* — submitted to **RSC Digital Discovery**.
- **Paper 2 (biosensing / multifunctional):** *Predicted Protein Binding and
  Optical Environment-Sensitivity Select Different Subsets of
  Organic-Semiconductor Candidates: A Multi-Objective Computational Evaluation* —
  submitted to **RSC Journal of Materials Chemistry C**,
  evaluating seven candidates on two independent axes.

---

## What Paper 1 shows

We screen 17,458 molecules with `PCE_SAScore = PCE − SAScore` (Scharber efficiency
minus the RDKit synthetic-accessibility score). The four physical gates admit **25
molecule–reference pairings** under the complementary absorber convention and 16
under the strict one — but the metric reads only computed descriptors, so what it
admits is not a candidate list:

- **14 of the 25 are not organic semiconductors**, and 11 of those are
  disconnected ion pairs: calcium thiosulfate (`24964`) ranks second under one
  convention, magnesium carbonate (`11029`) third under the other, plus a sodium
  dinitrophenolate, a folate disodium salt and six quaternary ammonium or iminium
  chlorides.
- **Molecular oxygen** (`977`, 1.965 eV) and a **quinhydrone cocrystal** (`7801`,
  1.467 eV) both absorb inside the useful window and are removed only by
  structural rules. A band-gap pre-filter removes **no** artifact — the energetic
  and structural criteria are independent. O₂'s tabulated gap is a spin-restricted
  singlet, not its triplet ground state.
- **Two admitted molecules carry a nitroso group**: `4550`, which ranks *first*
  among all admitted pairings under every convention tested, and `12919`.

The cascade is **25 → 11 → 9 → 4** (admitted → chemically valid → passes the
reactive-group screen → above threshold), or 16 → 5 → 4 → 1 under the strict
convention. **First place is convention-dependent** — 9168 / 9168 / 19598 / 22936
across the four (absorber × fill-factor) combinations — and only molecule `9168`
(C₂₆H₁₆, an unsubstituted polycyclic aromatic hydrocarbon) is positive under all
four, with a small margin.

Three further findings:

- **Molecule `1712` is not admitted.** Its −3.823 eV LUMO gives −0.123 eV of
  driving force against the PCBM reference. Earlier versions of this repository
  described it as the single surviving candidate; that is superseded.
- Three of the four molecules above threshold are ***o*-quinone tautomers**
  (`19598`, `21736`, `23450`) and the fourth is an **anthraquinone sulfonic-acid
  dye** (`22936`). All pass every rule in the filter.
- **The screen admits none of the nine established organic semiconductors in the
  same shard** (benzene through pentacene): computed Kohn–Sham gaps of
  2.397–6.781 eV against pentacene's ≈1.8 eV measured optical gap. This bounds how
  much weight the admitted set can carry.

The corrected four-point reorganization energy of the azido-triazines is
**negative** (unphysical; a sign of instability under ionization), and heteroatom
composition correlates only **weakly** with frontier-orbital energies
(|r| ≤ 0.34, n = 17,457).

**The methodological point:** a chemical-validity gate belongs as a standard step
in descriptor-based materials screening, and the gap between the pre- and post-gate
sets is itself a result.

### Reproducing Paper 1

```bash
cd paper1_opv_screening
python -m pip install -r requirements.txt
# the 12 MB property table is Zenodo-only; see data/README_DATA.txt
bash run_all.sh
```

`run_all.sh` regenerates every Paper-1 number, figure and table. Both filter
scripts re-derive their rules from the SMILES and assert agreement with the
pipeline's own columns on all 17,458 rows, so a broken edit fails loudly rather
than silently changing a result. No supervised model is trained anywhere in the
pipeline.

> Note: earlier versions of this repository presented a "top 7 candidates for
> synthesis" narrative and some uncorrected statistics (e.g. reorganization
> energies with a ~0.35 eV mean, and efficiencies from a photovoltaic stage whose
> incident power was 900.1 W m⁻² with no E_CT ≤ E_g check). Those are
> **superseded** by the results above.

---

## What Paper 2 shows

Seven candidates — the four admitted by the Paper-1 validity filter (stable donor
`1712` + reactive foils `17851`, `20778`, `4550`) plus three further stable donors
(`17574`, `18506`, `9168`) added on a distinct stability axis — are evaluated on
two independent axes. **Whole-protein blind docking** (smina, seed 42,
exhaustiveness 32) against four structurally distinct targets — HIV-1 protease
(1DMP), Hsp90 (2XJX), a neurodegenerative target (1SYH), SARS-CoV-2 main protease
(6Y2F) — separates the set by stability class: all four stable donors bind within
the native-inhibitor range ($-6.4$ to $-8.0$ kcal/mol, all-carbon aromatic `9168`
strongest); the three reactive structures bind weaker ($-4.2$ to $-6.4$ kcal/mol).
**TD-DFT solvatochromic shift** (TDA-TD-DFT, B3LYP + CAM-B3LYP/6-31G*, CPCM
toluene vs. water) gives an unrelated ranking: the stable donors barely move with
solvent polarity, while a large shift appears only in the reactive high-dipole
structures — and there in optically dark transitions. The two readouts select
different subsets and are reported side by side, not merged.

A **reference-ligand calibration** (native co-crystal ligands redocked under the
same protocol) places the stable donors within the $-5.5$ to $-7.7$ kcal/mol
known-inhibitor band; the blind search does not recover crystallographic poses
(RMSD 16–34 Å), so affinities are a coarse relative scale, not site-specific. See
`paper2_biosensing_screening/` for the full protocol and results; the raw ORCA
output is on the Zenodo deposit only.

---

## Data provenance

- **Source:** PubChemQC B3LYP/6-31G\*//PM6 dataset, `CHNOPSFClNaKMgCa500` subset
  (Nakata & Maeda 2023, doi:10.1021/acs.jcim.3c00899).
- **Working set:** 17,458 neutral, closed-shell singlet molecules with mass ≲500 Da.
  No heteroatom requirement and no lower mass bound were applied.
- Filenames/notebooks use legacy `GDB9`/`qm9` naming; the data are **PubChemQC**,
  not GDB-9/QM9 (they contain S, Cl, P, metals, and up to 38 heavy atoms).

---

## Repository structure

```
organic-semiconductor-discovery/
├── README.md
├── CITATION.md
├── LICENSE                              # MIT (code)
├── requirements.txt
├── COMPUTATIONAL_WORKFLOW.md            # end-to-end methodology
│
├── paper1_opv_screening/               # CURRENT Paper-1 material, self-contained
│   ├── run_all.sh                       # regenerates every number and artifact
│   ├── scharber_recompute.py            # screening stage, gates G1-G4
│   ├── filter_genuine_osc.py            # chemical-validity filter R1-R6
│   ├── stability_screen.py              # reactive-group screen R7
│   ├── regen_tables.py                  # manuscript tables
│   ├── regen_figures.py, regen_fig6.py, regen_fig7.py
│   ├── figstyle.py                      # shared figure style
│   ├── requirements.txt                 # 5 pinned versions, Python 3.12.3
│   ├── INPUT_CHECKSUM.txt               # SHA-256 of the PubChemQC shard archive
│   ├── MANIFEST.txt                     # SHA-256 + bytes for every file
│   ├── data/                            # ASTM G173 spectrum, corrected lambda table
│   │                                    # (12 MB property table: Zenodo only)
│   ├── outputs/recompute/               # screen_full.csv, summary, robustness, log
│   ├── figures/                         # manuscript figures (PDF; PNG on Zenodo)
│   └── tables/                          # manuscript tables (LaTeX)
│
├── paper2_biosensing_screening/        # CURRENT Paper-2 material
│   ├── scripts/
│   │   ├── redock_candidates.py         # whole-protein blind docking (smina)
│   │   ├── run_tddft_scan.py            # TD-DFT solvatochromic scan (ORCA)
│   │   ├── render_pose.py, crop_pose.py # docking-pose figure
│   │   ├── generate_figures{,2}.py, make_si_table.py, figstyle_p2.py
│   │   ├── dock_gap_candidates_apo.py   # apo redock, 3 gap candidates
│   │   ├── multiseed_variance.py        # seed-to-seed docking spread
│   │   ├── run_tddft_diffuse_basis.py   # 6-31+G* diffuse-basis rerun
│   │   ├── run_tddft_scan_gap_candidates.py
│   │   ├── optimize_gap_geometries.py   # ORCA opt, 3 gap candidates
│   │   ├── correlation_dipole_shift.py  # dipole vs solvatochromic shift, Spearman
│   │   ├── oscillator_window_bound.py   # f-ceiling of the 8-root window
│   │   └── generate_figures_JMC-C.py, generate_figures_modern-V2.py
│   ├── data/
│   │   ├── geometries/                  # 7 optimized B3LYP/6-31G* structures
│   │   └── geometries_gap/              # 3 gap candidates (5081, 20549, 23424)
│   ├── docking_outputs/                 # ligands, docked poses (PDBQT)
│   │   └── poses_multiseed/             # 84 poses from the multi-seed run
│   └── {DOCKING,TDDFT}_METHODS_NOTE.md  # full protocols
│       # raw ORCA output (145 MB) is on the Zenodo deposit only
│
├── notebooks/
│   ├── screening_pce_sascore.ipynb      # screening + PCE_SAScore workflow
│   └── scharber_pce_calculation.ipynb   # Scharber PCE calculation
│
├── scripts/                            # supporting analyses (transport, motifs)
├── figures/                            # legacy/supporting figures
└── data/                               # master dataset + analyses (via Zenodo)
```

The master dataset (`dataset_pubchemqc_opv_17458.csv`, ~12 MB) is distributed via
the Zenodo deposit; place it under `data/` to run the pipeline.

### What is here, and what is on Zenodo

This repository is the **code** record. Its `.gitignore` deliberately excludes bulk
data (`*.csv`, `*.log`) and rendered figures (`*.pdf`, `*.png`), so the result tables,
per-mode docking logs and manuscript figures are **not** tracked here even when they
exist in a local working copy. What is tracked: the analysis scripts, the optimized
molecular geometries (`.xyz`) and the docked poses (`.pdbqt`).

The **complete** record — every CSV, every docking log, the ORCA optimization outputs
and the figures — is archived on Zenodo under concept DOI
[10.5281/zenodo.18201812](https://doi.org/10.5281/zenodo.18201812). Reproducing the
manuscript's numbers requires the Zenodo archive; this repository alone is not
sufficient. Cite the DOI when the data is what matters.

---

## Quick start

```bash
git clone https://github.com/TchapetNjafa/organic-semiconductor-discovery.git
cd organic-semiconductor-discovery
pip install -r requirements.txt

# reproduce the Paper-1 screening result (7 -> 4 -> 1)
cd paper1_opv_screening/scripts
python filter_genuine_osc.py       # -> data/genuine_osc_ranked.csv
python stability_screen.py         # -> data/viable_view_a/b_*.csv

# regenerate the manuscript figures
python generate_fig2_scatter.py    # etc.

# reproduce the Paper-2 docking + TD-DFT results
cd ../../paper2_biosensing_screening/scripts
python redock_candidates.py        # -> ../data/docking_results_clean.csv
python run_tddft_scan.py           # -> ../data/solvatochromic_results_*.csv
```

---

## Citation

See [CITATION.md](CITATION.md). Cite the Zenodo DOI for the data/code and the
relevant paper for the analysis.
