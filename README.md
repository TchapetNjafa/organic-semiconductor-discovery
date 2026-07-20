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
  in preparation for the **Journal of Chemical Information and Modeling** (ACS),
  evaluating seven candidates on two independent axes.

---

## What Paper 1 shows

We screen 17,458 molecules with `PCE_SAScore = PCE − SAScore` (Scharber efficiency
minus the RDKit synthetic-accessibility score). Seven molecules score above zero —
but the metric reads only computed descriptors, so the raw "viable" set is not a
candidate list:

- **Three are not organic semiconductors:** molecular oxygen (`977`), a magnesium
  carbonate salt (`11029`), and a hydroquinone–benzoquinone cocrystal (`7801`).
- **Two are reactive azides** (`17851`, `20778`) and one carries a **nitroso**
  group (`4550`).

A chemical-validity filter and a reactive-group screen reduce the viable set
**7 → 4 → 1**. Only molecule `1712` — a stable, conventionally conjugated
molecule — survives, with a modest predicted PCE (8.5%). Two further findings:
the corrected four-point reorganization energy of the azido-triazines is
**negative** (unphysical; a sign of instability), and simple heteroatom-composition
descriptors correlate only **weakly** with frontier-orbital energies (|r| ≤ 0.34).

**The methodological point:** a chemical-validity gate belongs as a standard step
in descriptor-based materials screening, and the gap between the pre- and post-gate
sets is itself a result.

> Note: earlier versions of this repository presented a "top 7 candidates for
> synthesis" narrative and some uncorrected statistics (e.g. reorganization
> energies with a ~0.35 eV mean). Those are **superseded** by the results above.

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
├── paper1_opv_screening/               # CURRENT Paper-1 material
│   ├── scripts/
│   │   ├── filter_genuine_osc.py        # chemical-validity filter
│   │   ├── stability_screen.py          # reactive-group screen
│   │   ├── figstyle.py                  # shared figure style
│   │   └── generate_fig{1,2,3,4,6,7,8}*.py
│   ├── data/                            # genuine_osc_ranked, viable_view_a/b
│   └── figures/                         # the 8 manuscript figures (PDF + PNG)
│
├── paper2_biosensing_screening/        # CURRENT Paper-2 material
│   ├── scripts/
│   │   ├── redock_candidates.py         # whole-protein blind docking (smina)
│   │   ├── run_tddft_scan.py            # TD-DFT solvatochromic scan (ORCA)
│   │   ├── render_pose.py, crop_pose.py # docking-pose figure
│   │   ├── generate_figures{,2}.py, make_si_table.py, figstyle_p2.py
│   ├── data/                            # docking + solvatochromic result CSVs, geometries
│   ├── figures/                         # the 7 manuscript figures (PDF + PNG)
│   ├── docking_outputs/                 # ligands, per-mode logs, docked poses
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
