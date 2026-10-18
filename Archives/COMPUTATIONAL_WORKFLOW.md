# Computational Workflow Documentation

## Overview

This document provides a detailed description of the computational workflow used
across the two-paper programme built on this deposit: Paper 1 (OPV screening,
Stages 1-5, submitted to RSC Digital Discovery) and Paper 2 (biosensing
two-axis evaluation, Stages 6-7, in preparation for the Journal of Chemical
Information and Modeling).
Stages 6 and 7 below describe the docking and TD-DFT protocols **actually used
in Paper 2**; an earlier six-target/TADF-OLED docking plan was explored and
abandoned (see `paper2_biosensing_screening/DOCKING_METHODS_NOTE.md` and
`TDDFT_METHODS_NOTE.md` for the full protocol and the reasons for the change).

## Workflow Stages

### Stage 1: Data Acquisition and Preprocessing

**Input:** PubChemQC B3LYP/6-31G\*//PM6 dataset, `CHNOPSFClNaKMgCa500` subset
(Nakata & Maeda 2023, doi:10.1021/acs.jcim.3c00899). The subset name encodes its
scope: molecules composed of C, H, N, O, P, S, F, Cl, Na, K, Mg, Ca with mass up
to ~500 Da, provided as neutral closed-shell singlets.

**Process:**
1. Download molecular data from the PubChemQC subset above.
2. Extract HOMO, LUMO, and energy-gap values.
3. Assemble the working set: neutral, closed-shell singlet molecules with
   mass ≲ 500 Da. **No heteroatom requirement and no lower mass bound were
   applied** — the working set includes heteroatom-free molecules (benzene,
   methane, etc.). This permissive, chemistry-agnostic selection is deliberate;
   Paper 1 shows it is part of why the raw metric surfaces artifacts.
4. Total molecules in the working set: 17,458.

**Output:** Structured dataset with quantum-chemical descriptors
(`dataset_pubchemqc_opv_17458.csv`). Some original working notebooks/files use
legacy `GDB9`/`qm9` naming — the data are PubChemQC, not GDB-9/QM9.

**Notebook:** `dataset_extraction.ipynb`

---

### Stage 2: PCE Calculation (Scharber Model)

**Input:** HOMO, LUMO energies from Stage 1

**Process:**

1. **Frontier Molecular Orbital (FMO) Alignment:**
   - Donor candidates: HOMO > HOMO_PCBM and LUMO > LUMO_PCBM
   - Acceptor candidates: HOMO < HOMO_PCDTBT and LUMO < LUMO_PCDTBT
   
   Reference materials:
   - PCBM: HOMO = -6.1 eV, LUMO = -3.7 eV
   - PCDTBT: HOMO = -5.5 eV, LUMO = -3.6 eV

2. **Open-Circuit Voltage (V_OC):**
   ```
   V_OC = (1/e) × (|HOMO_donor| - |LUMO_acceptor|) - 0.3 V
   ```

3. **Short-Circuit Current Density (J_SC):**
   ```
   J_SC = ∫ P(λ) × EQE(λ) dλ
   ```
   where EQE is assumed constant at 0.65

4. **Fill Factor (FF):**
   Dynamic calculation based on V_OC:
   ```
   FF = (V_OC - ln(V_OC + 0.72)) / (V_OC + 1)
   ```

5. **Power Conversion Efficiency (PCE):**
   ```
   PCE = (V_OC × J_SC × FF) / P_in
   ```
   where P_in = 1000 W/m² (AM1.5G solar spectrum)

**Output:** 
- 17,334 potential donors (vs. PCBM)
- 16 potential acceptors (vs. PCDTBT)

**Notebook:** `scharber_pce_calculation.ipynb`

---

### Stage 3: Synthetic Accessibility Evaluation

**Input:** SMILES strings for all candidates

**Process:**

1. **SAScore Calculation:**
   - Uses RDKit implementation of Ertl & Schuffenhauer algorithm
   - Score range: 1 (easy) to 10 (difficult)
   - Based on fragment contributions and complexity penalties

2. **PCE_SAScore Metric:**
   ```
   PCE_SAScore = PCE - SAScore
   ```

3. **Filtering:**
   - Retain molecules with PCE_SAScore > 0

**Output (raw):**
- 3 molecules positive as donors (vs. PCBM) + 4 positive as acceptors (vs. PCDTBT)
- 7 unique molecules with PCE_SAScore > 0
- **Caution:** this raw set is dominated by chemically implausible species and is
  NOT a candidate list. It is filtered in Stage 3b.

**Notebook:** `screening_pce_sascore.ipynb`

---

### Stage 3b: Chemical-validity and stability screening (Paper 1)

**Input:** the 7 raw molecules with PCE_SAScore > 0 (and the full ranking)

**Process:**
1. **Chemical-validity filter** (`filter_genuine_osc.py`): retain a molecule only
   if it is a single covalent species (no salts/cocrystals), contains no metal,
   has at least one aromatic ring, ≥6 heavy atoms, and ≥6 conjugated atoms.
2. **Stability screen** (`stability_screen.py`): flag/remove reactive groups
   (azide, nitroso, diazo) via SMARTS.

**Output:**
- Raw viable set (7) → chemically valid (4: 17851, 20778, 4550, 1712)
  → stable (1: molecule 1712).
- Removed as artifacts: molecular oxygen (977), magnesium carbonate (11029),
  quinhydrone cocrystal (7801). Removed as reactive: azides 17851, 20778; nitroso 4550.

**Scripts:** `paper1_opv_screening/scripts/filter_genuine_osc.py`,
`stability_screen.py`

---

### Stage 4: Sensitivity Analysis

**Input:** PCE and SAScore values for all molecules

**Process:**

1. **Weighting Scenarios:**
   - 1×PCE - 1×SA (original)
   - 2×PCE - 1×SA (emphasize efficiency)
   - 1×PCE - 2×SA (emphasize synthesis)
   - 3×PCE - 1×SA
   - 1×PCE - 3×SA
   - 1.5×PCE - 1×SA
   - 1×PCE - 1.5×SA

2. **Analysis:**
   - Count candidates with score > 0 for each weighting
   - Identify top 7 candidates for each scenario
   - Assess stability of rankings

**Output:** 
- Sensitivity analysis results table
- Visualization of candidate stability
- Validation of 1:1 weighting choice

**Script:** `pce_sascore_sensitivity_analysis.py`

---

### Stage 5: Molecular Structure Analysis

**Input:** SMILES strings for top 7 candidates

**Process:**

1. **Functional Group Identification:**
   - Aromatic rings (benzene, thiophene, furan, pyridine)
   - Functional groups (carbonyl, carboxyl, amino, phenol, nitrile, azide)
   - Heteroatom counting

2. **Conjugation Analysis:**
   - Aromatic atom fraction
   - Aromatic bond count
   - Rotatable bonds (rigidity measure)

3. **Physicochemical Properties:**
   - Molecular formula
   - Heavy atom count
   - LogP, TPSA, QED

**Output:**
- Structure-property analysis table
- Molecular structure visualizations
- Functional group distribution

**Script:** `molecular_structure_analysis.py`

---

### Stage 6: Molecular Docking Simulations (Paper 2, as actually run)

**Input:** The four Paper-1 chemical-validity-filter survivors (1712, 17851,
20778, 4550) — no new candidate selection was performed for Paper 2.

**Process:**

1. **Protein target selection:** whole-protein blind docking (no assumed
   active site) against four structurally distinct targets: HIV-1 protease
   (PDB 1DMP), Hsp90 (2XJX), a neurodegenerative-disease target (1SYH), and
   the SARS-CoV-2 main protease (6Y2F). Targets were chosen to sample distinct
   binding environments, not to represent a single disease hypothesis.
2. **Ligand preparation:** RDKit `ETKDGv3` embedding (`randomSeed=42`) + UFF
   optimization on the reused B3LYP/6-31G* geometries, converted to PDBQT in a
   single Open Babel pass (Gasteiger charges); no second 3D-embedding step, so
   the seeded conformer is preserved.
3. **Receptor preparation:** prepared PDBQT files used as provided; search box
   defined from each receptor's heavy-atom bounding box plus 8 Å padding
   (rigid receptor, no induced fit).
4. **Docking:** `smina` (AutoDock Vina 1.1.2 fork), fixed seed = 42,
   exhaustiveness = 32, nine output modes; best-scoring mode reported per pair.
   Two independent runs reproduced all twenty-eight ligand-target affinities
   bit-identically.

**Output:**
- Best-mode affinity matrix (7 molecules × 4 targets, 28 pairs). Docking separates
  the set by stability class: all four stable donors (1712, 9168, 17574, 18506)
  bind within the native-inhibitor range ($-6.4$ to $-8.0$ kcal/mol, the
  all-carbon aromatic 9168 strongest); the three reactive structures (4550,
  17851, 20778) bind more weakly ($-4.2$ to $-6.4$ kcal/mol).
- **Reference-ligand calibration:** each target's native co-crystal ligand was
  redocked under the identical protocol (`calibrate_reference_ligands.py`),
  giving a known-inhibitor band of $-5.5$ to $-7.7$ kcal/mol. The blind search
  does not recover the crystallographic poses (RMSD 16–34 Å), so affinities are
  a coarse, relative binding-potential scale, not site-specific predictions.

**Scripts:** `paper2_biosensing_screening/scripts/redock_candidates.py` (original
four), `dock_new_candidates.py` (three added stable donors, identical protocol),
`calibrate_reference_ligands.py` (native-ligand calibration).
Full protocol and limitations: `paper2_biosensing_screening/DOCKING_METHODS_NOTE.md`.

> An earlier six-target plan (HIV-1 protease 1HVR, Proteasome 5LF3, SARS-CoV-2
> M$^{\mathrm{pro}}$ 6LU7, aromatase, thymidine kinase, DHFR; exhaustiveness=8)
> was explored and abandoned: the Proteasome target (7PG9) proved intractable
> for whole-protein docking at ~45,000 atoms, and two receptor conversions
> (4LDE, 5JWA) failed. The four-target protocol above is the one used in the
> manuscript.

---

### Stage 7: TD-DFT Solvatochromic Shift (Paper 2, as actually run)

**Input:** All seven candidates. The four original molecules reuse their
Paper-1 B3LYP/6-31G* geometries; the three added stable donors (17574, 18506,
9168) were geometry-optimized at the same B3LYP/6-31G* level
(`optimize_new_geometries.py`) so all seven share a common footing.

**Process:**

1. **Excited-state method:** TDA-TD-DFT (Tamm-Dancoff approximation), eight
   lowest singlet roots, computed with ORCA. Full linear-response TD-DFT was
   attempted first but failed to run reliably (CIS module instability on a
   converged reference); TDA was used throughout as the more stable choice.
2. **Functionals:** B3LYP/6-31G* (matches the Paper-1 screening level) and
   CAM-B3LYP/6-31G* as a range-separated cross-check, since global hybrids are
   known to underestimate charge-transfer excited-state energies.
3. **Environment:** CPCM implicit solvent, two solvents per molecule —
   toluene (nonpolar, organic-blend proxy) and water (polar, aqueous-biosensing
   proxy).
4. **Solvatochromic descriptor:** $\Delta E_{S_1} = E_{S_1}(\text{water}) -
   E_{S_1}(\text{toluene})$, computed at both functional levels.

**Output:**
- Molecule 1712 (low dipole, 2.56 D): shift within numerical noise and changes
  sign between functionals (+0.006 eV B3LYP, -0.071 eV CAM-B3LYP) — direction
  not resolved, magnitude far smaller than the other three.
- The three high-dipole molecules (7.42-9.74 D): large, same-sign shifts under
  both functionals (e.g. 20778: +0.58 eV B3LYP, +1.10 eV CAM-B3LYP).
- **Caveat:** under CAM-B3LYP the tracked S$_1$ state is nearly oscillator-dark
  for 17851 and 20778 ($f \sim 10^{-5}$-$10^{-6}$), which bears on how directly
  the shift could be exploited as an absorption-based optical readout.

**Script:** `paper2_biosensing_screening/scripts/run_tddft_scan.py`. Full
eight-state tables: `paper2_biosensing_screening/data/solvatochromic_results_*.csv`
and `tab_S1_absorption_states.tex`. Raw ORCA output (146 MB) is archived
separately as `paper2_biosensing_screening/raw_tddft_output_B3LYP_CAMB3LYP.zip`
in the Zenodo deposit (not included in the GitHub copy of this folder). Full
protocol: `paper2_biosensing_screening/TDDFT_METHODS_NOTE.md`.

> An earlier plan targeted TADF/OLED photophysics (singlet-triplet gap,
> excited-state lifetime). This was abandoned: no experimental or computed
> S$_1$/T$_1$ data existed for the real candidate set to support that framing.
> The solvatochromic-shift protocol above is the one used in the manuscript.

---

### Stage 8: Clustering Analysis

**Input:** All candidates with calculated properties

**Process:**

1. **Feature Engineering:**
   - ECFP4 fingerprints (Extended-Connectivity Fingerprints, radius 2, 2048 bits)
   - Molecular descriptors (V_OC, J_SC, PCE, SAScore)

2. **Dimensionality Reduction:**
   - PCA (Principal Component Analysis)
   - t-SNE (t-Distributed Stochastic Neighbor Embedding)

3. **Clustering:**
   - K-means algorithm
   - Optimal k determined by silhouette score, Calinski-Harabasz index, Davies-Bouldin index

4. **Analysis:**
   - Cluster validation metrics
   - Property distributions by cluster
   - Identification of high-performance clusters

**Output:**
- Cluster assignments for all molecules
- Visualization of chemical space
- Structure-property relationships by cluster

**Notebook:** `screening_pce_sascore.ipynb`

---

## Software and Tools

### Core Dependencies
- **Python:** 3.8+
- **RDKit:** 2022.03.1+ (cheminformatics)
- **NumPy:** 1.21.0+ (numerical computing)
- **Pandas:** 1.3.0+ (data manipulation)
- **Scikit-learn:** 1.0.0+ (machine learning)
- **Matplotlib/Seaborn:** Visualization

### Specialized Tools
- **AutoDock Vina:** 1.2.0+ (molecular docking)
- **Open Babel:** 3.1.1+ (chemical file format conversion)

### Computational Resources
- **CPU:** Standard desktop/laptop (no GPU required)
- **RAM:** 8 GB minimum, 16 GB recommended
- **Storage:** ~500 MB for datasets and results

---

## Reproducibility Notes

1. **Random Seeds:** All stochastic processes (k-means, t-SNE) use fixed random seeds for reproducibility
2. **Software Versions:** Exact versions specified in `requirements.txt`
3. **Data Provenance:** All input data traceable to PubChemQC database
4. **Parameter Documentation:** All hyperparameters explicitly stated in notebooks

---

## Validation and Quality Control

1. **PCE Model Validation:** Compared against Harvard CEP benchmarks
2. **SAScore Validation:** Verified against known synthetic complexity examples
3. **Docking Validation:** Cross-checked with experimental binding data where available
4. **Clustering Validation:** Multiple metrics (silhouette, CH, DB) used

---

## Contact

For questions about the computational workflow:
- **Email:** sorrel.mvoto@facsciences-uy1.cm
- **GitHub Issues:** [Repository URL]/issues

