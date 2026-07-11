# Computational Workflow Documentation

## Overview

This document provides a detailed description of the computational workflow used in the manuscript "Data-Driven Discovery of Synthetically Compatible Organic Semiconductors for Multifunctional Applications."

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

### Stage 6: Molecular Docking Simulations

**Input:** Top candidates with PCE_SAScore > 0

**Process:**

1. **Protein Target Selection:**
   - HIV-1 protease (PDB: 1HVR)
   - Proteasome (PDB: 5LF3)
   - SARS-CoV-2 M^pro (PDB: 6LU7)
   - Aromatase (PDB: 3EQM)
   - Thymidine kinase (PDB: 1KIM)
   - Dihydrofolate reductase (PDB: 1RX2)

2. **Ligand Preparation:**
   - Convert SMILES to 3D structures using RDKit
   - Energy minimization with UFF force field
   - Convert to PDBQT format for AutoDock Vina

3. **Receptor Preparation:**
   - Download PDB structures
   - Remove water molecules and heteroatoms
   - Add polar hydrogens
   - Define binding box around active site

4. **Docking:**
   - AutoDock Vina with exhaustiveness = 8
   - Generate 9 binding poses per ligand
   - Extract best binding affinity (kcal/mol)

**Output:**
- Binding affinity matrix (7 molecules × 6 targets)
- Docking poses and interaction profiles
- Molecular recognition capabilities assessment

**Notebook:** `docking_all_targets_with_download.ipynb`

---

### Stage 7: TDDFT Calculations (Optional)

**Input:** Top candidates

**Process:**

1. **Excited State Calculations:**
   - Time-dependent DFT (TDDFT)
   - B3LYP/6-31G* level of theory
   - Calculate first 10 excited states

2. **Properties Extracted:**
   - Singlet-triplet gap (ΔE_ST)
   - Oscillator strength (f)
   - Excited state lifetime (τ)

**Output:**
- Fluorescence properties
- OLED suitability assessment
- TADF potential identification

**Notebook:** `xtb_crest_validation_report.ipynb`

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

