# Data-Driven Discovery of Synthetically Compatible Organic Semiconductors

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18201813.svg)](https://doi.org/10.5281/zenodo.18201813)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the complete computational workflow, analysis scripts, and Jupyter notebooks for the manuscript:

**"Data-Driven Discovery of Synthetically Compatible Organic Semiconductors for Multifunctional Applications: A Computational Workflow"**

Submitted to: *Computational Materials Science* (Elsevier)

### Key Features

- High-throughput screening of 17,458 organic molecules from PubChemQC database
- Integrated PCE_SAScore metric balancing photovoltaic efficiency and synthetic accessibility
- Large-scale structural validation (1,000 molecules) with statistical significance (p < 0.001)
- Multifunctional analysis for bio-optoelectronic applications
- FAIR-compliant data and code

## Repository Structure

```
.
├── README.md                           # This file
├── LICENSE                             # MIT License
├── requirements.txt                    # Python dependencies
├── COMPUTATIONAL_WORKFLOW.md           # Detailed methodology
├── scripts/                            # Analysis scripts
│   ├── molecular_structure_analysis.py
│   ├── pce_sascore_sensitivity_analysis.py
│   └── structural_motif_analysis_enhanced.py
├── notebooks/                          # Jupyter notebooks
│   ├── PCE_GDB9.ipynb                 # Main PCE calculations
│   └── KAMENI2025.ipynb               # Complete workflow
└── figures/                            # Key figures
    ├── pce_sascore_sensitivity_analysis.pdf
    ├── top_molecules_structures.png
    ├── correlation_heatmap_top1000.pdf
    └── donor_acceptor_patterns_top1000.pdf
```

## Installation

### Prerequisites

- Python 3.8 or higher
- Jupyter Notebook
- RDKit

### Setup

1. Clone this repository:
```bash
git clone https://github.com/YOUR_USERNAME/organic-semiconductor-discovery.git
cd organic-semiconductor-discovery
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Install RDKit (recommended via conda):
```bash
conda install -c conda-forge rdkit
```

## Usage

### 1. PCE Calculation and Screening

Run the main PCE calculation notebook:
```bash
cd notebooks
jupyter notebook PCE_GDB9.ipynb
```

This notebook:
- Loads PubChemQC molecular data (HOMO, LUMO, gap)
- Applies Scharber model for PCE prediction
- Calculates SAScore for synthetic accessibility
- Computes PCE_SAScore metric

### 2. Sensitivity Analysis

Validate the PCE_SAScore metric:
```bash
cd scripts
python pce_sascore_sensitivity_analysis.py
```

### 3. Large-Scale Structural Analysis

Analyze 100, 500, or 1000 molecules:
```bash
cd scripts
python structural_motif_analysis_enhanced.py 1000
```

Output: Structure-property correlations, donor/acceptor patterns, statistical validation

### 4. Complete Workflow

For the full analysis pipeline:
```bash
cd notebooks
jupyter notebook KAMENI2025.ipynb
```

## Data Availability

All datasets are available on Zenodo:

**DOI:** [10.5281/zenodo.18201813](https://doi.org/10.5281/zenodo.18201813)

The Zenodo repository includes:
- PCE calculations for 17,458 molecules (`PCE_paper_GDB9.csv`)
- Top 7 candidate molecules (`predictions_molecules_cibles.csv`)
- Large-scale structural analysis of 1,000 molecules (`structural_motif_analysis_top1000.csv`)
- Sensitivity analysis results (`pce_sascore_sensitivity_results.csv`)
- All analysis scripts and notebooks

Original PubChemQC database: [https://pubchemqc.riken.jp/](https://pubchemqc.riken.jp/)

## Key Results

From 17,458 molecules screened:
- **7 candidates** satisfy PCE_SAScore > 0 criterion
- **Top molecule (17851)**: PCE_PCDTBT = 36.11%, SAScore = 7.62
- **40%** show strong fluorescence (f > 0.5)
- **7%** exhibit TADF potential (ΔE_ST < 0.3 eV)
- **Design principles validated** across 1,000 molecules (p < 0.001)

### Structural Insights

**Donor-like molecules (n=356):**
- N-rich (avg 1.13 N atoms)
- Lower O content (avg 1.34 O atoms)
- Moderate aromaticity (avg 0.77 rings)

**Acceptor-like molecules (n=500):**
- O-rich (avg 2.52 O atoms)
- Higher N content (avg 1.48 N atoms)
- Higher aromaticity (avg 0.88 rings)

## Citation

If you use this code or data, please cite:

```bibtex
@article{mvoto2026organic,
  title={Data-Driven Discovery of Synthetically Compatible Organic Semiconductors for Multifunctional Applications: A Computational Workflow},
  author={Mvoto Kongo, Patrick Sorrel and Teguia Kouam, Steve Cabrel and Tchapet Njafa, Jean-Pierre and Nana Engo, Serge Guy},
  journal={Computational Materials Science},
  year={2026},
  publisher={Elsevier}
}

@dataset{mvoto2026data,
  author={Mvoto Kongo, Patrick Sorrel and Teguia Kouam, Steve Cabrel and Tchapet Njafa, Jean-Pierre and Nana Engo, Serge Guy},
  title={Data for: Data-Driven Discovery of Synthetically Compatible Organic Semiconductors for Multifunctional Applications},
  year={2026},
  publisher={Zenodo},
  doi={10.5281/zenodo.18201813},
  url={https://doi.org/10.5281/zenodo.18201813}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

- **Patrick Sorrel MVOTO KONGO** (First Author) - *Original code development, computational analyses, initial draft*
- **Jean-Pierre TCHAPET NJAFA** (Corresponding Author) - *Project supervision, methodology refinement*
- **Steve Cabrel TEGUIA KOUAM** - *Data analysis, manuscript revision*
- **Serge Guy NANA ENGO** - *Data analysis, manuscript revision*

**Corresponding Author:** Jean-Pierre TCHAPET NJAFA  
**Email:** jean-pierre.tchapet@facsciences-uy1.cm  
**Institution:** Department of Physics, Faculty of Science, University of Yaounde 1, Po. Box 812, Yaounde, Cameroon

## Acknowledgments

- PubChemQC database (RIKEN)
- Harvard Clean Energy Project for benchmarking data
- RDKit and AutoDock Vina development teams

## Contributing

This repository contains the code for a published research article. For questions or collaborations, please contact the corresponding author.

## Funding

[Add funding information if applicable]

---

**Repository maintained by:** Jean-Pierre TCHAPET NJAFA  
**Last updated:** 2026-01-09
