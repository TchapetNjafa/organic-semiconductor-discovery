#!/usr/bin/env bash
# Regenerate every Paper-1 number, figure and table from the deposited inputs.
# Needs the pinned versions in requirements.txt:  pip install -r requirements.txt
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
echo "[1/7] screening stage: gates G1-G4, both absorber and both FF conventions"
python scharber_recompute.py
echo "[2/7] chemical-validity filter R1-R6 (cross-checks the pipeline's own rules)"
python filter_genuine_osc.py
echo "[3/7] reactive-group screen R7 (cross-checks the pipeline's own SMARTS)"
python stability_screen.py
echo "[4/7] LaTeX tables"
python regen_tables.py
echo "[5/7] figures 1,2,3,8,9"
python regen_figures.py
echo "[6/7] figure 6 (heteroatom correlations)"
python regen_fig6.py
echo "[7/7] figure 7 (reorganization energies)"
python regen_fig7.py
echo
echo "done. outputs/recompute/ holds the numbers, figures/ and tables/ the artifacts."
