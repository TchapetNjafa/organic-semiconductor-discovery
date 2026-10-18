---
name: docking_methods_note
description: Reproducible docking methodology and results for Paper 2 (replaces the unreliable prior runs)
---

# Docking re-run — methodology and results (2026-07-10)

## Why this exists

Two prior docking runs (`MVOTOOPV/PUBCHEMQC_MULTIFONCTIONAL_2025/DATASET/docking_matrix.csv`
and `dock_by_target.csv`) disagreed by 3-4 kcal/mol for identical ligand/target pairs. Root
causes identified: `smina` called unqualified (silently failing in some cells), a hardcoded
`(0,0,0)`/20x20x20 Å box unrelated to any receptor's actual coordinates, an unseeded second
3D-conformer-generation pass (OpenBabel `make3D()` overwriting a seeded RDKit conformer), and
two broken 0-byte receptor conversions (4LDE, 5JWA). Both prior notebooks are considered
unreliable and are **not** used as data sources going forward.

## Scope

- **Ligands (4):** the genuine Paper 1 OPV candidates — 1712, 17851, 20778, 4550 (977, 7801,
  11029 are screening artifacts per Paper 1 and are excluded here).
- **Targets (4):** HIV1protease (PDB 1DMP), Hsp90 (2XJX), a neurodegenerative-disease target
  (1SYH), COVID-19 main protease (6Y2F). Proteasome (7PG9) was excluded — it is the full 28-chain,
  ~45,000-atom 20S core particle (180x146x120 Å), not a tractable single-domain target for blind
  docking. 4LDE and 5JWA were excluded — their PDB->PDBQT conversion failed (0-byte files) and
  were not re-attempted (out of scope for this candidate set).

## Method

Whole-protein **blind docking** (no assumed active site or co-crystallized-ligand pocket) via
smina (`/home/tchapet/SMINA/smina.static`, based on AutoDock Vina 1.1.2).

- **Ligand prep:** RDKit (`ETKDGv3`, `randomSeed=42`) 3D embedding + UFF optimization -> PDB ->
  single-pass conversion to PDBQT via OpenBabel (`--partialcharge gasteiger`, no `--gen3d`, so the
  RDKit conformer is preserved rather than overwritten).
- **Receptor prep:** existing prepared PDBQT files reused as-is (`receptors/{1DMP,2XJX,1SYH,6Y2F}.pdbqt`).
- **Box:** per-receptor heavy-atom bounding box (min/max over all `ATOM`/`HETATM` coordinates in the
  `.pdb`) + 8 Å padding on each side, center = box midpoint. Exact per-target center/size values are
  recorded in `docking_results_clean.csv`.
- **Search parameters:** `exhaustiveness=32`, `num_modes=9`, `seed=42` (fixed for every run).
- **Script:** `redock_candidates.py` (this directory). Verified reproducible: two independent full
  runs produced bit-identical affinities for all 16 ligand/target pairs.

## Results (best-mode affinity, kcal/mol)

| mol_id | HIV1protease | Hsp90 | Neurodegenerative (1SYH) | COVID-19 (6Y2F) |
|--------|:---:|:---:|:---:|:---:|
| 1712   | -6.5 | -7.0 | -6.5 | -7.0 |
| 4550   | -6.4 | -5.8 | -5.8 | -5.9 |
| 17851  | -5.5 | -4.5 | -4.4 | -4.4 |
| 20778  | -4.7 | -4.4 | -4.5 | -4.2 |

**1712 (the sole View-B-viable OPV candidate from Paper 1) is the strongest and most consistent
binder across all four targets.** All values are blind-docking estimates, not validated against
any experimental binding assay — this should be stated plainly in the manuscript as a computational
prediction, not compared to literature inhibitor affinities without a same-protocol reference
ligand run for calibration.

## What is NOT yet done

- No reference/positive-control ligand was docked under the same protocol, so "comparable to known
  inhibitors" cannot be claimed — would need a same-protocol docking of a literature inhibitor for
  each target to calibrate the numbers.
- No flexible side chains / induced fit; rigid-receptor blind docking only.
- No replicate runs with different seeds to estimate result variance (the reproducibility check
  only confirms *determinism*, not that the funnel search actually converges on the true global
  minimum for these box sizes).
