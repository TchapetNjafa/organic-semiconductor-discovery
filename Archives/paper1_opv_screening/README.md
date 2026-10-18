# Paper 1 — chemical-validity filtering in high-throughput OPV screening

Screening pipeline, inputs and regenerated artifacts for the Digital Discovery
manuscript *Efficiency and Accessibility Are Not Enough: Chemical-Validity
Filtering in High-Throughput Screening of Organic Photovoltaic Materials*.

Everything here runs standalone. `bash run_all.sh` regenerates every number,
figure and table in the manuscript from the deposited inputs.

## What the pipeline does

17458 molecules from the first identifier shard of the PubChemQC
`CHNOPSFClNaKMgCa500` subset (B3LYP/6-31G*//PM6) are passed through:

1. **Photovoltaic stage** (`scharber_recompute.py`) — Scharber construction with
   the tabulated ASTM G173-03 spectrum. A molecule–reference pairing is admitted
   only if it satisfies four conditions:
   G1 usable absorption (1.1–2.2 eV under the
   complementary convention, narrower-gap absorber under the strict one),
   G2 E_CT ≤ E_g (charge transfer not uphill),
   G3 V_OC > 0, G4 LUMO–LUMO offset ≥ 0.3 eV.
2. **Chemical-validity filter** (`filter_genuine_osc.py`) — rules R1–R6:
   single covalent species, metal-free, at least one aromatic ring, ≥ 6 heavy
   atoms, ≥ 6 conjugated atoms. 10215 of 17458
   molecules pass.
3. **Reactive-group screen** (`stability_screen.py`) — rule R7, three SMARTS
   motifs (azide, nitroso, diazo). 10162 molecules pass both.
   **This is an inspection diagnostic, not a safety filter** — the script prints
   the motifs it does *not* test for.
4. **Artifacts** (`regen_tables.py`, `regen_figures.py`, `regen_fig7.py`) — the
   figures and tables of the manuscript, from `screen_full.csv`.

## Headline results, as deposited

| quantity | value |
|---|---|
| molecules screened | 17458 |
| chemically valid (R1–R6) | 10215 |
| valid and stable (R1–R7) | 10162 |
| spectrum integral | 1000.371 W m⁻² |
| J_SC at E_g = 1.4 eV | 213.73 A m⁻² |
| max PCE, FF = 0.65 | 13.672 % |
| max PCE, Green FF (n = 1) | 18.877 % |
| admitted → valid → above threshold (complementary, FF 0.65) | 13 → 6 → 4 |
| admitted → valid → above threshold (strict, FF 0.65) | 8 → 2 → 1 |
| artifacts inside the gap window | 3 |

**First place is convention-dependent** — this is a result, not a nuisance:

| absorber convention | fill factor | first place |
|---|---|---|
| strict | 0.65 | 9168 |
| strict | Green | 9168 |
| complementary | 0.65 | 19598 |
| complementary | Green | 22936 |

Intersection of all four: **9168**
(C₂₆H₁₆, an unsubstituted polycyclic aromatic hydrocarbon). Union:
9168, 19598, 21736, 22936.

Molecule 1712 is **not** admitted: LUMO -3.823 eV gives -0.123 eV driving force vs PCBM (-3.70) and +0.223 eV vs PCDTBT (-3.60); both below the 0.3 eV requirement

## Reproducing

```bash
python -m pip install -r requirements.txt
bash run_all.sh
```

Provenance of the input shard is fixed by SHA-256, not by a retrieval date:
`6f02d7b8d39bfc97885a52ca2b55f6599e09a27b601bf5da092d94ed4b517a69` (see `INPUT_CHECKSUM.txt` for the archive name, byte count
and entry count).

## Layout

```
paper1_opv_screening/
  run_all.sh                     regenerate everything
  scharber_recompute.py          photovoltaic stage, gates G1-G4
  filter_genuine_osc.py          chemical-validity filter R1-R6
  stability_screen.py            reactive-group screen R7
  regen_tables.py                LaTeX tables
  regen_figures.py  regen_fig6.py  regen_fig7.py  figures
  figstyle.py                    shared plot style
  requirements.txt               pinned versions
  INPUT_CHECKSUM.txt             SHA-256 of the PubChemQC shard archive
  data/
    dataset_pubchemqc_opv_17458.csv    the property table (Zenodo only)
    ASTMG173.csv                       reference solar spectrum
    reorganization_energies_corrected.csv
  outputs/recompute/
    screen_full.csv              per-molecule results, all four conventions
    summary.json                 the numbers in the table above
    robustness.json              convention comparison, intersection, 1712
    recompute_log.txt            full run log with every intermediate check
  genuine_osc_ranked.csv         validity-filter output
  viable_view_{a,b}*.csv         before / after the reactive-group screen
  figures/                       manuscript figures (PDF, and PNG on Zenodo)
  tables/                        manuscript tables (LaTeX)
```

## Scope and limitations

The Scharber construction is a thermodynamic upper bound from frontier-orbital
energies. It has no morphology, no mobility, no exciton diffusion and no
recombination, and the efficiencies here rank molecules — they do not predict
device performance. B3LYP/6-31G*//PM6 Kohn–Sham gaps are not optical gaps: the
nine established organic semiconductors present in this shard (benzene through
pentacene) are **all rejected** by G1 for that reason, which bounds how much
weight the admitted set can carry. See `tables/tab_comparators.tex`.

## Citation

If you use this deposit, cite the manuscript and the Zenodo record
(concept DOI 10.5281/zenodo.18201812, resolves to the latest version). Mirror:
<https://github.com/TchapetNjafa/organic-semiconductor-discovery>.

Licence: see `LICENSE` in the repository root (code) — the PubChemQC source data
retains its own terms and is not redistributed here beyond the derived property
table.
