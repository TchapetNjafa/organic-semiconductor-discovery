# Paper 1 (Adv. Theory Simul.) — Analysis Ledger

**Cap: ≤ 400 lines / ~5,000 tok.** Live entries only; supersede, never overwrite. IDs `L-YYYY-MM-DD-NN`.
Created 2026-10-17. Entries with IDs earlier than 2026-10-17 (e.g. `L-2026-09-02-23`) belong to an earlier project
ledger that is not part of this deposit.
Primary data: `ZENODO_UPLOAD/dataset_pubchemqc_opv_17458.csv` (never read raw) · sTDA-xTB spectra
`NTO_and_related_study/data/xtb_results/<id>/tda.dat` · extracts in `outputs/extract/`.
Environment: venv `~/VirtualEnv`; `~/xtb4stda/xtb4stda`, `~/stda_v1-6/stda`; Lean via physlib Mathlib (v4.34.1).

## L-2026-10-17-01 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (sTDA-xTB S1 of the nine reference PAHs in the shard (calibration input)). Verified against: see archived entry.

## L-2026-10-17-02 — Gap observable decides the survivor set: KS gap vs sTDA-xTB S1, full shard

**Method.** Join of `dataset_pubchemqc_opv_17458.csv` (KS GAP) with `outputs/extract/stda_all.csv` on mol_id
(15,952 molecules with an sTDA spectrum; 1,506 without). Window [1.1, 2.2] eV applied to each observable alone,
no other gate. Script inline (pandas 3.0.2) → `outputs/extract/ks-vs-stda-summary.txt`, `ks_vs_stda_joined.csv`.
**Result.** Pearson r(KS, S1) = 0.677, Spearman 0.802. S1 − KS: median −1.058 eV, IQR [−1.527, −0.688],
5–95 % [−3.006, −0.168], mean −1.237 ± 0.993 (SD). In window: **KS 33 · S1 559 · bright state 407 ·
both 4**. Of the 33 KS-window molecules, 9 have S1 > KS + 1 eV (10657 C15H24, 16955, 19598, 9168, 21736,
977 O₂, 23450, 4988, 9104) and 20 have S1 below KS, 15 of them with S1 < 1 eV (salts, ion pairs, e.g. 1712).
**Interpretation.** Swapping the gap observable replaces the survivor set almost entirely (4 of 33 shared),
so the choice of observable, not the ranking metric, determines what a Scharber screen returns — the
central hypothesis holds at the window step. The KS window is enriched in two artifact classes visible
only when the two observables are compared: (i) records whose KS gap is anomalously small relative to S1
(a saturated C15H24 at 1.66 eV, triplet O₂, *o*-quinone tautomers, 9168), and (ii) ionic species where
sTDA-xTB collapses to near-zero CT states. Class (ii) means sTDA-xTB has its own failure mode, so the
validity gate must precede the optical gate.
**Verified against.** `dataset_pubchemqc_opv_17458.csv` + `stda_all.csv` (built from raw `tda.dat`); independently recomputed by copilot (-13): exact match.
**Caveats.** 1,506 molecules lack spectra (failed sTDA); bias of that subset unknown. sTDA on xTB geometry.
559 is before any validity gate and contains unknown artifacts. Correlations are descriptive, n = 15,952.
**Claim link.** C-central (hypothesis) · C-order · C-artifacts.

## L-2026-10-17-03 — 9168's PubChemQC Kohn–Sham gap does not reproduce: the apparent lead is a record artifact

**Method.** RKS B3LYP/6-31G* (PySCF 2.13.0, libxc B3LYP id 402) single points on the xTB geometries
(`NTO_and_related_study/data/xtb_results/<id>/mol.xyz`) for 9168 and two controls of the same class,
pentacene 8671 and chrysene 9171; UKS triplet probe at the same geometry. Script
`Paper1_AdvTheorySimul/ks_gap_check.py` → `outputs/extract/ks-gap-recheck-9168.txt`. All SCFs converged.
**Result.** Gap recomputed vs PubChemQC [eV]: pentacene **1.896 vs 2.397 (−0.501)**, chrysene **4.018 vs 4.365
(−0.347)**, **9168 3.784 vs 1.771 (+2.013)**. 9168 HOMO/LUMO −5.258/−1.474 (PubChemQC −4.980/−3.208).
E(T1)−E(S0): pentacene 0.653, chrysene 2.699, 9168 2.676 eV — 9168 is a closed-shell singlet with a
chrysene-like triplet energy. sTDA-xTB S1 of 9168 (3.441 eV, `-01`) sits 0.34 eV below the recomputed gap, the
normal ordering.
**Interpretation.** The two controls shift down by 0.35–0.50 eV when the geometry changes from the PubChemQC one
to xTB; 9168 shifts **up by 2.0 eV**, the opposite sign and four times the size, while its LUMO moves by 1.7 eV.
A geometry effect cannot do that, so the 1.771 eV value in the dataset belongs to a defective record (wrong
structure or state), not to the molecule its SMILES describes. 9168 therefore never belonged in the 1.1–2.2 eV
window: the molecule that an uncalibrated screen keeps as "the only candidate robust to all four conventions"
was admitted by a database error. The same comparison flags the 31 molecules with S1 > KS + 1 eV (`-02`)
as candidates for the same defect, which turns the KS/S1 comparison into a record-validity test.
**Verified against.** fresh PySCF single points on the geometries (raw output); challenged by two same-class controls and by B3LYP re-optimisation (-09).
**Caveats.** Geometry is xTB, not the PubChemQC B3LYP one (the PubChemQC geometry is not on disk); the control
offset (−0.35 to −0.50 eV) bounds the geometry effect only for this class. Pentacene's 1.896 eV is lower than
typical B3LYP/6-31G* literature values [CITATION NEEDED]. n = 3; no uncertainty model. Re-optimising 9168 at
B3LYP/6-31G* would close the geometry caveat (cost ≈ hours) — do before the claim reaches the abstract.
**Claim link.** C-9168 (not a lead) · C-artifacts · supersedes the "robust lead" reading of `L-2026-09-02-19`.

## L-2026-10-17-04 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (CORRECTION: sTDA/PySCF geometries are RDKit ETKDG + UFF (seed 42), not xTB). Verified against: see archived entry.

## L-2026-10-17-05 — Halogens at the admission frontier are chloride counter-ions, not covalent halogens

**Method.** `ks_vs_stda_joined.csv` (n = 15,952) + dataset SMILES. Halogen flags from formula; multi-fragment =
SMILES containing "."; Cl⁻ counter-ion = `[Cl-]` or bare `Cl` fragment; covalent Cl = Cl in formula and not a
counter-ion. Fisher exact test of halogen enrichment per window vs the rest. → `outputs/extract/halogen-frontier.txt`.
**Result.** Shard: 21.7 % halogenated (Cl 2,753, F 854, **no Br or I**); 4.5 % multi-fragment; 2.9 % Cl⁻
counter-ion; 14.4 % covalent Cl. Halogen enrichment (odds ratio, p): KS window 36.4 % (2.06, 0.055); S1 window
22.2 % (1.03, 0.79); bright window 37.3 % (2.20, 5e-13); S1 < 1 eV 61.5 % (6.46, 2e-127). Multi-fragment share:
**KS window 48.5 %**, S1 window 10.0 %, bright window 39.8 %, **S1 < 1 eV 84.4 %**. Covalent Cl: KS window
**0.0 %**, S1 window 13.2 %, shard 14.4 %. Single-fragment molecules only: KS window 17, S1 window 503.
**Interpretation.** The halogen excess near the frontier is carried entirely by salts: chloride counter-ions
make up a third of the KS window and over half of the sTDA collapse class, while covalently bound Cl and F —
the motifs used to tune NFA energy levels — occur in the S1 window at their shard frequency. So the answer to
"do halogenated species cluster at the frontier" is no for covalent halogens and yes for ion pairs, which are
artifacts for both observables. Removing multi-fragment records is the single most effective validity rule here:
it halves the KS window and removes 84 % of the sTDA collapse cases.
**Verified against.** dataset SMILES/formula columns; challenged by multiple-comparison correction (-12).
**Caveats.** Fragment flag is syntactic (SMILES "."); no Br/I in shard so heavy-halogen effects untestable here;
the bright-window 1e-13 enrichment is a multiple-comparison-safe signal but its cause is the same salt class.
**Claim link.** C-validity-order · C-halogens · C-artifacts.

## L-2026-10-17-06 — HOPV15: raw sTDA-xTB S1 overestimates OPV optical gaps by 0.75–1.0 eV; after calibration all observables tie

**Method.** HOPV15 v4 (Lopez 2016, CC BY 4.0, md5 verified; 350 records, 343 sTDA OK, 7 ETKDG embed failures).
Shard protocol `Paper1_AdvTheorySimul/hopv_stda.py` (SMILES→ETKDG/UFF→xtb4stda→stda) and the same sTDA on HOPV's
own DFT conformer-1 geometry (`hopv_stda_dftgeom.py`). KS = HOPV B3LYP/def2-SVP conformer 1. Reference =
experimental optical gap column. Classes: small molecule (complement ≠ TiO₂), polymer (oligomer fragment), DSSC.
Leave-one-out (LOO) linear calibration. → `outputs/extract/hopv15-calibration.txt`, `outputs/extract/hopv15/`.
**Result.** OPV small molecules with an experimental gap, **n = 34** (exp 1.49–2.58 eV):
KS bias **+0.293**, MAE 0.293, r 0.922 · S1 (UFF) bias **+0.999**, r 0.907 · S1 (DFT geometry) bias **+0.754**,
r 0.894; UFF − DFT-geometry S1 = +0.245 ± 0.183 eV. Window [1.1, 2.2] eV recall of the 29 molecules whose
experimental gap lies in it: **KS 21/29, S1(UFF) 1/29, S1(DFT geom) 0/29**. LOO-calibrated MAE: KS 0.086,
S1(UFF) 0.097, S1(DFT geom) 0.105 eV (RMSE 0.115–0.129). Polymers (n = 174): r ≤ 0.26 for every observable
(fragments ≠ polymer); DSSC n = 3, not interpretable.
**Interpretation.** On real donor–acceptor OPV molecules, an *uncalibrated* sTDA-xTB S1 is not an optical gap:
it sits ~1 eV too high, so a raw optical-gap window rejects almost every working OPV molecule — the failure
commonly attributed to the KS gap occurs here for sTDA, worse. A quarter of the error is the UFF geometry;
the rest is the method on push–pull chromophores. Once each observable is linearly calibrated, the three tie
at ≈0.1 eV: the gap step's error comes from calibration, not from which observable is chosen.
**Verified against.** HOPV15 v4 file (md5-checked) + raw sTDA outputs; independently recomputed by copilot (-13): match, n definition clarified.
**Caveats.** n = 34; in-class LOO only (no out-of-class test here); KS basis differs from PubChemQC (def2-SVP vs
6-31G*); experimental gaps are literature values of mixed method (film/solution onset); polymers excluded.
**Claim link.** C-central (hypothesis revised) · C-validation vs experiment · C-calib.

## L-2026-10-17-07 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (The raw bias is class-dependent for sTDA-xTB (0 vs +0.75–1.0 eV), less so for the KS gap). Verified against: see archived entry.

## L-2026-10-17-08 — Gap step on the full shard: validity gate, raw vs HOPV-calibrated KS, raw S1

**Method.** All 17,458 records; validity = `filter_genuine_osc.reject_reason` (R1–R6, copied to
`Paper1_AdvTheorySimul/pipeline/`). KS calibration fitted on HOPV15 OPV molecules (n = 34, `-06`):
E_opt = 0.681·KS + 0.391. Window [1.1, 2.2] eV on each observable. → `outputs/extract/gap-step-screen.txt`,
`gap_step_screen.csv`.
**Result.** Validity rejects 7,243 (no ring 4,034; no aromatic ring 2,289; multi-fragment 838; low conjugation 70;
<6 heavy atoms 9; metal 2; unparsable 1); valid 10,215. Window-in (all / valid): KS raw 41 / **13**; KS calibrated
98 / **46**; S1 raw 559 / 391. Valid overlaps: KS raw ⊂ KS calibrated (13 shared, 33 new); KS raw vs S1 raw
share 3; KS calibrated vs S1 raw share 28. Calibrated KS: pentacene **2.023** (in), tetracene 2.387, anthracene
2.917 (out); 9168 1.597 (in — still a defective record, `-03`).
**Interpretation.** Calibration, not the observable, moves the survivor set most usefully: it admits pentacene,
which the raw KS window rejected, and triples the valid survivors (13 → 46) while keeping every raw survivor.
Raw S1 admits 391 valid molecules but, by `-06`, it is ~1 eV high for D–A chromophores and ≈0 for acenes, so its
survivors are a class-biased set. Calibration does not remove defective records: 9168 stays in, so a
record-validity check is needed in addition to the validity gate.
**Verified against.** dataset + validity filter code; challenged by calibration-noise sensitivity (-12). Not yet independently recomputed by a second agent.
**Caveats.** Calibration fitted on 34 D–A molecules (def2-SVP KS) and applied to 6-31G* KS of a different chemical
space — extrapolation; window edges are conventions; no Scharber offset gate applied yet (gap step only).
**Claim link.** C-central (v2 findings i, iii) · C-calib · C-order.

## L-2026-10-17-09 — The S1-vs-KS flag finds defective PubChemQC records: 21/29 flagged vs 0/29 controls; 9168 closed

**Method.** Flag = single-fragment record with S1(sTDA-xTB) − KS(PubChemQC) > 1 eV (29 records). Controls = 29
random single-fragment records with S1 − KS in [−1.5, 0.5] eV (seed 7). B3LYP/6-31G* RKS single point on the
UFF geometry (`ks_gap_check.py`, TRIPLET=0) for all 58; defect = recomputed gap exceeds the PubChemQC gap by
> 0.5 eV. Separately, 9168 and pentacene re-optimised at B3LYP/6-31G* (geomeTRIC, `b3lyp_reopt.py`).
→ `outputs/extract/record-validity-{sample.csv,recompute.txt,result.csv,summary.txt}`, `b3lyp-reopt/`.
**Result.** All 58 SCFs converged. Recomputed − PubChemQC: flagged mean **+1.309 eV** (median |dev| 1.521),
**21/29 > +0.5 eV, 18/29 > +1.0 eV**; controls mean −0.263 (median |dev| 0.326), **0/29 > +0.5 eV** (one-sided
Fisher p from the summary file; flag precision 0.72, Clopper–Pearson 95 % CI in file; control defect rate
< 0.12 at 95 %). Largest defects: 10657 C15H24 1.663 → 6.580; 17600 C20H30O2 3.684 → 6.587; several steroid-like
C23–C25 oxygenates 3.8 → 5.3. Two flagged *o*-quinone-type records (19598, 21736) move the other way (−0.8 eV):
the SMILES-built tautomer differs from the PubChemQC one. **B3LYP-optimised 9168: gap 4.021 eV** (PubChemQC
1.771); pentacene 2.209 eV (PubChemQC 2.397).
**Interpretation.** A KS gap far below the cheap S1 identifies records whose stored KS gap belongs to something
other than the molecule in the SMILES: 72 % of flagged records are wrong by > 0.5 eV, while unflagged records show
only the −0.3 eV geometry shift seen in `-03`. The test costs one sTDA-xTB run per molecule. It is a filter for
*suspects*, not a verdict: confirmation needs a recompute, and tautomer mismatches give false alarms of the
opposite sign. 9168's record defect is now closed without the geometry caveat (2.25 eV too low).
**Verified against.** raw PySCF outputs for all 58 records; independently recomputed by copilot (-13): 21/29 and 0/29 match.
**Caveats.** Recompute uses UFF geometries (controls bias −0.26 eV), single-point only except 9168/pentacene;
flag recall unknown (only flagged and in-range controls sampled); 31 flags in 15,952 = 0.2 % of records.
**Claim link.** C-central (v2 finding ii) · C-9168 (closed) · C-artifacts. Supersedes caveat of `-03`.

## L-2026-10-18-10 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (Lean F1-F3 machine-checked, FAITHFUL). Verified against: see archived entry.

## L-2026-10-18-11 — Full calibrated cascade (donor vs PCBM): 37 survivors, pentacene admitted, 9168 removed

**Method.** `Paper1_AdvTheorySimul/gate_cascade.py`, reusing `pipeline/scharber_recompute.py` (AM1.5G J_SC,
EQE 0.65, V_loss 0.3 V, PCBM −6.10/−3.70 eV, Green FF n = 1, 300 K). Order: G0 validity R1–R6 → G1 record flag
(drop S1 − KS > 1 eV, `-09`) → G2 window [1.1, 2.2] eV on the HOPV-calibrated optical gap (E_opt = 0.681·KS + 0.391,
`-08`) → G3 V_OC > 0 → G4 LUMO offset ≥ 0.3 eV. V_OC and offset use stored KS frontier energies; J_SC uses E_opt.
Strict convention only (molecule sets the onset). → `outputs/extract/gate_cascade.{txt,csv}`.
**Result.** Cumulative survivors: G0 10,215 → G1 10,208 → G2 40 → G3 40 → G4 **37**. PCE (FF 0.65): max
**10.14 %**, median 5.52 %; Green FF max 14.32 %. Pentacene **admitted** (E_opt 2.023 eV, PCE 4.06 %); tetracene
out (2.387 eV); **9168 removed by G1**; 19598 (*o*-quinone) removed by G1; anthraquinone dye 22936 still admitted
(PCE 9.37 %). 7/37 survivors have S1 more than 0.5 eV below E_opt (one, 16356, at 0.369 eV — the sTDA
collapse class); all 37 have an sTDA spectrum.
**Interpretation.** With calibration and the record flag in place the screen admits the one acene whose optical
gap lies in the window and drops the database artifact the uncalibrated screen ranked first. The survivor count
is small and the best PCE (10 %) is about half the uncalibrated maximum, so most of the earlier headroom came
from uncalibrated gaps. The cascade still keeps a dye that a chemist would reject (22936) and molecules whose two
observables disagree downward; the record flag is one-sided by design and catches only one failure direction.
**Verified against.** dataset + `stda_all.csv` + HOPV calibration via `gate_cascade.py`; reference molecules checked by hand. **Not yet independently recomputed** — to do before any number reaches the manuscript.
**Caveats.** Calibration extrapolated from 34 D–A molecules; KS frontier energies uncalibrated for V_OC; no
reactive-group (R7) or stability screen applied here; PCE is a Scharber upper bound, not a prediction; flag
threshold 1 eV can remove genuine D–A chromophores (HOPV S1−KS up to +1.31 eV, `-06` data).
**Claim link.** C4 (extends), C6, new C11 (calibrated screen outcome).

## L-2026-10-18-12 — Statistical hardening after the kiro-cli overclaim pass (amends -02, -05, -06, -08, -09)

**Method.** Independent overclaim review by kiro-cli (effort high) → `outputs/cli-agent/L1-kiro-overclaim.md`
(verdict OVERCLAIMS on -02, -05, -06, -07, -08). Fixes computed here → `outputs/extract/stat-hardening.txt`:
Holm correction of the four Fisher tests of -05; paired bootstrap (10,000 resamples, seed 2026) of LOO absolute
errors of -06; window-edge shift by ±1 SD of the LOO calibration error (0.097 / 0.115 eV) for -08.
**Result.** (-05) Holm-adjusted p: KS window 0.11 (not significant), S1 window 0.79, bright window 1.5e-12,
S1 < 1 eV 7.1e-127. (-06) LOO MAE [95 % CI]: KS 0.086 [0.061, 0.113]; S1 UFF 0.097 [0.074, 0.122]; S1 DFT-geom
0.105 [0.080, 0.132] eV; paired differences vs KS +0.011 [−0.021, +0.041] and +0.019 [−0.004, +0.041] eV.
(-08) valid calibrated-KS survivors 46 at the nominal edges, **34–58** for ±1 SD calibration shifts.
(-09) one-sided Fisher p = 1.28e-9 (inline).
**Interpretation.** Amended readings, which replace the stronger wording of the earlier entries:
-02's "the central hypothesis holds" → the window counts are *consistent with* an observable-controlled survivor
set (they are a census of the shard, so no sampling CI applies, but one library is one library); -05's halogen
excess in the KS window is a non-significant trend after correction, while the salt/collapse enrichment survives;
-06's "the three tie" is supported — no paired difference excludes zero at n = 34; -08's 46 survivors carry a
calibration-noise range of 34–58, so the survivor count is known to about ±25 %, not to the unit.
-07's class-dependence claim awaits the def2-SVP acene check (running, `outputs/extract/acene-def2svp.txt`).
**Verified against.** recomputed from the data CSVs (bootstrap, Holm), not from the ledger values.
**Caveats.** Bootstrap at n = 34 is coarse; window shift assumes a uniform calibration error.
**Claim link.** C1, C3, C4, C7, C5 (wording) — manuscript must use the amended phrasing.

## L-2026-10-18-13 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (copilot independent recompute of -02/-06/-09: MATCH; n = 36 vs 34 = 2 ETKDG failures). Verified against: see archived entry.

## L-2026-10-18-14 — DFT ground state (std2) removes sTDA's D–A bias but adds scatter; calibrated KS stays best (C10)

**Method.** PySCF 2.13 RKS B3LYP/6-31G* (Cartesian) on HOPV15 DFT conformer-1 geometries → Molden (`[Atoms] Angs`)
→ std2 2.0.1 sTDA and XsTDA (`-ax 0.2 -e 7`), runner `Paper1_AdvTheorySimul/std2_b3lyp.py`; 36 OPV small molecules
with an experimental optical gap (DFT geometries need no ETKDG, hence 36 not 34). Chain validated on pentacene vs full
TDA-B3LYP (`outputs/literature/std2-assessment.md`). LOO linear calibration; paired bootstrap (10,000, seed 2026)
on the 34 complete cases. → `outputs/extract/hopv15-std2-vs-exp.txt`, `hopv15/std2_reparsed.csv`.
**Verified against.** S1 re-extracted from every raw `std2/*/{stda,xstda}.out` with the corrected parser (header
required): 36/36 parsed, 0 mismatches with the run's summary file. Experimental gaps from `hopv15_records.csv`.
**Result.** vs experiment [eV] (bias / MAE / SD / r / raw-window recall / LOO MAE):
KS B3LYP/def2-SVP +0.288 / 0.288 / 0.155 / 0.921 / 23/31 / 0.083 · **sTDA@B3LYP −0.021 / 0.211 / 0.256 / 0.747 /
28/31 / 0.141** · XsTDA@B3LYP +0.237 / 0.247 / 0.289 / 0.692 / 23/31 / 0.157 · sTDA-xTB (DFT geom) +0.758 / 0.758 /
0.121 / 0.891 / 0/31 / 0.103 · sTDA-xTB (UFF, n 34) +0.999 / … / 1/29 / 0.097. Paired vs KS (n 34): raw MAE
sTDA@B3LYP −0.082 [−0.138, −0.030], XsTDA −0.035 [−0.098, +0.037]; LOO MAE sTDA@B3LYP **+0.060 [+0.031, +0.091]**,
XsTDA +0.076 [+0.046, +0.109], sTDA-xTB +0.011–0.019 (CIs include 0).
**Interpretation.** Moving the excited-state calculation onto a B3LYP ground state removes the ~0.75 eV donor–acceptor
bias of sTDA-xTB — so that bias comes from the tight-binding ground state, not from the sTDA step — and makes sTDA
the most accurate *uncalibrated* gap observable here (28/31 recalled by the raw window). It does not make it the best
screening observable: its molecule-to-molecule scatter is larger than the KS gap's, so once each observable gets a
linear calibration the stored KS gap is significantly better (by 0.06 eV), and it is already in the database at no
cost. The exact-integral XsTDA variant brings no gain at this level. For a screen, cost and calibrated error point
to the KS gap plus calibration; std2 is useful as an uncalibrated cross-check tier.
**Caveats.** n = 36, one conformer, mixed experimental methods; basis differs (6-31G* std2 vs def2-SVP KS);
S1 is the lowest state regardless of brightness; std2 shares B3LYP with the KS gap, so the two are not independent.
**Claim link.** C10 (answered), C3.

## L-2026-10-18-15 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (basis vs geometry vs class in the KS offset). Verified against: see archived entry.

## L-2026-10-18-16 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (codex independent re-implementation of the -11 cascade: MATCH (12919 PCE +0.0014 pt)). Verified against: see archived entry.

## L-2026-10-18-17 — Applicability domain: the calibration interpolates in gap but extrapolates in chemistry (death-angle 2)

**Method.** `Paper1_AdvTheorySimul/applicability.py`: six RDKit descriptors (heavy atoms, conjugated fraction, aromatic
rings, heteroatom fraction, S count, rotatable bonds); Mahalanobis distance to the 34 HOPV calibration molecules
(standardised); "inside" = below the 95th percentile of HOPV leave-one-out distances (4.29). Sets: HOPV calibration,
37 cascade survivors (-11), 2,000 random single-fragment shard records (seed 1). → `outputs/extract/applicability.{txt,csv}`.
**Verified against.** SMILES columns of `hopv15_records.csv` and `gate_cascade.csv` (primary data); counts re-read
from `applicability.csv`. Independently re-implemented by codex (effort high, detached): **MATCH** on every quantity
(`outputs/cli-agent/L1-codex-domain-recompute.md`); the random-shard row was outside that check.
**Result.** Median [p5, p95]: heavy atoms HOPV 54.5 [24.7, 68.7] vs survivors 24 [10, 33.2] vs shard 15 [6, 29];
aromatic rings 8 / 2 / 1; S atoms 5.5 / 0 / 0. Inside the HOPV domain: **survivors 27 % (10/37)**, shard 18 %.
KS range of the calibration fit [1.605, 3.347] eV; survivors' KS [1.736, 2.629] eV — **0/37 outside the fit range**.
**Interpretation.** The calibration is applied inside the gap range it was fitted on, so the problem is not
numerical extrapolation; it is chemical: three survivors in four are smaller, sulfur-free and less ring-rich than
any typical HOPV donor–acceptor molecule, i.e. outside the class the calibration describes. Combined with -18, their
calibrated gaps carry an error that the in-class LOO MAE (≈ 0.09 eV) does not represent.
**Caveats.** Descriptor choice and the 95 % threshold are conventions; Mahalanobis on n = 34 in 6 dimensions is coarse.
**Claim link.** new C12 (domain), C11 (qualified).

## L-2026-10-18-18 — Out-of-sample test: calibration transfers across papers within D–A, fails across classes (death-angle 7)

**Method.** (i) Leave-one-paper-out (10 source DOIs, largest group 12) within the 34 HOPV molecules, linear calibration
of KS, sTDA-xTB (UFF) and sTDA@B3LYP (std2). (ii) External class: shard acenes (PubChemQC KS 6-31G*, sTDA-xTB S1) vs
experimental bright La midpoints (Malloci 2007), KS calibrated with the HOPV fit (0.681·KS + 0.391).
→ `outputs/extract/out-of-sample.txt`.
**Verified against.** HOPV CSVs, `std2_reparsed.csv`, `ks_vs_stda_joined.csv` (primary data); Malloci 2007 tables. Independently recomputed by copilot (gpt-5.4, effort high, detached run): **MATCH** on every stated quantity
(`outputs/cli-agent/L1-copilot-oos-recompute.md`).
**Result.** Leave-one-paper-out MAE (max): KS **0.102** (0.331), sTDA-xTB 0.105 (0.307), sTDA@B3LYP 0.168 (0.550) eV;
LOO-molecule 0.086 / 0.097 / 0.146. Acenes (n = 4): KS raw MAE 0.315 (bias +0.315); **KS calibrated 0.406 (bias −0.406)**,
error growing with gap (pentacene −0.19 → naphthalene −0.72); sTDA-xTB raw **0.125** (bias −0.118).
Pentacene experimental La 2.12–2.31 eV straddles the 2.2 eV window edge.
**Interpretation.** Within the donor–acceptor class the calibration survives the removal of whole source papers (MAE
rises only from 0.086 to 0.102 eV), so it is not fitted to one lab's molecules. Across classes it fails: applied to
acenes it makes the KS gap *worse* than no calibration, because its slope (0.68) compresses gaps that, for acenes,
need almost no correction. This is the central claim in its sharpest form — no single calibration of a cheap gap
observable transfers between chemical classes — and it means a screen must calibrate per class and flag molecules
outside every calibrated class (-17). "Calibration admits pentacene" (C4) is not evidence of a better screen:
pentacene's measured gap sits on the window edge.
**Caveats.** n = 4 external molecules, one class; experimental values are a range across media; the acene KS values
come from PubChemQC geometries (6-31G*), the fit from HOPV (def2-SVP; basis effect 0.01 eV, -15).
**Claim link.** C2 (strengthened), C3 (qualified: in-class), C4 (weakened), new C12.

## L-2026-10-18-19 — Table values (T1 complete-case statistics, T2 cascade) generated from the data

**Method.** `Paper1_AdvTheorySimul/make_tables.py`: T1 on the 34 HOPV15 molecules where all four observables exist
(same molecules in every column); LOO and leave-one-paper-out linear calibration; bootstrap 10,000, seed 2026. T2
recomputes the cascade gates from `gate_cascade.csv` (asserted equal to `keep`). → `manuscript/tables/T{1,2}.tex`,
`tables_values.json`.
**Verified against.** HOPV15 CSVs, `std2_reparsed.csv`, `ks_vs_stda_joined.csv`, `gate_cascade.csv`, `applicability.csv`.
Values that also appear in -06, -12, -18 match to the printed digit.
**Result.** New on the complete-case set: sTDA@B3LYP bias −0.010 eV, MAE 0.211 eV, r 0.748, raw-window recall 26/29,
LOO MAE 0.146 [0.107, 0.191] eV; sTDA-xTB (DFT geom) leave-one-paper-out MAE 0.122 eV. Unchanged: KS bias +0.293,
LOO 0.086 [0.061, 0.113], LOPO 0.102; sTDA-xTB (UFF) +0.999, LOO 0.097 [0.074, 0.122], LOPO 0.105; acenes KS 0.32 raw /
0.41 calibrated, sTDA-xTB 0.12. T2 as in -11/-16/-17. New in T2: survivors' PCE range 1.86–10.14 % (FF 0.65; median
5.52 %) and 1.71–14.32 % with Green's FF (median 7.62 %).
**Interpretation.** Restricting to the same 34 molecules leaves every conclusion of -14 unchanged (std2 nearly unbiased,
best raw recall, worst calibrated error); the n = 36 figures in -14 and the n = 34 figures here differ only in the third
decimal, and the manuscript tables use the complete-case set throughout so that columns compare the same molecules.
**Caveats.** As in -06, -14, -18.
**Claim link.** C2, C3, C10, C11, C12 (table values).

## L-2026-10-18-20 — CORRECTION: stored gaps are B3LYP/6-31G* single points on PM6 geometries (PubChemQC B3LYP//PM6)

**Method.** Re-read the dataset provenance in the earlier project ledger (`L-2026-09-01-01`): the 17,458 records are the
`.B3LYP@PM6.S0.json` entries of PubChemQC block `CHNOPSFClNaKMgCa500`, IDs 1–24,999 (7,541 IDs without a record),
deterministic, no sampling, no element restriction.
**Verified against.** Earlier ledger entry and its provenance file (grep); dataset documentation: Nakata & Maeda,
J. Chem. Inf. Model. 2023 (10.1021/acs.jcim.3c00899, B3LYP/6-31G*//PM6, 86 M molecules).
**Result.** Every "PubChemQC KS gap" in -01…-19 is a B3LYP/6-31G* HOMO–LUMO gap at the PM6-optimised geometry.
**Interpretation.** Wording only changes, no number: the comparison geometry in -03 ("the PubChemQC geometry is not on
disk") is a PM6 geometry, not a B3LYP one. It also gives the most likely cause of the defective records found in -09:
PubChemQC documents PM6 optimisations that end on a structure different from the PubChem identifier (Nakata 2020,
10.1021/acs.jcim.0c00740; OGB PCQM4Mv2 ~10 % SMILES↔gap mismatch). The manuscript cites the 2023 B3LYP//PM6 dataset.
**Caveats.** The PM6 geometries themselves are not on disk, so the mechanism is inferred, not shown, for 9168.
**Claim link.** C5, C6 (mechanism), Methods (data source).

## L-2026-10-18-21 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (codex Methods check, 5 mismatches fixed). Verified against: see archived entry.

## L-2026-10-18-22 — ROTATED → `outputs/analysis/archive/ledger-2026-10.md` (copilot Results recompute, one count error (17 not 15)). Verified against: see archived entry.

## L-2026-10-18-23 — CORRECTION from the independent SI check (codex): chlorine classes, naphthalene midpoint, window-shift label, order wording

**Method.** codex (gpt-6.1-sol, effort high) recomputed every SI number from the data with its own code →
`outputs/cli-agent/L2-codex-si-check.md` (VERDICT MISMATCH; ST1–ST3 reproduce). Each DIFF re-derived here from the data:
RDKit Cl atoms of degree 0 (isolated fragment) vs degree > 0 (covalent), non-exclusive, on `ks_vs_stda_joined.csv` ⋈
`gate_cascade.csv` SMILES; HOPV LOO residuals (n 34); Malloci ranges in `malloci2007-La.txt`.
**Verified against.** dataset SMILES; `hopv15_*.csv`; `malloci2007-La.txt`; `std2-assessment.md`, `pentacene-fullTDA-validation.txt`.
**Result.** (-05) covalent Cl (non-exclusive): shard 14.7 %, KS window **6.1 % (2/33: 20579, 23804, both hydrochloride
salts)**, S1 window 13.4 %, S1 < 1 eV 12.0 %; the old "Cl⁻ counter-ion" column = any isolated Cl fragment (incl. neutral
HCl encodings): 2.9 / 36.4 / 5.9 / 53.9 %; strictly anionic [Cl-] 1.5 / 36.4 / 2.9 / 27.4 %. (-12) the ±0.097/0.115 eV window
shifts are the sTDA-xTB LOO MAE and the KS LOO RMSE; the KS LOO RMSE (0.115 eV) alone gives 34–58 around 46 (counts
unchanged). Naphthalene La midpoint of 4.44–4.45 = **4.445** (was 4.45): calibrated-KS acene MAE **0.40** (was 0.41),
naphthalene calibrated error −0.71 (was −0.72); T1, F3 regenerated. std2 vs full TDA on pentacene: 2.094 vs 2.203 eV.
F1 (Lean): order never changes the final survivor set, also for gates whose observable is invalid on some records.
**Interpretation.** No conclusion changes: covalent Cl is not enriched near either window, and the halogen excess still
follows the salts. Three statements were wrong as printed and are fixed in the manuscript: "covalent Cl absent from the KS
window", "order matters where an observable is invalid" (Results 2.5 now: order changes only per-gate credit; the validity
gate's presence is what matters), and "std2 within published accuracy" (Methods now gives the numbers). The geometry-corrected
sTDA-xTB class contrast (~0.5 eV, -15) is an inference, not a calculation, and is no longer stated as a number.
**Caveats.** Bootstrap upper bound 0.091 (seed 2026, 10⁴) vs 0.090 (codex, 2×10⁵): Monte-Carlo variation, kept.
**Claim link.** C7 (salts), C9 (F1 wording), C10, C12 (acene MAE 0.40).

## L-2026-10-18-24 — L3 verification of the death-angle queue (A4, A6) and of the checkers' findings (codex, copilot)

**Method.** Direct computation on `gate_cascade.csv` (+ `pipeline/scharber_recompute.py` AM1.5G J_SC, FF 0.65) for A6 →
`outputs/extract/window-redundancy.txt`; SMILES read for A4; `record-validity-result.csv` recounted (signed and absolute
deviation); Crossref API for every DOI a CLI returned; `applicability.py:60` read for the background-sample definition.
**Verified against.** the CSVs above and the script line cited (primary data), not the checkers' summaries.
**Result.** (A6) Without the window the cascade (G0, G1, G3, G4) keeps 10,153 molecules; the window keeps 37. Of the 37 the
PCE range is 1.86–10.14 %; **566 molecules with PCE ≥ 1.86 % lie outside the window** (31 with PCE ≥ 5 %, E_opt 2.23–2.58 eV,
max 8.34 %), and only 18 of the 37 highest-PCE molecules without a window sit inside it. (A4) The top PCE survivors 17574
(10.14 %), 22936, 22916 are sulfonated amino-anthraquinones (acid dyes); 12919 and 12913 are nitroso/azo-quinone
structures. (Record check) `record-validity-summary.txt` prints a two-sided count (23/29 flag, 8/29 control, |dev| > 0.5 eV);
the signed criterion used in the text (recomputed − stored > +0.5 eV) gives 21/29 and 0/29, p 1.28e-9, reproduced; the 8
controls are all *negative* deviations (UFF geometry lowers gaps). (Priors) Pyzer-Knapp 2016, Lopez 2017, Zhao 2021, Greenstein
2023, Omar 2022 resolve in Crossref with matching titles/authors. (Mislabel) `applicability.py:60` draws the 2,000 background records from single-fragment records, not from G0-valid ones.
**Interpretation.** The window is not redundant with the Scharber efficiency: it is the dominant selection step, it excludes
hundreds of molecules with equal or higher model PCE, and its edge at 2.2 eV is a convention, so "37 survivors" is a property of
that convention (already stated with the ±25 % calibration range in `-12`; the exclusion of higher-PCE molecules is new and goes to
SI S5). The best three survivors being dyes confirms that the cascade removes record errors, not chemistry that a chemist would
reject; the text now says so. The flag precision stands; only the summary file's two-sided wording differs and is not quoted.
**Caveats.** A6 uses the strict convention and the calibrated optical gap; PCE is a model upper bound. A4 availability of the dyes
was not looked up (off the central question). Prior-art search by codex is bounded, not proof of absence.
**Claim link.** C4, C11 (wording), novelty card (P4 Pyzer-Knapp 2016 added), SI S5.

## L-2026-10-18-25 — L3 corrections applied to the manuscript (checker findings verified before use)

**Method.** Each finding of `outputs/cli-agent/L3-*.md` re-checked against the manuscript line and the data; accepted or rejected below.
**Verified against.** `main.tex`, `sections/*.tex`, `cover_letter.tex` (lines quoted in `outputs/critical-reviews/L3-findings.md`); abstract
word count recomputed (193 words with hyphenated compounds as one word, 200 with en-dashes split).
**Result.** Applied: abstract names what is new and scopes the negative result ("for the classes and block tested"); 'every
observable about 0.1 eV' → '0.09–0.15 eV'; 'close to exact/nearly exact' → 'accurate' for the four acenes; 'before any chemistry
enters' removed; 1,506 records without spectra stated in the Discussion; figure caption 'within 0.1 eV' → LOO MAE and '2,000 valid
records' → 'single-fragment records'; Pyzer-Knapp 2016 (+ Lopez 2017, Omar 2022) cited in the Introduction; dye wording; cover
letter 'every high-throughput screen' and the categorical 'no' scoped. Rejected: copilot's '202-word abstract' (recount 193–200,
within the 200 limit); kiro's quote 'Three elements are new' located in the *abstract* (it is in the Introduction, line 44 —
the checker mis-attributed it; the abstract now carries the statement); kiro-cli overclaim pass void (monthly quota) and re-run on copilot.
**Interpretation.** The triple-read reviewer FAIL was right in substance: the abstract did not state what is new. Two checker
errors (word count, quote location) were caught only because each finding was re-opened against the file (rule 18).
**Caveats.** The attrition 39 → 34 HOPV molecules (3 DSSC, 2 UFF embed failures) was not stated; now in Results 2.2 (recount from hopv15_records.csv + hopv15_stda.csv).
**Claim link.** Abstract, Introduction, Discussion, Results 2.1/2.2/2.5, Fig. 5 caption, cover letter.
