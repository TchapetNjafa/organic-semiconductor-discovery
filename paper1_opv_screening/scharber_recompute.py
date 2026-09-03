#!/usr/bin/env python3
"""
scharber_recompute.py — Paper 1, single documented re-runnable Scharber stage.

WHY THIS SCRIPT EXISTS
----------------------
The pipeline used for the original submission (DD-ART-07-2026-000487) was recovered
from `ZENODO_UPLOAD/scharber_pce_calculation.ipynb` and reproduces the deposited
columns to 1e-14 (see outputs/extract/referee-verification-10.txt, [V44]). It does
NOT implement what the manuscript's Methods section describes:

  quantity   Methods claims                          notebook actually does
  ---------  --------------------------------------  ------------------------------------
  J_SC       AM1.5G photon flux integrated above     min(433.116 * exp(-Eg^2/2.335),
             the gap, EQE = 0.65                     415.225) -- a fitted Gaussian; no
                                                     spectrum, no integral, no EQE
  P_in       100 mW/cm^2 = 1000 W/m^2                900.139 W/m^2  (+11.09 % on every PCE)
  FF         "empirical V_OC-dependent expression,   Green's FF with ideality n = 2
             following the original Scharber          hard-coded (Scharber 2006 fixes
             parameterization"                        FF = 0.65)
  donor      HOMO > -6.1 AND LUMO > -3.7             LUMO >= -4.0 only; no HOMO test;
             (LUMO_PCBM = -3.7)                       LUMO_PCBM = -4.3
  acceptor   HOMO < -5.5 AND LUMO < -3.6             LUMO <= -3.9 only; no HOMO test

This script replaces that stage with an implementation that matches its own
documentation, line by line, and computes every number the manuscript quotes.

CONSTRUCTION (every choice stated; nothing implicit)
---------------------------------------------------
Reference materials, from the manuscript (unchanged):
    PCBM    HOMO = -6.10 eV   LUMO = -3.70 eV   E_g(opt) = 1.80 eV [electrochemical gap]
    PCDTBT  HOMO = -5.50 eV   LUMO = -3.60 eV   E_g(opt) = 1.90 eV

Two pairings are evaluated for every molecule:
    scenario D : molecule = DONOR    , PCBM   = acceptor
    scenario A : molecule = ACCEPTOR , PCDTBT = donor

Interfacial charge-transfer energy (Scharber 2006, eq. 1):
    E_CT = |HOMO_donor| - |LUMO_acceptor|                                   [eV]
Open-circuit voltage, with the empirical 0.3 V CT-state loss:
    V_OC = E_CT/e - 0.3                                                     [V]

Short-circuit current density, from the real spectrum:
    J_SC = q * EQE * INT_{0}^{lambda_g} Phi(lambda) d lambda                [A/m^2]
    lambda_g = h c / E_g(absorber) ;  EQE = 0.65 (constant, as the manuscript states)
    Phi = AM1.5G spectral photon flux, from ASTM G173-03 global tilt irradiance.
    Spectrum file: data/ASTMG173.csv (provenance recorded in that file's header note).

Fill factor -- reported BOTH ways, because the paper's own two statements disagree
and the sign of the conclusion depends on the choice:
    FF_scharber = 0.65                       (fixed, as in Scharber 2006)
    FF_green    = (v - ln(v + 0.72))/(v + 1) , v = V_OC/(n kT/q), n = 1, T = 300 K

    PCE = 100 * V_OC * J_SC * FF / P_in ,     P_in = 1000.0 W/m^2   [AM1.5G, exact]

PHYSICAL GATES (each is one line and each is reported separately)
    G2  exciton dissociation must not be uphill:    E_CT <= E_g(mol)
    G3  positive photovoltage:                      V_OC > 0
    G4  minimum driving force for charge separation, as required by Scharber 2006:
        molecule as donor    : LUMO_mol - LUMO_PCBM   >= 0.3 eV
        molecule as acceptor : LUMO_PCDTBT - LUMO_mol >= 0.3 eV
        This gate exists in the original notebook (`lumo_offset >= 0.3`) but with
        LUMO_PCBM = -4.3 eV instead of the -3.7 eV the Methods section states. Here it
        uses the stated value. Without G4 the construction returns PCE up to 14.7 %,
        above the ~11 % Scharber single-junction optimum, because V_OC is then allowed
        to reach E_g - 0.3 rather than E_g - 0.6.
    plus one of two ABSORBER CONVENTIONS, both computed, neither privileged:

    STRICT        G1s: E_g(mol) <= E_g(partner)      molecule IS the narrower-gap
                                                     component and sets the onset
                  J_SC integrated above E_g(mol)
    COMPLEMENTARY G1c: EG_LO <= E_g(mol) <= EG_HI    molecule absorbs visible light;
                                                     the blend onset is the narrower of
                                                     the two components
                  J_SC integrated above min(E_g(mol), E_g(partner))

    Why both. STRICT is the conservative reading: only credit a molecule for photocurrent
    it can itself generate. COMPLEMENTARY is how real bulk heterojunctions work (a
    wide-gap donor with a narrow-gap acceptor, or the reverse) -- but it must be paired
    with a visible-window requirement on the molecule, because integrating at min(E_g)
    with no such requirement lets the fullerene supply the current and hands high scores
    to molecules transparent across the solar spectrum (verified: 4,016 such artifacts,
    outputs/extract/referee-verification-9.txt [V39]).
    A molecule that ranks first under both conventions and both FF choices is robust to
    this modelling freedom; one that does not, is not.

CHEMICAL GATES (unchanged from the submitted manuscript)
    validity R1-R6 : single covalent species, no metal/metalloid, >=1 aromatic ring,
                     >=6 heavy atoms, >=6 conjugated atoms
    stability      : no azide / nitroso / diazo (SMARTS, Table S3)

OUTPUT
    outputs/recompute/screen_full.csv     one row per molecule, every intermediate
    outputs/recompute/summary.json        every headline number the manuscript quotes
    stdout                                a human-readable audit trail

Run:  python scharber_recompute.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# ----------------------------------------------------------------- paths
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                      # ARTICLE_MVOTO/


def _first_existing(*candidates):
    """Return the first path that exists, else the first candidate (for the error)."""
    for c in candidates:
        if os.path.exists(c):
            return c
    return candidates[0]


# The property table lives in the working tree at ../ZENODO_UPLOAD/; in the
# published deposit the same file sits beside this script under data/. Both are
# accepted so the deposit runs standalone.
DATASET = _first_existing(
    os.path.join(HERE, "data", "dataset_pubchemqc_opv_17458.csv"),
    os.path.join(ROOT, "ZENODO_UPLOAD", "dataset_pubchemqc_opv_17458.csv"),
)
SPECTRUM = os.path.join(HERE, "data", "ASTMG173.csv")
OUTDIR = os.path.join(HERE, "outputs", "recompute")
os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------------------------------------------- constants
Q = 1.602176634e-19          # C            (SI 2019 exact)
H = 6.62607015e-34           # J s          (SI 2019 exact)
C = 2.99792458e8             # m/s          (SI 2019 exact)
KB = 1.380649e-23            # J/K          (SI 2019 exact)
HC_EV_NM = H * C / Q * 1e9   # 1239.842 eV nm

EQE = 0.65                   # constant external quantum efficiency, as stated
P_IN = 1000.0                # W/m^2, AM1.5G
V_LOSS = 0.30                # V, empirical CT-state loss (Scharber 2006)
DRIVING_MIN = 0.30           # eV, minimum LUMO-LUMO offset for charge separation
FF_FIXED = 0.65              # Scharber 2006
IDEALITY = 1.0               # Green's FF ideality factor
TEMP = 300.0                 # K
EG_LO, EG_HI = 1.1, 2.2      # eV, visible-harvesting window (gate G4)

REFERENCES = {
    # name:        (HOMO, LUMO, optical gap used as the partner's absorption onset)
    "PCBM":   (-6.10, -3.70, 1.80),
    "PCDTBT": (-5.50, -3.60, 1.90),
}

METALS = set(
    "Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge Rb Sr Y Zr Nb Mo Tc Ru "
    "Rh Pd Ag Cd In Sn Sb Cs Ba La Ce Pt Au Hg Tl Pb Bi B Si As Te".split()
)
REACTIVE_SMARTS = {
    "azide":   ["[#7]-[#7+]#[#7]", "[#7]=[#7+]=[#7-]", "[#7-]-[#7+]#[#7]"],
    "nitroso": ["[#6][NX2]=[OX1]"],
    "diazo":   ["[#6]=[#7+]=[#7-]"],
}


# ----------------------------------------------------------------- spectrum
def load_am15g(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (wavelength_nm, spectral photon flux in photons m^-2 s^-1 nm^-1)."""
    df = pd.read_csv(path, skiprows=1)
    lam = df["wavelength"].to_numpy(float)               # nm
    irr = df["global"].to_numpy(float)                   # W m^-2 nm^-1
    photon_energy = HC_EV_NM / lam * Q                   # J per photon
    return lam, irr / photon_energy


def jsc_above_gap(lam: np.ndarray, flux: np.ndarray, eg_ev: np.ndarray) -> np.ndarray:
    """q * EQE * cumulative photon flux for lambda <= hc/Eg, vectorised over eg_ev.

    Uses the exact cumulative trapezoid of the tabulated spectrum, then linearly
    interpolates at the cut-off wavelength. Monotone increasing in lambda_g, so
    monotone DECREASING in E_g -- the physically expected behaviour.
    """
    cum = np.concatenate([[0.0], np.cumsum(np.diff(lam) * 0.5 * (flux[1:] + flux[:-1]))])
    lam_g = np.where(np.isfinite(eg_ev) & (eg_ev > 0), HC_EV_NM / eg_ev, np.nan)
    integrated = np.interp(lam_g, lam, cum, left=0.0, right=cum[-1])
    return Q * EQE * integrated


def ff_green(voc: np.ndarray, n: float = IDEALITY, T: float = TEMP) -> np.ndarray:
    vt = KB * T / Q
    with np.errstate(divide="ignore", invalid="ignore"):
        v = voc / (n * vt)
        ff = (v - np.log(v + 0.72)) / (v + 1.0)
    return np.clip(np.nan_to_num(ff, nan=0.0), 0.0, 1.0)


# ----------------------------------------------------------------- chemistry
def validity_reason(smiles: str) -> str:
    """'' if the molecule is a plausible organic semiconductor, else the failing rule."""
    if not isinstance(smiles, str) or not smiles:
        return "unparsable"
    if "." in smiles:
        return "R1_multifragment"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unparsable"
    if {a.GetSymbol() for a in mol.GetAtoms()} & METALS:
        return "R2_metal"
    if mol.GetRingInfo().NumRings() == 0:
        return "R3_no_ring"
    if not any(a.GetIsAromatic() for a in mol.GetAtoms()):
        return "R4_no_aromatic_ring"
    if mol.GetNumHeavyAtoms() < 6:
        return "R5_too_small"
    conjugated = sum(
        1 for a in mol.GetAtoms()
        if a.GetIsAromatic() or a.GetHybridization() == Chem.HybridizationType.SP2
    )
    if conjugated < 6:
        return "R6_low_conjugation"
    return ""


_PATTERNS = {k: [Chem.MolFromSmarts(s) for s in v] for k, v in REACTIVE_SMARTS.items()}


def reactive_groups(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return ""
    return ",".join(
        name for name, patts in _PATTERNS.items()
        if any(mol.HasSubstructMatch(p) for p in patts)
    )


# ----------------------------------------------------------------- pipeline
def main() -> int:
    log: list[str] = []

    def say(msg: str = "") -> None:
        log.append(msg)
        print(msg)

    say("=" * 78)
    say("Paper 1 -- Scharber stage, recomputed from a documented implementation")
    say("=" * 78)

    # --- spectrum -----------------------------------------------------
    lam, flux = load_am15g(SPECTRUM)
    total_power = np.trapezoid(flux * (HC_EV_NM / lam * Q), lam)
    say(f"\n[1] AM1.5G spectrum: {len(lam)} points, {lam.min():.0f}-{lam.max():.0f} nm")
    say(f"    integrated irradiance = {total_power:.2f} W/m^2   "
        f"(ASTM G173-03 nominal 1000.37; deviation {abs(total_power-1000.37):.2f})")
    for eg in (1.1, 1.4, 1.9, 2.2, 3.0):
        say(f"    J_SC(E_g = {eg:.1f} eV, EQE = {EQE}) = "
            f"{float(jsc_above_gap(lam, flux, np.array([eg]))[0]):7.2f} A/m^2")
    say("    (Scharber 2006 reports ~200 A/m^2 at E_g = 1.4 eV for EQE = 0.65)")

    # --- dataset ------------------------------------------------------
    df = pd.read_csv(DATASET, low_memory=False)
    n_all = len(df)
    say(f"\n[2] dataset: {DATASET.split(os.sep)[-1]}  ->  {n_all} molecules")
    homo = pd.to_numeric(df["HOMO(eV)"], errors="coerce").to_numpy(float)
    lumo = pd.to_numeric(df["LUMO(eV)"], errors="coerce").to_numpy(float)
    gap = pd.to_numeric(df["GAP(eV)"], errors="coerce").to_numpy(float)
    sas = pd.to_numeric(df["sas1(%)"], errors="coerce").to_numpy(float)

    # --- two scenarios x two absorber conventions ----------------------
    say("\n[3] scenarios, gates and photovoltaic quantities")
    out: dict[str, np.ndarray] = {}
    for conv, conv_label in (("strict", "STRICT: molecule is the narrower-gap absorber"),
                             ("compl", "COMPLEMENTARY: blend onset = min(E_g), molecule in window")):
        say(f"\n    === {conv_label} ===")
        for tag, role in (("D", "donor"), ("A", "acceptor")):
            if role == "donor":                                # molecule donates
                partner, (p_homo, p_lumo, p_eg) = "PCBM", REFERENCES["PCBM"]
                e_ct = np.abs(homo) - abs(p_lumo)
                driving = lumo - p_lumo                        # LUMO_mol - LUMO_PCBM
            else:                                              # molecule accepts
                partner, (p_homo, p_lumo, p_eg) = "PCDTBT", REFERENCES["PCDTBT"]
                e_ct = abs(p_homo) - np.abs(lumo)
                driving = p_lumo - lumo                        # LUMO_PCDTBT - LUMO_mol

            if conv == "strict":
                g1 = gap <= p_eg                               # molecule sets the onset
                eg_absorb = gap
                g1_text = f"G1s E_g(mol) <= {p_eg} eV"
            else:
                g1 = (gap >= EG_LO) & (gap <= EG_HI)           # molecule harvests visible
                eg_absorb = np.minimum(gap, p_eg)              # blend onset
                g1_text = f"G1c {EG_LO} <= E_g(mol) <= {EG_HI} eV"

            g2 = e_ct <= gap                                   # dissociation not uphill
            voc = e_ct - V_LOSS
            g3 = voc > 0
            g4 = driving >= DRIVING_MIN                        # LUMO-LUMO offset
            keep = g1 & g2 & g3 & g4

            jsc = jsc_above_gap(lam, flux, eg_absorb)
            ffg = ff_green(voc)
            pce_fixed = np.where(keep, 100.0 * voc * jsc * FF_FIXED / P_IN, np.nan)
            pce_green = np.where(keep, 100.0 * voc * jsc * ffg / P_IN, np.nan)

            say(f"    scenario {tag} (molecule as {role}, partner {partner}, "
                f"E_g partner {p_eg} eV)")
            say(f"      {g1_text:<38} : {int(np.sum(g1)):>6}")
            say(f"      G2 E_CT <= E_g(mol)                    : {int(np.sum(g2)):>6}")
            say(f"      G3 V_OC > 0                            : {int(np.sum(g3)):>6}")
            say(f"      G4 LUMO offset >= {DRIVING_MIN} eV            : {int(np.sum(g4)):>6}")
            say(f"      all gates                              : {int(np.sum(keep)):>6}")
            out.update({
                f"{conv}_e_ct_{tag}": e_ct,
                f"{conv}_voc_{tag}": np.where(keep, voc, np.nan),
                f"{conv}_jsc_{tag}": np.where(keep, jsc, np.nan),
                f"{conv}_pce_fixed_{tag}": pce_fixed,
                f"{conv}_pce_green_{tag}": pce_green,
                f"{conv}_pass_{tag}": keep,
            })

        # best scenario per molecule, for each FF convention
        for ff_key in ("fixed", "green"):
            d = out[f"{conv}_pce_{ff_key}_D"]
            a = out[f"{conv}_pce_{ff_key}_A"]
            best = np.fmax(np.nan_to_num(d, nan=-np.inf), np.nan_to_num(a, nan=-np.inf))
            best = np.where(np.isfinite(best), best, np.nan)
            out[f"{conv}_pce_{ff_key}"] = best
            out[f"{conv}_ps_{ff_key}"] = best - sas
        passes = out[f"{conv}_pass_D"] | out[f"{conv}_pass_A"]
        out[f"{conv}_pass"] = passes
        out[f"{conv}_role"] = np.where(
            passes,
            np.where(np.nan_to_num(out[f"{conv}_pce_fixed_A"], nan=-np.inf)
                     > np.nan_to_num(out[f"{conv}_pce_fixed_D"], nan=-np.inf),
                     "acceptor/PCDTBT", "donor/PCBM"),
            "gated out")
        say(f"      passing in at least one scenario        : {int(passes.sum()):>6}"
            f"  ({100.0*passes.sum()/n_all:.2f} %)")
        for ff_key, ff_label in (("fixed", "FF = 0.65"), ("green", "FF = Green n=1")):
            v = out[f"{conv}_pce_{ff_key}"]
            say(f"      {ff_label:<16} max PCE {np.nanmax(v):6.2f} %   "
                f"n(PCE > 10 %) = {int(np.nansum(v > 10)):>4}")

    # --- chemistry ----------------------------------------------------
    say("\n[4] chemical gates")
    smiles = df["SMILES"].astype(str).tolist()
    reject = np.array([validity_reason(s) for s in smiles])
    react = np.array([reactive_groups(s) for s in smiles])
    valid, stable = reject == "", react == ""
    say(f"    chemically valid                          : {int(valid.sum()):>6}"
        f"   rejected {int((~valid).sum())}")
    say(f"    valid and free of reactive groups         : {int((valid & stable).sum()):>6}")
    vc = pd.Series(reject[~valid]).value_counts()
    for rule, count in vc.items():
        say(f"      {rule:<22} {count:>6}")

    # --- assemble -----------------------------------------------------
    cols: dict[str, object] = {
        "mol_id": df["mol_id"], "formula": df["formula"], "SMILES": df["SMILES"],
        "mass": df["mass"], "HOMO_eV": homo, "LUMO_eV": lumo, "gap_eV": gap,
        "dipole_D": pd.to_numeric(df["dipole_moment"], errors="coerce"),
        "SAScore": sas,
        "E_CT_donor_eV": out["strict_e_ct_D"], "E_CT_acceptor_eV": out["strict_e_ct_A"],
        "validity_reject": reject, "reactive_groups": react,
        "chemically_valid": valid, "stable": stable,
        "in_visible_window": (gap >= EG_LO) & (gap <= EG_HI),
    }
    # four (absorber convention x FF convention) variants, each fully traceable
    VARIANTS = [("strict", "fixed", "strict_FF065"), ("strict", "green", "strict_FFgreen"),
                ("compl", "fixed", "compl_FF065"), ("compl", "green", "compl_FFgreen")]
    for conv, ff_key, label in VARIANTS:
        cols[f"PCE_{label}"] = out[f"{conv}_pce_{ff_key}"]
        cols[f"PCE_SAScore_{label}"] = out[f"{conv}_ps_{ff_key}"]
    for conv in ("strict", "compl"):
        cols[f"VOC_donor_{conv}_V"] = out[f"{conv}_voc_D"]
        cols[f"VOC_acceptor_{conv}_V"] = out[f"{conv}_voc_A"]
        cols[f"JSC_donor_{conv}_Am2"] = out[f"{conv}_jsc_D"]
        cols[f"JSC_acceptor_{conv}_Am2"] = out[f"{conv}_jsc_A"]
        cols[f"passes_{conv}"] = out[f"{conv}_pass"]
        cols[f"role_{conv}"] = out[f"{conv}_role"]
    res = pd.DataFrame(cols)
    res.to_csv(os.path.join(OUTDIR, "screen_full.csv"), index=False)
    say(f"\n[5] wrote {os.path.join('outputs', 'recompute', 'screen_full.csv')} "
        f"({len(res)} rows, {len(res.columns)} columns)")

    # --- ranking (the paper's claim is a RANKING, not a threshold) ----
    say("\n[6] RANKING of chemically valid, stable, physically admissible molecules")
    say("    Four variants: 2 absorber conventions x 2 fill-factor conventions.")
    summary: dict[str, object] = {
        "n_molecules": n_all,
        "n_valid": int(valid.sum()), "n_valid_stable": int((valid & stable).sum()),
        "spectrum_integrated_W_m2": round(float(total_power), 3),
        "jsc_at_1p4eV_Am2": round(float(jsc_above_gap(lam, flux, np.array([1.4]))[0]), 2),
        "constants": {"EQE": EQE, "P_in_W_m2": P_IN, "V_loss_V": V_LOSS,
                      "driving_force_min_eV": DRIVING_MIN,
                      "FF_fixed": FF_FIXED, "FF_green_ideality": IDEALITY,
                      "Eg_window_eV": [EG_LO, EG_HI]},
        "references": REFERENCES,
        "variants": {},
    }
    first_places: dict[str, int] = {}
    for conv, ff_key, label in VARIANTS:
        pool = res[res[f"passes_{conv}"] & res.chemically_valid & res.stable].copy()
        col = f"PCE_SAScore_{label}"
        top = pool.nlargest(10, col)
        say(f"\n    --- {label}   (absorber: {conv}, FF: "
            f"{'0.65 fixed' if ff_key == 'fixed' else 'Green n=1'}) ---")
        say(f"    admissible pool = {len(pool)}   "
            f"max PCE over all molecules = {np.nanmax(res[f'PCE_{label}']):.2f} %")
        say("    rank  id      formula          PCE_SAScore   PCE     SAScore  E_g")
        for rank, (_, r) in enumerate(top.iterrows(), 1):
            say(f"    {rank:>4}  {int(r.mol_id):<7} {str(r.formula):<16} "
                f"{r[col]:+11.2f} {r[f'PCE_{label}']:7.2f}  {r.SAScore:7.2f}  {r.gap_eV:5.3f}")
        entry: dict[str, object] = {
            "pool": len(pool),
            "max_pce_all": round(float(np.nanmax(res[f"PCE_{label}"])), 3),
            "top10": [int(x) for x in top.mol_id],
        }
        if len(top):
            first_places[label] = int(top.iloc[0].mol_id)
            entry["first"] = int(top.iloc[0].mol_id)
        r1712 = pool[pool.mol_id == 1712]
        if len(r1712):
            rr = r1712.iloc[0]
            rank_1712 = int((pool[col] > rr[col]).sum()) + 1
            say(f"    molecule 1712 (submitted lead): rank {rank_1712} / {len(pool)}   "
                f"PCE {rr[f'PCE_{label}']:.2f} %   PCE_SAScore {rr[col]:+.2f}")
            entry.update({"rank_1712": rank_1712,
                          "pce_1712": round(float(rr[f"PCE_{label}"]), 3),
                          "ps_1712": round(float(rr[col]), 3)})
        else:
            say("    molecule 1712 (submitted lead): NOT in the admissible pool")
            entry["rank_1712"] = None
        summary["variants"][label] = entry

    say("\n    first place by variant:")
    for label, mid_ in first_places.items():
        say(f"      {label:<16} -> {mid_}")
    robust = len(set(first_places.values())) == 1
    say(f"    -> ranking is {'ROBUST' if robust else 'NOT robust'} to the two modelling "
        f"freedoms ({len(set(first_places.values()))} distinct first places)")
    summary["ranking_robust"] = robust
    summary["first_places"] = first_places

    # --- the cascade, for reference and for the orthogonality argument -
    say("\n[7] threshold cascade (for continuity with the submitted version;")
    say("    the paper's claim is now the ranking above, not this threshold)")
    for conv, ff_key, label in VARIANTS:
        col = f"PCE_SAScore_{label}"
        raw = res[res[col] > 0]
        g_valid = raw[raw.chemically_valid]
        g_stable = g_valid[g_valid.stable]
        say(f"    {label:<16} PCE_SAScore > 0: {len(raw):>4} -> valid {len(g_valid):>4}"
            f" -> stable {len(g_stable):>4}   "
            f"ids={sorted(int(x) for x in g_stable.mol_id)[:12]}")
        summary["variants"][label]["cascade"] = [len(raw), len(g_valid), len(g_stable)]

    say("\n[8] orthogonality control: does a gap pre-filter remove the artifacts?")
    win = res[res.in_visible_window]
    say(f"    molecules with E_g in [{EG_LO}, {EG_HI}] eV: {len(win)}")
    arts = res[res.mol_id.isin([977, 11029, 7801])][
        ["mol_id", "formula", "gap_eV", "in_visible_window", "validity_reject"]]
    for _, r in arts.iterrows():
        say(f"      {int(r.mol_id):<6} {str(r.formula):<10} E_g {r.gap_eV:5.3f} eV  "
            f"inside window: {bool(r.in_visible_window)}   ({r.validity_reject})")
    n_art_in = int(arts.in_visible_window.sum())
    say(f"    -> {n_art_in} of 3 non-semiconductor artifacts lie INSIDE the window;")
    say("       an energetic pre-filter cannot substitute for a structural one.")
    summary["artifacts_inside_gap_window"] = n_art_in

    # --- robustness across the four conventions -------------------------
    # regen_tables.py and regen_figures.py consume this file; it must be written
    # here so the deposit regenerates every artifact from the inputs alone.
    say("\n[9] robustness of the surviving set across the four conventions")
    sets = {}
    for conv, ff_key, label in VARIANTS:
        col = f"PCE_SAScore_{label}"
        g = res[(res[col] > 0) & res.chemically_valid & res.stable]
        sets[label] = {int(x) for x in g.mol_id}
        say(f"    {label:<16} above threshold, valid and stable: "
            f"{sorted(sets[label])}")
    inter = sorted(set.intersection(*sets.values())) if sets else []
    union = sorted(set.union(*sets.values())) if sets else []
    say(f"    intersection of all four: {inter}")
    say(f"    union of all four       : {union}")

    # Why the previously reported lead is not admitted, stated from the numbers.
    r1712 = res[res.mol_id == 1712]
    if len(r1712):
        L = float(r1712.iloc[0]["LUMO_eV"])
        d_pcbm = L - REFERENCES["PCBM"][1]
        d_pcdtbt = REFERENCES["PCDTBT"][1] - L
        reason_1712 = (
            f"LUMO {L:.3f} eV gives {d_pcbm:+.3f} eV driving force vs PCBM "
            f"({REFERENCES['PCBM'][1]:.2f}) and {d_pcdtbt:+.3f} eV vs PCDTBT "
            f"({REFERENCES['PCDTBT'][1]:.2f}); both below the "
            f"{DRIVING_MIN} eV requirement")
    else:
        reason_1712 = "molecule 1712 is not present in this dataset"
    say(f"    molecule 1712: {reason_1712}")

    robustness = {
        "intersection_all_variants": inter,
        "union_all_variants": union,
        "mol_1712_excluded_reason": reason_1712,
        "robust_lead": inter,
        "per_variant_stable_viable": {k: sorted(v) for k, v in sets.items()},
    }
    with open(os.path.join(OUTDIR, "robustness.json"), "w") as fh:
        json.dump(robustness, fh, indent=2)

    with open(os.path.join(OUTDIR, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    with open(os.path.join(OUTDIR, "recompute_log.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    say("\n[10] wrote outputs/recompute/{summary.json, robustness.json, recompute_log.txt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
