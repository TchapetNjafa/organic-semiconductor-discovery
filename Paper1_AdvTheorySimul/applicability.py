"""Applicability domain of the HOPV15 gap calibration (death-angle angle 2).

Compares the HOPV15 calibration molecules with the cascade survivors and the valid shard in a small, chemically
readable descriptor space, and checks whether the calibration is used inside the KS-gap range it was fitted on.
Domain test: Mahalanobis distance to the HOPV set (standardised descriptors); "inside" = below the 95th
percentile of HOPV's own leave-one-out distances.

usage: python -I applicability.py <project_root>
writes outputs/extract/applicability.{txt,csv}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
COLS = ["heavy_atoms", "frac_conj", "arom_rings", "frac_hetero", "n_S", "rot_bonds"]


def descriptors(smiles):
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return [np.nan] * len(COLS)
    n = m.GetNumHeavyAtoms()
    conj = sum(a.GetIsAromatic() or a.GetHybridization() == Chem.HybridizationType.SP2 for a in m.GetAtoms())
    het = sum(a.GetSymbol() not in ("C", "H") for a in m.GetAtoms())
    return [n, conj / n, rdMolDescriptors.CalcNumAromaticRings(m), het / n,
            sum(a.GetSymbol() == "S" for a in m.GetAtoms()), Descriptors.NumRotatableBonds(m)]


def frame(smiles):
    return pd.DataFrame([descriptors(s) for s in smiles], columns=COLS)


def mahalanobis(ref, x):
    mu, sd = ref.mean(), ref.std(ddof=1)
    z_ref, z_x = (ref - mu) / sd, (x - mu) / sd
    inv = np.linalg.pinv(np.cov(z_ref.T))
    return np.sqrt(np.einsum("ij,jk,ik->i", z_x, inv, z_x))


def loo_threshold(ref, q=95):
    d = [mahalanobis(ref.drop(index=i), ref.loc[[i]])[0] for i in ref.index]
    return np.percentile(d, q)


if __name__ == "__main__":
    root = Path(sys.argv[1])
    h = root / "outputs/extract/hopv15"
    hp = pd.read_csv(h / "hopv15_records.csv").merge(pd.read_csv(h / "hopv15_stda.csv"), on="idx")
    hp["gap_opt"] = pd.to_numeric(hp.gap_opt, errors="coerce")
    hp = hp[(hp.construction == "molecule") & ~hp.complement.str.contains("TiO2", na=False)
            & (hp.status == "success")].dropna(subset=["gap_opt"]).reset_index(drop=True)
    cas = pd.read_csv(root / "outputs/extract/gate_cascade.csv")
    from_flt = cas.SMILES.map(lambda s: "." not in s if isinstance(s, str) else False)
    sets = {"HOPV calibration": hp.smiles, "cascade survivors": cas.loc[cas.keep, "SMILES"],
            "valid shard (sample 2000)": cas.loc[from_flt].sample(2000, random_state=1).SMILES}
    D = {k: frame(v.tolist()).dropna().reset_index(drop=True) for k, v in sets.items()}
    ref = D["HOPV calibration"]
    thr = loo_threshold(ref)
    lines = [f"Mahalanobis 95th-pct LOO threshold (HOPV, n={len(ref)}): {thr:.2f}"]
    lines.append("set | n | " + " | ".join(f"{c} median [p5,p95]" for c in COLS) + " | inside domain")
    out = []
    for k, d in D.items():
        dist = mahalanobis(ref, d)
        inside = (dist <= thr).mean()
        stats = " | ".join(f"{d[c].median():.2f} [{d[c].quantile(.05):.2f},{d[c].quantile(.95):.2f}]" for c in COLS)
        lines.append(f"{k} | {len(d)} | {stats} | {100 * inside:.0f} %")
        out.append(d.assign(set=k, mahal=dist))
    ks_lo, ks_hi = hp.gap_b3lyp_eV.min(), hp.gap_b3lyp_eV.max()
    sv = cas.loc[cas.keep, "GAP(eV)"]
    lines.append(f"\nKS range of the calibration fit: [{ks_lo:.3f}, {ks_hi:.3f}] eV; survivors' KS "
                 f"[{sv.min():.3f}, {sv.max():.3f}] eV; survivors outside fit range: {((sv < ks_lo) | (sv > ks_hi)).sum()}/{len(sv)}")
    (root / "outputs/extract/applicability.txt").write_text("\n".join(lines) + "\n")
    pd.concat(out).to_csv(root / "outputs/extract/applicability.csv", index=False)
    print("\n".join(lines))
