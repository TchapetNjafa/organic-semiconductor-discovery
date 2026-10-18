"""Full gate cascade, donor vs PCBM, strict convention (molecule sets the absorption onset).

Order: validity (R1-R6) -> record-validity flag (S1 - KS > 1 eV) -> optical-gap window on the
HOPV15-calibrated KS gap -> Scharber gates (V_OC > 0, LUMO offset >= 0.3 eV) -> PCE.
V_OC and the LUMO offset use the stored KS frontier orbitals; J_SC uses the calibrated optical gap.
Reuses pipeline/scharber_recompute.py (constants, AM1.5G J_SC, Green FF) — no physics re-implemented.

usage: python -I gate_cascade.py <project_root>
writes outputs/extract/gate_cascade.{txt,csv}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "pipeline"))
import scharber_recompute as sc  # noqa: E402
from filter_genuine_osc import reject_reason  # noqa: E402

FLAG_EV = 1.0  # S1 - KS above this marks a suspect record (L-2026-10-17-09)


def calibration(root):
    h = root / "outputs/extract/hopv15"
    r = pd.read_csv(h / "hopv15_records.csv").merge(pd.read_csv(h / "hopv15_stda.csv"), on="idx")
    r["gap_opt"] = pd.to_numeric(r.gap_opt, errors="coerce")
    r = r[(r.construction == "molecule") & ~r.complement.str.contains("TiO2", na=False)
          & (r.status == "success")].dropna(subset=["gap_opt"])
    return np.polyfit(r.gap_b3lyp_eV, r.gap_opt, 1)


def cascade(root):
    d = pd.read_csv(sc.DATASET, usecols=["mol_id", "SMILES", "formula", "HOMO(eV)", "LUMO(eV)", "GAP(eV)"],
                    dtype={"mol_id": str})
    d["mol_id"] = pd.to_numeric(d.mol_id)
    s = pd.read_csv(root / "outputs/extract/stda_all.csv", dtype=str).apply(pd.to_numeric, errors="coerce")
    d = d.merge(s[["mol_id", "E_S1_eV"]], on="mol_id", how="left")
    k, b = calibration(root)
    d["eg_opt"] = k * d["GAP(eV)"] + b
    p_homo, p_lumo, _ = sc.REFERENCES["PCBM"]
    d["voc"] = (-d["HOMO(eV)"]) - (-p_lumo) - sc.V_LOSS        # |HOMO_D| - |LUMO_A| - 0.3
    d["drive"] = d["LUMO(eV)"] - p_lumo
    gates = {
        "G0 validity (R1-R6)": d.SMILES.map(reject_reason) == "",
        "G1 record flag (S1-KS <= 1 eV)": ~((d.E_S1_eV - d["GAP(eV)"]) > FLAG_EV),
        f"G2 optical gap in [{sc.EG_LO}, {sc.EG_HI}] eV (calibrated)": d.eg_opt.between(sc.EG_LO, sc.EG_HI),
        "G3 V_OC > 0": d.voc > 0,
        f"G4 LUMO offset >= {sc.DRIVING_MIN} eV": d.drive >= sc.DRIVING_MIN,
    }
    lam, flux = sc.load_am15g(sc.SPECTRUM)
    keep = np.ones(len(d), bool)
    lines = [f"calibration E_opt = {k:.3f}*KS + {b:.3f} (HOPV15); n records {len(d)}"]
    for name, g in gates.items():
        keep &= g.to_numpy()
        lines.append(f"{name:42s} alone {int(g.sum()):6d} | cumulative {int(keep.sum()):6d}")
    d["keep"] = keep
    jsc = sc.jsc_above_gap(lam, flux, d.eg_opt.to_numpy())
    d["pce_ff065"] = np.where(keep, 100 * d.voc * jsc * sc.FF_FIXED / sc.P_IN, np.nan)
    d["pce_green"] = np.where(keep, 100 * d.voc * jsc * sc.ff_green(d.voc.to_numpy()) / sc.P_IN, np.nan)
    return d, lines


if __name__ == "__main__":
    root = Path(sys.argv[1])
    d, lines = cascade(root)
    top = d[d.keep].sort_values("pce_ff065", ascending=False)
    lines.append(f"\nsurvivors {len(top)}; PCE(FF 0.65) max {top.pce_ff065.max():.2f} %, "
                 f"median {top.pce_ff065.median():.2f} %")
    lines.append(top[["mol_id", "formula", "GAP(eV)", "eg_opt", "E_S1_eV", "voc", "drive",
                      "pce_ff065", "pce_green"]].head(25).round(3).to_string(index=False))
    out = root / "outputs/extract"
    (out / "gate_cascade.txt").write_text("\n".join(lines) + "\n")
    d.to_csv(out / "gate_cascade.csv", index=False)
    print("\n".join(lines))
