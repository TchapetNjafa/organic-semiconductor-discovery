"""Tables T1 and T2 of the ATS manuscript, computed from the primary data only.

T1: error of each gap observable vs experiment on the 34 HOPV15 molecules where every observable exists
    (complete cases), plus the external acene class. T2: gate counts of the calibrated cascade.
usage: python -I make_tables.py <project_root>   -> manuscript/tables/T{1,2}.tex + tables_values.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "pipeline"))
from filter_genuine_osc import reject_reason  # noqa: E402

OBS = [("gap_b3lyp_eV", r"$\Eks$ (B3LYP)"), ("E_S1_eV", r"sTDA-xTB, UFF"),
       ("E_S1_dftgeom_eV", r"sTDA-xTB, DFT geom."), ("S1_stda", r"sTDA@B3LYP (std2)")]
ACENES = {931: 4.445, 8418: 3.36, 7080: 2.66, 8671: 2.215}  # Malloci et al. 2007, bright La midpoints


def loo(x, y, groups=None):
    g = np.arange(len(x)) if groups is None else groups
    e = np.empty(len(x))
    for G in np.unique(g):
        k = g != G
        e[g == G] = y[g == G] - np.polyval(np.polyfit(x[k], y[k], 1), x[g == G])
    return np.abs(e)


def t1(root):
    h = root / "outputs/extract/hopv15"
    m = (pd.read_csv(h / "hopv15_records.csv")
         .merge(pd.read_csv(h / "hopv15_stda.csv")[["idx", "E_S1_eV"]], on="idx")
         .merge(pd.read_csv(h / "hopv15_stda_dftgeom.csv")[["idx", "E_S1_dftgeom_eV"]], on="idx")
         .merge(pd.read_csv(h / "std2_reparsed.csv")[["idx", "S1_stda"]], on="idx"))
    m["gap_opt"] = pd.to_numeric(m.gap_opt, errors="coerce")
    m = m[(m.construction == "molecule") & ~m.complement.str.contains("TiO2", na=False)]
    m = m.dropna(subset=["gap_opt"] + [c for c, _ in OBS]).reset_index(drop=True)
    y = m.gap_opt.values
    rng = np.random.default_rng(2026)
    boot = rng.integers(0, len(y), (10000, len(y)))
    j = pd.read_csv(root / "outputs/extract/ks_vs_stda_joined.csv").set_index("mol_id")
    k, b = np.polyfit(m.gap_b3lyp_eV, y, 1)
    ac_exp = np.array(list(ACENES.values()))
    rows, vals = [], {"n": len(m), "calibration": [k, b]}
    for c, lab in OBS:
        x = m[c].values
        err = x - y
        lo = loo(x, y)
        ci = np.percentile(lo[boot].mean(1), [2.5, 97.5])
        lopo = loo(x, y, m.doi.values).mean()
        rec = int(((x >= 1.1) & (x <= 2.2) & (y >= 1.1) & (y <= 2.2)).sum())
        n_in = int(((y >= 1.1) & (y <= 2.2)).sum())
        ac = "--"
        if c == "gap_b3lyp_eV":
            raw = np.array([j.loc[i, "GAP(eV)"] for i in ACENES])
            ac = f"{np.abs(raw - ac_exp).mean():.2f} / {np.abs(k * raw + b - ac_exp).mean():.2f}"
        elif c == "E_S1_eV":
            ac = f"{np.abs(np.array([j.loc[i, 'E_S1_eV'] for i in ACENES]) - ac_exp).mean():.2f}"
        v = dict(bias=err.mean(), mae=np.abs(err).mean(), r=np.corrcoef(x, y)[0, 1], recall=f"{rec}/{n_in}",
                 loo=lo.mean(), loo_ci=list(ci), lopo=lopo, acene=ac)
        vals[c] = v
        rows.append(f"{lab} & {v['bias']:+.2f} & {v['mae']:.2f} & {v['r']:.2f} & {v['recall']} & "
                    f"{v['loo']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] & {lopo:.3f} & {ac} \\\\")
    tex = r"""\begin{table*}[t]
\centering
\caption{Error of cheap gap observables against the experimental optical gap. Donor--acceptor small molecules of
HOPV15 for which every observable is available ($n = %d$). Bias and MAE are for the raw observable; recall counts the
molecules whose experimental and computed gaps both lie in the \qtyrange{1.1}{2.2}{eV} window. Calibrated errors use a
linear fit to experiment evaluated leave-one-out by molecule (95\%% bootstrap interval, 10\,000 resamples) and by
source publication (10 groups). Acenes: mean absolute error against measured $L_a$ energies for four acenes outside
the calibration class; for $\Eks$, raw / after the HOPV15 calibration.}
\label{tab:errors}
\begin{tabular}{lccccccc}
\hline
Observable & Bias [eV] & MAE [eV] & $r$ & Recall & LOO MAE [eV] & LOPO MAE [eV] & Acenes MAE [eV]\\
\hline
%s
\hline
\end{tabular}
\end{table*}
""" % (len(m), "\n".join(rows))
    return tex, vals


def t2(root):
    d = pd.read_csv(root / "outputs/extract/gate_cascade.csv")
    a = pd.read_csv(root / "outputs/extract/applicability.csv")
    g = {"G0 validity (R1--R6)": d.SMILES.map(reject_reason) == "",
         r"G1 record check, $\Esone-\Eks \le \qty{1}{eV}$": ~((d.E_S1_eV - d["GAP(eV)"]) > 1.0),
         r"G2 calibrated $\Eopt$ in \qtyrange{1.1}{2.2}{eV}": d.eg_opt.between(1.1, 2.2),
         r"G3 $\Voc > 0$ vs PCBM": d.voc > 0,
         r"G4 LUMO offset $\ge \qty{0.3}{eV}$": d.drive >= 0.3}
    keep = np.ones(len(d), bool)
    rows, vals = [], {"records": len(d)}
    for lab, m in g.items():
        keep &= m.to_numpy()
        rows.append(f"{lab} & {int(m.sum()):,} & {int(keep.sum()):,} \\\\".replace(",", "\\,"))
        vals[lab] = [int(m.sum()), int(keep.sum())]
    assert keep.sum() == d.keep.sum()
    s = d[d.keep]
    sv = a[a.set == "cascade survivors"]
    inside = int((sv.mahal <= 4.29).sum())
    vals.update(inside=inside, pce065=[s.pce_ff065.min(), s.pce_ff065.median(), s.pce_ff065.max()],
                pcegreen=[s.pce_green.min(), s.pce_green.median(), s.pce_green.max()])
    fate = {8671: "pentacene", 7080: "tetracene", 9168: "9168", 22936: "22936 (anthraquinone dye)"}
    ftxt = "; ".join(f"{n}: {'passes' if bool(d.loc[d.mol_id == i, 'keep'].iloc[0]) else 'removed'}"
                     for i, n in fate.items())
    tex = r"""\begin{table}[t]
\centering
\caption{Calibrated cascade on %s PubChemQC records. ``Alone'': records passing that gate on its own;
``cumulative'': records passing it and every gate above. Of the %d molecules left, %d lie inside the chemical domain
of the calibration set. Scharber efficiency of the survivors (upper bound, PCBM acceptor): %.1f--%.1f\%% (median
%.1f\%%) with $\FF = 0.65$, %.1f--%.1f\%% with Green's fill factor. Reference molecules: %s.}
\label{tab:cascade}
\begin{tabular}{lrr}
\hline
Gate & Alone & Cumulative\\
\hline
%s
\hline
\end{tabular}
\end{table}
""" % (f"{len(d):,}".replace(",", "\\,"), int(keep.sum()), inside, vals["pce065"][0], vals["pce065"][2],
       vals["pce065"][1], vals["pcegreen"][0], vals["pcegreen"][2], ftxt, "\n".join(rows))
    return tex, vals


if __name__ == "__main__":
    root = Path(sys.argv[1])
    out = root / "Paper1_AdvTheorySimul/manuscript/tables"
    out.mkdir(parents=True, exist_ok=True)
    t1tex, v1 = t1(root)
    t2tex, v2 = t2(root)
    (out / "T1.tex").write_text(t1tex)
    (out / "T2.tex").write_text(t2tex)
    (out / "tables_values.json").write_text(json.dumps({"T1": v1, "T2": v2}, indent=1, default=float))
    print(json.dumps({"T1": v1, "T2": v2}, indent=1, default=lambda o: round(float(o), 4)))
