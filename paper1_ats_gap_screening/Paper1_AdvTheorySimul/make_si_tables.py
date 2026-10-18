"""Supporting-Information tables, computed from the primary data only.

S-T1 validity rules, S-T2 std2 vs experiment (n = 36), S-T3 survivors of the calibrated cascade (n = 37).
usage: python -I make_si_tables.py <project_root>   -> manuscript/si/tables/ST{1,2,3}.tex + si_values.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "pipeline"))
from applicability import frame, loo_threshold, mahalanobis  # noqa: E402
from filter_genuine_osc import reject_reason  # noqa: E402
from make_tables import loo  # noqa: E402

RULES = [("R1_multifragment", "R1", "more than one fragment (salts, mixtures, co-crystals)"),
         ("R2_metal", "R2", "contains a metal"),
         ("R3_no_ring", "R3", "no ring"),
         ("R4_no_aromatic_ring", "R4", "no aromatic ring"),
         ("R5_too_small", "R5", "fewer than six heavy atoms"),
         ("R6_low_conjugation", "R6", "fewer than six conjugated atoms"),
         ("unparsable", "--", "SMILES not parsable by RDKit")]


def st1(cas):
    r = cas.SMILES.map(reject_reason)
    n = r.value_counts()
    rows = [f"{lab} & {txt} & {int(n.get(key, 0)):,} \\\\".replace(",", "\\,") for key, lab, txt in RULES]
    vals = {lab: int(n.get(key, 0)) for key, lab, _ in RULES} | {"valid": int((r == "").sum()), "records": len(cas)}
    assert sum(int(n.get(k, 0)) for k, _, _ in RULES) + vals["valid"] == len(cas)
    tex = r"""\begin{table}[h]
\centering
\caption{Validity rules of gate G0, applied in the order listed; each record is counted under the first rule it fails.
%s of %s records pass.}
\label{tab:si-validity}
\begin{tabular}{llr}
\toprule
Rule & Record removed when the structure has & Records\\
\midrule
%s
\bottomrule
\end{tabular}
\end{table}
""" % (f"{vals['valid']:,}".replace(",", "\\,"), f"{len(cas):,}".replace(",", "\\,"), "\n".join(rows))
    return tex, vals


def st2(root):
    h = root / "outputs/extract/hopv15"
    m = (pd.read_csv(h / "hopv15_records.csv")
         .merge(pd.read_csv(h / "hopv15_stda_dftgeom.csv")[["idx", "E_S1_dftgeom_eV"]], on="idx")
         .merge(pd.read_csv(h / "std2_reparsed.csv")[["idx", "S1_stda", "S1_xstda"]], on="idx"))
    m["gap_opt"] = pd.to_numeric(m.gap_opt, errors="coerce")
    m = m[(m.construction == "molecule") & ~m.complement.str.contains("TiO2", na=False)]
    obs = [("gap_b3lyp_eV", r"$\Eks$, B3LYP/def2-SVP"), ("S1_stda", "sTDA@B3LYP/6-31G*"),
           ("S1_xstda", "XsTDA@B3LYP/6-31G*"), ("E_S1_dftgeom_eV", "sTDA-xTB")]
    m = m.dropna(subset=["gap_opt"] + [c for c, _ in obs]).reset_index(drop=True)
    y = m.gap_opt.values
    win = (y >= 1.1) & (y <= 2.2)
    rows, vals = [], {"n": len(m)}
    for c, lab in obs:
        x = m[c].values
        e = x - y
        v = dict(bias=e.mean(), mae=np.abs(e).mean(), sd=e.std(ddof=1), r=np.corrcoef(x, y)[0, 1],
                 recall=f"{int((win & (x >= 1.1) & (x <= 2.2)).sum())}/{int(win.sum())}", loo=loo(x, y).mean())
        vals[c] = v
        rows.append(f"{lab} & {v['bias']:+.2f} & {v['mae']:.2f} & {v['sd']:.2f} & {v['r']:.2f} & {v['recall']} & "
                    f"{v['loo']:.3f} \\\\")
    tex = r"""\begin{table}[h]
\centering
\caption{Gap observables computed on the DFT geometries distributed with HOPV15, against the experimental optical gap
($n = %d$ donor--acceptor small molecules). SD: standard deviation of the signed error. Recall: molecules whose
experimental and computed gaps both lie in the \qtyrange{1.1}{2.2}{eV} window. LOO MAE: leave-one-out error of a linear
calibration against experiment.}
\label{tab:si-std2}
\begin{tabular}{lrrrrrr}
\toprule
Observable & Bias [eV] & MAE [eV] & SD [eV] & $r$ & Recall & LOO MAE [eV]\\
\midrule
%s
\bottomrule
\end{tabular}
\end{table}
""" % (len(m), "\n".join(rows))
    return tex, vals


def st3(root, cas):
    h = root / "outputs/extract/hopv15"
    hp = pd.read_csv(h / "hopv15_records.csv").merge(pd.read_csv(h / "hopv15_stda.csv"), on="idx")
    hp["gap_opt"] = pd.to_numeric(hp.gap_opt, errors="coerce")
    hp = hp[(hp.construction == "molecule") & ~hp.complement.str.contains("TiO2", na=False)
            & (hp.status == "success")].dropna(subset=["gap_opt"]).reset_index(drop=True)
    ref = frame(hp.smiles.tolist()).dropna().reset_index(drop=True)
    thr = loo_threshold(ref)
    s = cas[cas.keep].copy().sort_values("pce_ff065", ascending=False).reset_index(drop=True)
    s["mahal"] = mahalanobis(ref, frame(s.SMILES.tolist()))
    s["inside"] = s.mahal <= thr
    s.to_csv(root / "outputs/extract/si_survivors.csv", index=False)
    rows = []
    for _, r in s.iterrows():
        rows.append(f"{int(r.mol_id)} & {r.formula} & \\texttt{{\\seqsplit{{{r.SMILES}}}}} & {r['GAP(eV)']:.2f} & "
                    f"{r.eg_opt:.2f} & {r.E_S1_eV:.2f} & {r.voc:.2f} & {r.pce_ff065:.1f} & {r.pce_green:.1f} & "
                    f"{r.mahal:.1f}{'' if r.inside else '$^{*}$'} \\\\")
    vals = {"n": len(s), "inside": int(s.inside.sum()), "threshold": thr}
    tex = r"""{\footnotesize
\begin{longtable}{rlp{4.1cm}rrrrrrr}
\caption{The %d molecules that pass every gate of the calibrated cascade, ordered by Scharber efficiency (upper bound,
PCBM acceptor). $\Eks$: stored PubChemQC gap; $\Eopt$: calibrated optical gap; $\Esone$: sTDA-xTB singlet energy;
PCE$_{0.65}$ and PCE$_\mathrm{G}$: efficiency with $\FF = 0.65$ and with Green's fill factor; $d_\mathrm{M}$:
Mahalanobis distance to the calibration set, $^{*}$ marking the %d molecules beyond the domain threshold
($d_\mathrm{M} > %.2f$). Energies in eV, $\Voc$ in V, efficiencies in \%%.}\label{tab:si-survivors}\\
\toprule
ID & Formula & SMILES & $\Eks$ & $\Eopt$ & $\Esone$ & $\Voc$ & PCE$_{0.65}$ & PCE$_\mathrm{G}$ & $d_\mathrm{M}$\\
\midrule
\endfirsthead
\toprule
ID & Formula & SMILES & $\Eks$ & $\Eopt$ & $\Esone$ & $\Voc$ & PCE$_{0.65}$ & PCE$_\mathrm{G}$ & $d_\mathrm{M}$\\
\midrule
\endhead
%s
\bottomrule
\end{longtable}}
""" % (len(s), len(s) - vals["inside"], thr, "\n".join(rows))
    return tex, vals


if __name__ == "__main__":
    root = Path(sys.argv[1])
    out = root / "Paper1_AdvTheorySimul/manuscript/si/tables"
    out.mkdir(parents=True, exist_ok=True)
    cas = pd.read_csv(root / "outputs/extract/gate_cascade.csv")
    res = {}
    for name, (tex, v) in {"ST1": st1(cas), "ST2": st2(root), "ST3": st3(root, cas)}.items():
        (out / f"{name}.tex").write_text(tex)
        res[name] = v
    (out / "si_values.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=lambda o: round(float(o), 4)))
