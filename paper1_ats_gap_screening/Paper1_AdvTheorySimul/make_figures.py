"""Figures F2–F5 of the ATS manuscript, built only from the primary data files listed in PROVENANCE.

usage: python -I make_figures.py <project_root>
writes Paper1_AdvTheorySimul/manuscript/figures/F{2,3,4,5}.{pdf,png} + figures_provenance.json
Okabe–Ito colours, every colour doubled by marker shape; dense scatters rasterised inside vector PDFs.
"""
import json
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OK = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
      "yellow": "#F0E442", "blue": "#0072B2", "verm": "#D55E00", "purple": "#CC79A7", "grey": "#999999"}
MM = 1 / 25.4
SINGLE, DOUBLE = 85 * MM, 178 * MM
WIN = (1.1, 2.2)
ACENES = {931: ("naphthalene", 4.445), 8418: ("anthracene", 3.36), 7080: ("tetracene", 2.66), 8671: ("pentacene", 2.215)}
mpl.rcParams.update({"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7, "legend.fontsize": 6,
                     "xtick.labelsize": 6, "ytick.labelsize": 6, "pdf.fonttype": 42, "axes.linewidth": 0.6,
                     "savefig.facecolor": "white"})


def save(fig, out, name):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}", dpi=600)
    plt.close(fig)


def window(ax, axis="both"):
    kw = dict(color=OK["grey"], alpha=0.15, lw=0, zorder=0)
    if axis in ("x", "both"):
        ax.axvspan(*WIN, **kw)
    if axis in ("y", "both"):
        ax.axhspan(*WIN, **kw)


def hopv(root):
    h = root / "outputs/extract/hopv15"
    m = (pd.read_csv(h / "hopv15_records.csv")
         .merge(pd.read_csv(h / "hopv15_stda.csv")[["idx", "E_S1_eV"]], on="idx", how="left")
         .merge(pd.read_csv(h / "hopv15_stda_dftgeom.csv")[["idx", "E_S1_dftgeom_eV"]], on="idx", how="left")
         .merge(pd.read_csv(h / "std2_reparsed.csv")[["idx", "S1_stda"]], on="idx", how="left"))
    m["gap_opt"] = pd.to_numeric(m.gap_opt, errors="coerce")
    return m[(m.construction == "molecule") & ~m.complement.str.contains("TiO2", na=False)].dropna(subset=["gap_opt"])


def f2(root, out):
    m = pd.read_csv(root / "outputs/extract/ks_vs_stda_joined.csv")
    d = pd.read_csv(root / "ZENODO_UPLOAD/dataset_pubchemqc_opv_17458.csv", usecols=["mol_id", "SMILES"], dtype=str)
    d["mol_id"] = pd.to_numeric(d.mol_id)
    m = m.merge(d, on="mol_id")
    salt = m.SMILES.str.contains(".", regex=False)
    flag = (m.E_S1_eV - m["GAP(eV)"] > 1.0) & ~salt
    fig, ax = plt.subplots(figsize=(SINGLE, SINGLE), layout="constrained")
    window(ax)
    base = ~salt & ~flag
    ax.scatter(m.loc[base, "GAP(eV)"], m.loc[base, "E_S1_eV"], s=1.5, c=OK["grey"], marker="o", lw=0,
               alpha=0.5, rasterized=True, label=f"single fragment (n = {base.sum():,})")
    ax.scatter(m.loc[salt, "GAP(eV)"], m.loc[salt, "E_S1_eV"], s=3, c=OK["blue"], marker="s", lw=0,
               alpha=0.7, rasterized=True, label=f"multi-fragment (n = {salt.sum():,})")
    ax.scatter(m.loc[flag, "GAP(eV)"], m.loc[flag, "E_S1_eV"], s=12, facecolors="none", edgecolors=OK["verm"],
               marker="^", lw=0.7, label=rf"flag: $E_{{S_1}}-E_{{KS}}>1$ eV (n = {flag.sum()})")
    for i, (n, _) in ACENES.items():
        r = m[m.mol_id == i].iloc[0]
        ax.scatter(r["GAP(eV)"], r.E_S1_eV, s=14, c=OK["green"], marker="D", lw=0, zorder=4)
        ax.annotate(n, (r["GAP(eV)"], r.E_S1_eV), xytext=(3, -6), textcoords="offset points", fontsize=5.5)
    r = m[m.mol_id == 9168].iloc[0]
    ax.annotate("9168", (r["GAP(eV)"], r.E_S1_eV), xytext=(-18, 8), textcoords="offset points", fontsize=6,
                color=OK["verm"], arrowprops=dict(arrowstyle="-", lw=0.5, color=OK["verm"]))
    lim = (0, 9)
    ax.plot(lim, lim, ls="--", lw=0.6, c=OK["black"])
    ax.set(xlim=lim, ylim=lim, xlabel=r"$E_\mathrm{KS}$ (PubChemQC B3LYP/6-31G*) [eV]",
           ylabel=r"$E_{S_1}$ (sTDA-xTB) [eV]", aspect="equal")
    ax.legend(loc="upper left", frameon=False, markerscale=1.5, handletextpad=0.3)
    save(fig, out, "F2")


def f3(root, out):
    m = hopv(root)
    panels = [("gap_b3lyp_eV", "KS gap (B3LYP/def2-SVP)", OK["blue"], "o"),
              ("E_S1_eV", "sTDA-xTB, UFF geometry", OK["orange"], "s"),
              ("E_S1_dftgeom_eV", "sTDA-xTB, DFT geometry", OK["purple"], "^"),
              ("S1_stda", "sTDA@B3LYP (std2)", OK["green"], "D")]
    fig, axs = plt.subplots(1, 5, figsize=(DOUBLE, DOUBLE / 4.3), layout="constrained")
    lim = (1.2, 3.9)
    for ax, (c, lab, col, mk) in zip(axs, panels):
        a = m.dropna(subset=[c])
        window(ax, "x")
        ax.plot(lim, lim, ls="--", lw=0.6, c=OK["black"])
        ax.scatter(a.gap_opt, a[c], s=9, c=col, marker=mk, lw=0)
        e = a[c] - a.gap_opt
        ax.set(xlim=lim, ylim=lim, aspect="equal", title=lab, xlabel=r"$E_\mathrm{opt}$ exp. [eV]")
        ax.text(0.96, 0.04, f"n = {len(a)}\nbias {e.mean():+.2f} eV\nMAE {e.abs().mean():.2f} eV",
                transform=ax.transAxes, va="bottom", ha="right", fontsize=5.5)
    axs[0].set_ylabel("computed gap [eV]")
    fit = m.dropna(subset=["gap_b3lyp_eV", "E_S1_eV"])  # same 34 complete cases as the cascade calibration
    k, b = np.polyfit(fit.gap_b3lyp_eV, fit.gap_opt, 1)
    j = pd.read_csv(root / "outputs/extract/ks_vs_stda_joined.csv").set_index("mol_id")
    ax = axs[4]
    lim = (1.7, 4.8)
    ax.plot(lim, lim, ls="--", lw=0.6, c=OK["black"])
    window(ax, "x")
    for i, (n, ex) in ACENES.items():
        ks, s1 = j.loc[i, "GAP(eV)"], j.loc[i, "E_S1_eV"]
        ax.scatter(ex, ks, s=11, c=OK["blue"], marker="o", lw=0)
        ax.scatter(ex, k * ks + b, s=11, facecolors="none", edgecolors=OK["blue"], marker="o", lw=0.7)
        ax.scatter(ex, s1, s=11, c=OK["orange"], marker="s", lw=0)
        ax.annotate(n, (ex, min(ks, k * ks + b, s1)), xytext=(4, -3), textcoords="offset points", fontsize=5,
                    ha="left", va="top")
    ax.scatter([], [], c=OK["blue"], marker="o", s=9, label="KS raw")
    ax.scatter([], [], facecolors="none", edgecolors=OK["blue"], marker="o", s=9, label="KS calibrated")
    ax.scatter([], [], c=OK["orange"], marker="s", s=9, label="sTDA-xTB raw")
    ax.set(xlim=lim, ylim=lim, aspect="equal", title="acenes (external class)", xlabel=r"$E_\mathrm{opt}$ exp. [eV]")
    ax.legend(loc="upper left", frameon=False, fontsize=5, handletextpad=0.2, borderaxespad=0.2)
    save(fig, out, "F3")


def f4(root, out):
    s = pd.read_csv(root / "outputs/extract/record-validity-result.csv")
    fig, ax = plt.subplots(figsize=(SINGLE, SINGLE * 0.75), layout="constrained")
    rng = np.random.default_rng(7)
    for x, (g, col, mk, lab) in enumerate([("control", OK["blue"], "o", "controls"), ("flag", OK["verm"], "^", "flagged")]):
        a = s[s.group == g]
        ax.boxplot(a.dev, positions=[x], widths=0.45, showfliers=False,
                   medianprops=dict(color=OK["black"], lw=0.8), boxprops=dict(lw=0.6),
                   whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        ax.scatter(x + rng.uniform(-0.12, 0.12, len(a)), a.dev, s=10, c=col, marker=mk, lw=0, zorder=3,
                   label=f"{lab} (n = {len(a)})")
    for i, txt, dy in ((9168, "9168", 0), (19598, "19598 (tautomer)", -7), (21736, "21736 (tautomer)", 7)):
        r = s[s.mol_id == i].iloc[0]
        ax.annotate(txt, (1.15, r.dev), xytext=(10, dy), textcoords="offset points", fontsize=5.5, va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.4))
    ax.axhline(0, lw=0.6, c=OK["black"])
    ax.axhline(0.5, lw=0.6, ls=":", c=OK["black"])
    ax.set(xticks=[0, 1], xticklabels=["controls\n" r"$-1.5\leq\Delta\leq0.5$ eV", "flagged\n" r"$\Delta>1$ eV"],
           ylabel=r"$E_\mathrm{KS}$ recomputed $-$ stored [eV]", xlim=(-0.6, 1.9))
    ax.legend(loc="upper left", frameon=False)
    save(fig, out, "F4")


def f5(root, out):
    a = pd.read_csv(root / "outputs/extract/applicability.csv")
    thr = 4.29
    fig, ax = plt.subplots(figsize=(SINGLE, SINGLE * 0.8), layout="constrained")
    rng = np.random.default_rng(1)

    def jit(n):  # integer ring counts; jitter disclosed in caption
        return rng.uniform(-0.2, 0.2, n)

    sh = a[a.set == "valid shard (sample 2000)"]
    ax.scatter(sh.heavy_atoms, sh.arom_rings + jit(len(sh)), s=2, c=OK["grey"], marker="o", lw=0, alpha=0.5,
               rasterized=True, label="shard sample (n = 2,000)")
    hp = a[a.set == "HOPV calibration"]
    ax.scatter(hp.heavy_atoms, hp.arom_rings + jit(len(hp)), s=12, c=OK["blue"], marker="o", lw=0,
               label=f"HOPV calibration (n = {len(hp)})")
    sv = a[a.set == "cascade survivors"]
    ins = sv.mahal <= thr
    ax.scatter(sv.loc[ins, "heavy_atoms"], sv.loc[ins, "arom_rings"] + jit(ins.sum()), s=16, c=OK["green"],
               marker="D", lw=0, label=f"survivors, inside domain (n = {ins.sum()})")
    ax.scatter(sv.loc[~ins, "heavy_atoms"], sv.loc[~ins, "arom_rings"] + jit((~ins).sum()), s=18,
               facecolors="none", edgecolors=OK["verm"], marker="^", lw=0.7,
               label=f"survivors, outside domain (n = {(~ins).sum()})")
    ax.set(xlabel="heavy atoms", ylabel="aromatic rings (jittered ±0.2)")
    ax.legend(loc="lower right", frameon=False)
    save(fig, out, "F5")


def f1(root, out):
    """Cascade schematic; every count recomputed from gate_cascade.csv / applicability.csv."""
    import sys as _s
    _s.path.insert(0, str(root / "Paper1_AdvTheorySimul/pipeline"))
    from filter_genuine_osc import reject_reason
    d = pd.read_csv(root / "outputs/extract/gate_cascade.csv")
    g0 = d.SMILES.map(reject_reason) == ""
    g1 = g0 & ~((d.E_S1_eV - d["GAP(eV)"]) > 1.0)
    g2 = g1 & d.eg_opt.between(*WIN)
    g3 = g2 & (d.voc > 0)
    g4 = g3 & (d.drive >= 0.3)
    assert g4.sum() == d.keep.sum(), "cascade counts drifted from gate_cascade.csv"
    a = pd.read_csv(root / "outputs/extract/applicability.csv")
    sv = a[a.set == "cascade survivors"]
    inside = int((sv.mahal <= 4.29).sum())
    steps = [(f"{len(d):,} PubChemQC records", None),
             ("G0  chemical validity (R1–R6)", int(g0.sum())),
             (r"G1  record check: drop $E_{S_1}-E_\mathrm{KS}>1$ eV", int(g1.sum())),
             (r"G2  $E_\mathrm{opt}$ = calibrated $E_\mathrm{KS}$ in 1.1–2.2 eV", int(g2.sum())),
             (r"G3  $V_\mathrm{OC}>0$ (vs PCBM)", int(g3.sum())),
             ("G4  LUMO offset ≥ 0.3 eV", int(g4.sum()))]
    fig, ax = plt.subplots(figsize=(SINGLE, SINGLE * 1.05), layout="constrained")
    ax.set_axis_off()
    y = 0.95
    for i, (lab, n) in enumerate(steps):
        txt = lab if n is None else f"{lab}\n→ {n:,} remain"
        ax.text(0.5, y, txt, ha="center", va="top", fontsize=6.5, transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.4", fc="white" if i else "#EEEEEE", ec=OK["black"], lw=0.6))
        if i < len(steps) - 1:
            ax.annotate("", xy=(0.5, y - 0.115), xytext=(0.5, y - 0.075), xycoords="axes fraction",
                        arrowprops=dict(arrowstyle="-|>", lw=0.7, color=OK["black"]))
        y -= 0.135
    ax.text(0.27, y, f"inside HOPV domain\n{inside} assessable", ha="center", va="top", fontsize=6.5,
            transform=ax.transAxes, bbox=dict(boxstyle="round,pad=0.4", fc=OK["green"], alpha=0.25, ec=OK["green"]))
    ax.text(0.73, y, f"outside the calibrated class\n{int(g4.sum()) - inside} not assessable", ha="center",
            va="top", fontsize=6.5, transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.4", fc=OK["verm"], alpha=0.2, ec=OK["verm"]))
    save(fig, out, "F1")


if __name__ == "__main__":
    root = Path(sys.argv[1])
    out = root / "Paper1_AdvTheorySimul/manuscript/figures"
    out.mkdir(parents=True, exist_ok=True)
    for f in (f1, f2, f3, f4, f5):
        f(root, out)
    prov = {
        "F1": ["outputs/extract/gate_cascade.csv (counts recomputed, asserted equal to keep)",
               "outputs/extract/applicability.csv (Mahalanobis threshold 4.29)", "ledger -11, -16, -17"],
        "F2": ["outputs/extract/ks_vs_stda_joined.csv", "ZENODO_UPLOAD/dataset_pubchemqc_opv_17458.csv",
               "ledger L-2026-10-17-02, -05, -09"],
        "F3": ["outputs/extract/hopv15/{hopv15_records,hopv15_stda,hopv15_stda_dftgeom,std2_reparsed}.csv",
               "outputs/extract/ks_vs_stda_joined.csv", "Malloci 2007 exp. La midpoints",
               "ledger -06, -14, -18; calibration on the 34 complete cases"],
        "F4": ["outputs/extract/record-validity-result.csv", "ledger -09"],
        "F5": ["outputs/extract/applicability.csv (Mahalanobis threshold 4.29)", "ledger -17; ring jitter ±0.2, seed 1"],
    }
    (out / "figures_provenance.json").write_text(json.dumps(prov, indent=2))
    print("figures written:", sorted(p.name for p in out.glob("F*.pdf")))
