#!/usr/bin/env python3
"""
regen_figures.py — regenerate every data-driven Paper-1 figure from the
rebuilt screen (outputs/recompute/screen_full.csv).

WHY: the submitted figures were drawn from the old pipeline, whose J_SC was a
fitted Gaussian, P_in was 900.14 W m^-2 and FF used a hard-coded ideality n = 2
(ledger L-2026-09-01-05/-06/-13/-16). Every number below now comes from
scharber_recompute.py. No figure carries a number that is not in screen_full.csv,
summary.json or robustness.json.

Figures regenerated here:
  figure1_workflow                  the honest cascade, both gate conventions
  figure2_pce_sascore_scatter       landscape + the artifacts that still pass
  figure3_structures                the promoted molecules, R7 flags visible
  figure8_sensitivity               viable-set size vs SAScore weight, recomputed
  figure9_validity_recursion        NEW: what each filter layer removes (L-21)

NOT regenerated (documented reasons, not oversights):
  figure4_orbitals      molden wavefunctions exist only for 17851 and 1712;
                        there is none for 9168. Kept as the case-study pair.
  figure6_heteroatom    composition-vs-frontier correlations do not involve the
                        Scharber model, so the recomputation cannot change them.
  figure7_reorg         separate data source (reorganization_energies_corrected.csv);
                        handled by generate_fig7_reorg.py, which needs the 19531
                        counterexample added (T0-5), not a Scharber rerun.

Run from the Paper-1 folder:  python regen_figures.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

import figstyle as fs

HERE = Path(__file__).resolve().parent
REC = HERE / "outputs" / "recompute"
FIGDIR = HERE / "figures"
FIGDIR.mkdir(exist_ok=True)

df = pd.read_csv(REC / "screen_full.csv")
summary = json.loads((REC / "summary.json").read_text())
robust = json.loads((REC / "robustness.json").read_text())

VARIANTS = ["strict_FF065", "strict_FFgreen", "compl_FF065", "compl_FFgreen"]
LEAD = 9168                                    # robust across all four variants
PAH = [9168, 9104, 9139, 9406]                 # chemically defensible family
QUINONE = [19598, 21736, 23450]                # R7-flagged tautomer artifacts
DYES = [22936]                                 # anthraquinone sulfonic-acid dye
ARTIFACTS = {977: "O$_2$", 11029: "MgCO$_3$", 7801: "quinhydrone",
             24964: "CaO$_3$S$_2$", 22986: "Na$_2$ salt"}

INK = {"grey": "#475569", "blue": "#1d4ed8", "red": "#b91c1c",
       "green": "#047857", "gold": "#b45309", "purple": "#6b21a8"}
FILL = {"grey": "#f1f5f9", "blue": "#eff6ff", "red": "#fef2f2",
        "green": "#ecfdf5", "gold": "#fffbeb", "purple": "#faf5ff"}
TEXT = "#1e293b"


def row(mid):
    return df.loc[df.mol_id == mid].iloc[0]


def cascade(kind):
    """pool -> chemically valid -> valid+stable -> PCE_SAScore > 0, for one gate set."""
    pool = df[f"passes_{kind}"].astype(bool)
    ps = f"PCE_SAScore_{kind}_FF065" if False else f"PCE_SAScore_{kind}_FF065"
    valid = pool & df.chemically_valid.astype(bool)
    stable = valid & df.stable.astype(bool)
    positive = stable & (df[ps] > 0)
    return int(pool.sum()), int(valid.sum()), int(stable.sum()), int(positive.sum())


# =============================================================================
# FIGURE 1 — the honest cascade
# =============================================================================
def figure1():
    n_tot = summary["n_molecules"]
    n_win = int(df.in_visible_window.sum())
    s_pool, s_val, s_sta, s_pos = cascade("strict")
    c_pool, c_val, c_sta, c_pos = cascade("compl")
    # The strict pool is a SUBSET of the complementary pool (16 of 25, both inside
    # the 41-molecule window). 16 + 25 = 41 is a coincidence, so the figure must
    # state the nesting explicitly or a reader will read it as a partition.
    nested = int((df.passes_strict.astype(bool) & df.passes_compl.astype(bool)).sum())
    assert nested == s_pool, "strict pool is no longer a subset of the complementary pool"

    LW, RW, H = 3.85, 4.05, 1.06
    LX, RX = 0.25, 5.15
    GX = (LX + LW + RX) / 2                    # gutter centre, for the elbow route
    FSZ = 6.8
    STEP = 1.30
    L_TOP, R_TOP = 9.05, 8.40

    pipeline = [
        "PubChemQC, first shard\nIDs 1\u201324,999",
        f"with a B3LYP@PM6 record\n{n_tot:,} molecules",
        "B3LYP/6-31G*//PM6\nHOMO, LUMO, gap, dipole",
        "Scharber model, AM1.5G\n$J_{\\mathrm{SC}}$ from the real spectrum",
        "physical gates\n$E_{\\mathrm{CT}} \\leq E_g$, $V_{\\mathrm{OC}} > 0$, $\\Delta E \\geq 0.3$ eV",
        f"admissible pairs\n{c_pool} complementary, of which {s_pool} strict",
    ]
    pipe_keys = ["grey", "grey", "blue", "blue", "blue", "grey"]
    gate = [
        (f"ranked by PCE$_{{\\mathrm{{SAScore}}}}$\ntop ranks include MgCO$_3$, CaO$_3$S$_2$", "red"),
        ("chemical-validity filter\nR1\u2013R6: fragments, metals, rings,\naromaticity, size, conjugation", "green"),
        (f"chemically valid\n{c_val} complementary, {s_val} strict", "green"),
        ("stability screen\nazide / nitroso / diazo", "green"),
        (f"valid, stable, score $>0$\n{c_pos} complementary, {s_pos} strict", "gold"),
        ("R7 diagnostic (this work)\n$o$-quinone tautomers flagged\nrobust lead: molecule 9168", "purple"),
    ]

    fig, ax = plt.subplots(figsize=(fs.DOUBLE * 0.68, 14.2 / 2.54))
    ax.set_xlim(0, 9.25); ax.set_ylim(0.55, 10.55); ax.axis("off")
    fig.patch.set_facecolor("white")

    pad = 0.24
    l_bot = L_TOP - STEP * (len(pipeline) - 1)
    r_bot = R_TOP - STEP * (len(gate) - 1)
    ax.add_patch(FancyBboxPatch((LX - pad, l_bot - pad), LW + 2 * pad,
                                (L_TOP - l_bot) + H + 2 * pad,
                                boxstyle="round,pad=0,rounding_size=0.16",
                                linewidth=0, facecolor="#f8fafc", zorder=0))
    ax.add_patch(FancyBboxPatch((RX - pad, r_bot - pad), RW + 2 * pad,
                                (R_TOP - r_bot) + H + 2 * pad,
                                boxstyle="round,pad=0,rounding_size=0.16",
                                linewidth=0, facecolor="#f8fafc", zorder=0))

    def box(x, y, w, text, key, emph=False):
        ax.add_patch(FancyBboxPatch((x + 0.035, y - 0.045), w, H,
                                    boxstyle="round,pad=0.02,rounding_size=0.09",
                                    linewidth=0, facecolor="#0f172a", alpha=0.06, zorder=1))
        ax.add_patch(FancyBboxPatch((x, y), w, H,
                                    boxstyle="round,pad=0.02,rounding_size=0.09",
                                    linewidth=1.6 if emph else 1.1,
                                    edgecolor=INK[key], facecolor=FILL[key], zorder=2))
        ax.text(x + w / 2, y + H / 2, text, ha="center", va="center",
                fontsize=FSZ + (0.3 if emph else 0), color=TEXT,
                fontweight="bold" if emph else "normal", zorder=3)
        return x + w / 2, y, y + H

    def arrow(x1, y1, x2, y2, color, scale=6.5):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                     arrowstyle="-|>,head_length=2.4,head_width=1.4",
                                     mutation_scale=scale, linewidth=1.1, color=color,
                                     zorder=1.5, shrinkA=0, shrinkB=0))

    pc = [box(LX, L_TOP - i * STEP, LW, t, k)
          for i, (t, k) in enumerate(zip(pipeline, pipe_keys))]
    for a, b in zip(pc, pc[1:]):
        arrow(a[0], a[1] - 0.02, b[0], b[2] + 0.02, INK["grey"])

    gc = [box(RX, R_TOP - i * STEP, RW, t, k, emph=(k in ("gold", "purple")))
          for i, (t, k) in enumerate(gate)]
    for i, (a, b) in enumerate(zip(gc, gc[1:])):
        colour = INK["purple"] if i == len(gc) - 2 else INK["green"]
        arrow(a[0], a[1] - 0.02, b[0], b[2] + 0.02, colour)

    # hand-off: an ELBOW, not a diagonal — a diagonal reads as a drawing error
    y_from = (L_TOP - STEP * (len(pipeline) - 1)) + H / 2
    y_to = R_TOP + H / 2
    ax.plot([LX + LW, GX], [y_from, y_from], color=INK["red"], lw=1.1,
            solid_capstyle="round", zorder=1.5)
    ax.plot([GX, GX], [y_from, y_to], color=INK["red"], lw=1.1,
            solid_capstyle="round", zorder=1.5)
    arrow(GX, y_to, RX - 0.005, y_to, INK["red"], scale=6.5)
    ax.text(GX, (y_from + y_to) / 2, "ranked output", rotation=90,
            ha="center", va="center", fontsize=6.0, color=INK["red"], style="italic",
            bbox=dict(boxstyle="round,pad=0.14", facecolor="white",
                      edgecolor="none", alpha=0.95), zorder=4)

    def header(cx, y, text, color):
        ax.text(cx, y, text, fontsize=9, weight="bold", ha="center", color=color)
        ax.plot([cx - 1.25, cx + 1.25], [y - 0.22, y - 0.22], color=color,
                linewidth=1.6, solid_capstyle="round")

    header(LX + LW / 2, L_TOP + H + 0.42, "Screening pipeline", INK["grey"])
    header(RX + RW / 2, R_TOP + H + 0.42, "Validity gates (this work)", INK["green"])
    ax.text(LX + LW / 2, l_bot - 0.52,
            f"{n_win} of {n_tot:,} molecules ({100*n_win/n_tot:.2f}%) absorb in\n"
            f"1.1\u20132.2 eV; {c_pool} of those clear the physical gates,\n"
            f"and the {s_pool} strict pairs are a subset of those {c_pool}",
            ha="center", va="top", fontsize=6.2, color=INK["grey"], style="italic")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figure1_workflow.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure1_workflow      cascade strict {s_pool}->{s_val}->{s_sta}->{s_pos}, "
          f"compl {c_pool}->{c_val}->{c_sta}->{c_pos}; strict subset of compl: {nested}=={s_pool}")


# =============================================================================
# FIGURE 2 — PCE vs SAScore, with the artifacts that survive the physical gates
# =============================================================================
def figure2():
    var = "compl_FF065"          # widest admissible pool; stated in the caption
    pce, ps = f"PCE_{var}", f"PCE_SAScore_{var}"
    d = df.dropna(subset=[pce])

    fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE * 0.95))
    ax.scatter(df.SAScore, np.where(df[pce].isna(), np.nan, df[pce]),
               s=4, c="#cbd5e0", alpha=0.55, linewidths=0, rasterized=True,
               label=None)

    lim = [0, 10]
    ax.plot(lim, lim, color="#2d3748", ls="--", lw=1)

    # artifacts that pass the physical gates but fail validity/stability
    shown_art = [m for m in (11029, 24964, 22986, 4550)
                 if m in set(d.mol_id) and not (row(m).chemically_valid and row(m).stable)]
    for mid in shown_art:
        r = row(mid)
        if pd.isna(r[pce]):
            continue
        ax.scatter(r.SAScore, r[pce], marker="x", s=58, c=INK["red"], lw=1.9, zorder=5)
        name = ARTIFACTS.get(mid, str(mid))
        ax.annotate(name, (r.SAScore, r[pce]), textcoords="offset points",
                    xytext=(6, -3), fontsize=7, color=INK["red"])

    # R7-flagged quinone tautomers
    for mid in QUINONE:
        r = row(mid)
        if pd.isna(r[pce]):
            continue
        ax.scatter(r.SAScore, r[pce], marker="s", s=34, facecolor="none",
                   edgecolor=INK["purple"], lw=1.3, zorder=5)
    r = row(QUINONE[0])
    ax.annotate("$o$-quinone\ntautomers", (r.SAScore, r[pce]),
                textcoords="offset points", xytext=(-32, 8), fontsize=6.6,
                color=INK["purple"], ha="center")

    # the robust lead
    r = row(LEAD)
    ax.scatter(r.SAScore, r[pce], marker="o", s=52, facecolor=INK["green"],
               edgecolor="black", lw=0.7, zorder=6)
    ax.annotate(f"{LEAD}\nC$_{{26}}$H$_{{16}}$", (r.SAScore, r[pce]),
                textcoords="offset points", xytext=(7, 2), fontsize=7,
                color=INK["green"], weight="bold")

    ax.set_xlabel("SAScore")
    ax.set_ylabel("predicted PCE (%)")
    ax.set_xlim(0.5, 10)
    ax.set_ylim(-0.6, 15.5)
    handles = [
        Line2D([0], [0], ls="--", color="#2d3748", label="PCE$_{\\mathrm{SAScore}} = 0$"),
        Line2D([0], [0], marker="o", ls="", mfc=INK["green"], mec="black",
               label="robust lead (1)"),
        Line2D([0], [0], marker="s", ls="", mfc="none", mec=INK["purple"],
               label=f"R7-flagged tautomer ({len(QUINONE)})"),
        Line2D([0], [0], marker="x", ls="", color=INK["red"],
               label=f"passes physics, fails validity ({len(shown_art)})"),
    ]
    ax.legend(handles=handles, fontsize=6.4, loc="upper left", framealpha=0)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figure2_pce_sascore_scatter.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure2_scatter       {len(d)} molecules with a finite PCE, "
          f"{len(shown_art)} artifacts marked, lead {LEAD} at PCE {row(LEAD)[pce]:.2f} %")


# =============================================================================
# FIGURE 8 — sensitivity of the viable set to the SAScore weight, recomputed
# =============================================================================
def figure8():
    weights = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE * 0.86))

    for var, colour, marker, label in (
            ("strict_FF065", INK["blue"], "o", "strict, $FF = 0.65$"),
            ("compl_FF065", INK["gold"], "s", "complementary, $FF = 0.65$")):
        kind = var.split("_")[0]
        pool = df[f"passes_{kind}"].astype(bool)
        pce = df[f"PCE_{var}"]
        raw, kept = [], []
        for w in weights:
            score = pce - w * df.SAScore
            raw.append(int((pool & (score > 0)).sum()))
            kept.append(int((pool & df.chemically_valid.astype(bool)
                             & df.stable.astype(bool) & (score > 0)).sum()))
        ax.plot(weights, raw, marker=marker, ls="-", color=colour, lw=1.4, ms=4.5,
                label=f"{label}, before filters")
        ax.plot(weights, kept, marker=marker, ls="--", color=colour, lw=1.4, ms=4.5,
                mfc="white", label=f"{label}, after filters")

    ax.axvline(1.0, color="#94a3b8", ls=":", lw=0.9)
    ax.text(1.04, ax.get_ylim()[1] * 0.94, "Eq.~(6)", fontsize=6.4, color="#64748b")
    ax.set_xlabel("SAScore weight $w$ in PCE $-\\, w\\,\\times$ SAScore")
    ax.set_ylabel("molecules with score $> 0$")
    ax.set_xlim(0.1, 3.15)
    ax.legend(fontsize=5.9, loc="upper right", framealpha=0, ncol=1)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figure8_sensitivity.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)

    kind, var = "compl", "compl_FF065"
    pool = df[f"passes_{kind}"].astype(bool)
    out = []
    for w in weights:
        score = df[f"PCE_{var}"] - w * df.SAScore
        out.append((w, int((pool & (score > 0)).sum()),
                    int((pool & df.chemically_valid.astype(bool)
                         & df.stable.astype(bool) & (score > 0)).sum())))
    print("  figure8_sensitivity   complementary w:(raw,kept) " +
          " ".join(f"{w}:({a},{b})" for w, a, b in out))


# =============================================================================
# FIGURE 9 (new) — what each filter layer actually removes
# =============================================================================
def figure9():
    layers, counts, colours = [], [], []
    kind, var = "compl", "compl_FF065"
    pool = df[f"passes_{kind}"].astype(bool)
    pos = pool & (df[f"PCE_SAScore_{var}"] > 0)

    n0 = int(pos.sum())
    v = pos & df.chemically_valid.astype(bool)
    n1 = int(v.sum())
    s = v & df.stable.astype(bool)
    n2 = int(s.sum())
    r7 = s & ~df.mol_id.isin(QUINONE)
    n3 = int(r7.sum())
    # 22936 is correctly encoded, so R7 misses it; it is an anthraquinone
    # sulfonic-acid dye, excluded on a second, independent ground (ledger -21).
    r8 = r7 & ~df.mol_id.isin(DYES)
    n4 = int(r8.sum())

    layers = ["physical gates\n+ score $> 0$", "R1\u2013R6\nvalidity", "stability\nscreen",
              "R7 tautomer\ndiagnostic", "dye / acidic\nsolubiliser"]
    counts = [n0, n1, n2, n3, n4]
    colours = [INK["red"], INK["blue"], INK["green"], INK["purple"], INK["gold"]]

    fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE * 0.72))
    bars = ax.bar(range(len(counts)), counts, color=colours, width=0.62,
                  edgecolor="black", linewidth=0.5)
    for i, (b, c) in enumerate(zip(bars, counts)):
        ax.text(b.get_x() + b.get_width() / 2, c + max(counts) * 0.03, str(c),
                ha="center", fontsize=8, weight="bold", color=TEXT)
        if i:
            removed = counts[i - 1] - c
            ax.annotate(f"$-${removed}", (i - 0.5, max(counts) * 0.62),
                        ha="center", fontsize=6.8, color="#64748b")
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels(layers, fontsize=6.4)
    ax.set_ylabel("molecules surviving")
    ax.set_ylim(0, max(counts) * 1.22)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figure9_validity_recursion.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure9_recursion     {n0} -> {n1} -> {n2} -> {n3} -> {n4} (complementary, FF 0.65)")


# =============================================================================
# FIGURE 3 — structures of the molecules the corrected screen promotes
# =============================================================================
def figure3():
    """2D structures with the reason each molecule is or is not a candidate.

    The submitted version showed 17851, 20778, 4550 and 1712 — none of which
    survives the rebuilt screen (ledger -18, -20). This version shows the
    complementary-convention survivors plus the lead's PAH siblings, with the
    R7 tautomer flag drawn where it applies.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Draw import rdMolDraw2D
    from PIL import Image
    import io
    RDLogger.DisableLog("rdApp.*")

    var = "compl_FF065"
    panel = [
        (9168, "robust lead", "green"),
        (19598, "$o$-quinone (R7)", "purple"),
        (21736, "$o$-quinone (R7)", "purple"),
        (22936, "anthraquinone dye", "red"),
        (23450, "$o$-quinone (R7)", "purple"),
        (9104, "PAH sibling", "blue"),
    ]

    def render(smiles, size=(430, 340), highlight=False):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        Chem.rdDepictor.Compute2DCoords(mol)
        hl_atoms, hl_bonds = [], []
        if highlight:                      # mark double bonds on aromatic atoms
            for b in mol.GetBonds():
                a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
                if b.GetBondType() == Chem.BondType.DOUBLE and (
                        a1.GetIsAromatic() or a2.GetIsAromatic()):
                    hl_bonds.append(b.GetIdx())
                    hl_atoms += [a1.GetIdx(), a2.GetIdx()]
        drawer = rdMolDraw2D.MolDraw2DCairo(*size)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.bondLineWidth = 2
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer, mol, highlightAtoms=hl_atoms or None,
            highlightBonds=hl_bonds or None)
        drawer.FinishDrawing()
        return Image.open(io.BytesIO(drawer.GetDrawingText()))

    fig, axes = plt.subplots(2, 3, figsize=(fs.DOUBLE, fs.DOUBLE * 0.52))
    for ax, (mid, tag, key) in zip(axes.ravel(), panel):
        r = row(mid)
        img = render(r.SMILES, highlight=(key == "purple"))
        ax.axis("off")
        if img is not None:
            ax.imshow(img)
        pce = r[f"PCE_{var}"]
        ps = r[f"PCE_SAScore_{var}"]
        num = "n/a" if pd.isna(pce) else f"PCE {pce:.2f}%,  score {ps:+.2f}"
        ax.set_title(f"{mid}   {tag}", fontsize=7.6, color=INK[key], pad=2)
        ax.text(0.5, -0.04,
                f"{r.formula},  gap {r.gap_eV:.2f} eV,  SAScore {r.SAScore:.2f}\n{num}",
                transform=ax.transAxes, ha="center", va="top", fontsize=6.5, color=TEXT)
    fig.tight_layout(h_pad=2.6)
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figure3_structures.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure3_structures    {len(panel)} panels: lead + 3 R7-flagged + dye + PAH sibling")


if __name__ == "__main__":
    print("regenerating data-driven figures from outputs/recompute/screen_full.csv")
    figure1()
    figure2()
    figure3()
    figure8()
    figure9()
    print("done")
