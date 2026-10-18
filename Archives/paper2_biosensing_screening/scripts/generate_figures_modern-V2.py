#!/usr/bin/env python3
"""
Modern-minimal redesign of Paper 2 data figures (Fig 3, 5, 6, 7).
Rewritten 2026-07-17 at author request: same real data, new publication-modern
visual design. Numbers are read live from the same CSVs the previous versions
used -- NO value is hardcoded except the seven already-cited dipole magnitudes.

  Fig 3  docking affinities   -> per-molecule dot/range plot, shaded inhibitor band
  Fig 5  solvatochromic shift -> toluene->water dumbbell (slope), two functionals
  Fig 6  synthesis            -> two clean panels vs dipole; marker size = osc. strength
  Fig 7  structures           -> high-quality RDKit render on class-colored cards

Aesthetic: modern-minimal. Colorblind-safe two-hue scheme
(blue family = 4 stable donors, orange family = 3 reactive structures).
Data provenance unchanged from generate_figures.py / generate_figures2.py.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io

BASE = Path(__file__).resolve().parent
FIGDIR = BASE / "figures"

# ----------------------------------------------------------------------------
# Shared data (identical provenance to the previous scripts)
# ----------------------------------------------------------------------------
STABLE = ["1712", "9168", "17574", "18506"]
REACTIVE = ["4550", "17851", "20778"]
MOLECULES = STABLE + REACTIVE
DIPOLE = {  # Debye, dataset_pubchemqc_opv_17458.csv col dipole_moment
    "1712": 2.56, "9168": 3.18, "17574": 5.18, "18506": 1.60,
    "4550": 9.74, "17851": 9.26, "20778": 7.42,
}
TARGETS = ["HIV1protease", "Hsp90", "Neurodegenerative_1SYH", "COVID19_6Y2F"]
TARGET_LABELS = {
    "HIV1protease": "HIV-1 protease",
    "Hsp90": "Hsp90",
    "Neurodegenerative_1SYH": "Neurodeg. (1SYH)",
    "COVID19_6Y2F": "SARS-CoV-2 M$^{pro}$",
}
SMILES = {
    "1712": "O=C1N=c2c(=C1CNc1ccc(cc1)S(=O)(=O)Nc1ccccn1)c1scnc1cc2",
    "9168": "c1ccc2c(c1)c1ccc3c(c1cc2)ccc1c3ccc2c1cccc2",
    "17574": "CC(=O)N(c1ccc(cc1)Nc1cc(c(c2c1C(=O)c1ccccc1C2=O)N)S(=O)(=O)O)C",
    "18506": "CCCCn1c(=O)c2c(c1=O)c(N)c1c(c2N)c(=O)c2c(c1=O)cccc2",
    "4550": "O=NC1=c2cc(N(O)O)c3c(c2=NC1=O)CCCC3",
    "17851": "CC[N-]c1nc(N[N+]#N)nc(n1)NC(C)C",
    "20778": "N#[N+]Nc1nc([N-]C(C)C)nc(n1)SC",
}

# Colorblind-safe two-hue scheme: blues = stable, oranges/red = reactive.
CLR = {
    "1712":  "#0072B2",  # blue
    "9168":  "#009E73",  # teal-green
    "17574": "#56B4E9",  # light blue
    "18506": "#1b4965",  # deep navy
    "4550":  "#D55E00",  # vermillion
    "17851": "#E69F00",  # amber
    "20778": "#B30000",  # deep red
}
STABLE_INK = "#0072B2"
REACTIVE_INK = "#D55E00"
GRID = "#e3e7ee"
INKTXT = "#1f2933"

# Elsevier/ACS column widths (inches)
SINGLE = 9.0 / 2.54
DOUBLE = 19.0 / 2.54

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "axes.labelcolor": INKTXT,
    "text.color": INKTXT,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.color": INKTXT,
    "ytick.color": INKTXT,
    "legend.fontsize": 7.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#8a94a6",
    "lines.linewidth": 1.4,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
})


# ----------------------------------------------------------------------------
# Shared save helper
# ----------------------------------------------------------------------------
def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"{stem}.{ext}")
    plt.close(fig)
    print(f"wrote figures/{stem}.pdf / .png")


# ----------------------------------------------------------------------------
# Fig 2 -- workflow, modern-minimal flowchart matching the fig 3-7 palette
# ----------------------------------------------------------------------------
def _card(ax, cx, cy, w, h, lines, edge, fill, bold_first=True, fs_pt=8.0):
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.10",
        linewidth=1.3, edgecolor=edge, facecolor=fill, zorder=2))
    if isinstance(lines, str):
        lines = [lines]
    n = len(lines)
    for i, ln in enumerate(lines):
        yy = cy + (n - 1) / 2 * 0.34 - i * 0.34
        weight = "bold" if (i == 0 and bold_first) else "normal"
        col = edge if (i == 0 and bold_first) else INKTXT
        ax.text(cx, yy, ln, ha="center", va="center", fontsize=fs_pt,
                fontweight=weight, color=col, zorder=3)


def _flow_arrow(ax, x0, y0, x1, y1, color):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13,
        linewidth=1.6, color=color, zorder=1,
        shrinkA=2, shrinkB=2))


def fig2_workflow():
    blue, orange, teal = STABLE_INK, REACTIVE_INK, "#009E73"
    bluef, orangef, tealf = "#eef5fb", "#fdf3ec", "#e9f7f1"
    greyf = "#f1f4f8"

    fig, ax = plt.subplots(figsize=(DOUBLE * 0.78, DOUBLE * 0.70))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    # Column centres for the two branches. With a card width of 4.2 the left
    # card spans x=0.5..4.7 and the right card spans x=5.3..9.5, leaving a clean
    # 0.6-unit gutter between them so the branches never overlap.
    LX, RX = 2.6, 7.4
    BRANCH_W = 4.2

    # top: shared provenance
    _card(ax, 5.0, 9.1, 8.6, 1.25,
          ["Companion screen: 17,458 molecules",
           "chemical-validity + stability gate  →  4 admitted + 3 stable donors = 7"],
          "#4a5568", greyf, fs_pt=8.0)

    # fork from the shared node into the two parallel branches
    _flow_arrow(ax, 4.4, 8.45, LX, 7.15, blue)
    _flow_arrow(ax, 5.6, 8.45, RX, 7.15, orange)

    _card(ax, LX, 6.35, BRANCH_W, 1.45,
          ["Molecular recognition",
           "blind docking (smina, seeded)",
           "4 disease targets + calibration"], blue, bluef, fs_pt=7.6)
    _card(ax, RX, 6.35, BRANCH_W, 1.45,
          ["Optical sensitivity",
           "TD-DFT solvatochromic shift",
           "CPCM toluene/water, 2 functionals"], orange, orangef, fs_pt=7.6)

    _flow_arrow(ax, LX, 5.6, LX, 4.35, blue)
    _flow_arrow(ax, RX, 5.6, RX, 4.35, orange)

    _card(ax, LX, 3.6, BRANCH_W, 1.35,
          ["Result", "4 stable donors bind", "within inhibitor band"], blue, bluef, fs_pt=7.6)
    _card(ax, RX, 3.6, BRANCH_W, 1.35,
          ["Result", "large shift only in the", "reactive set (dark states)"], orange, orangef, fs_pt=7.6)

    # merge both branches back into the shared takeaway node
    _flow_arrow(ax, LX, 2.85, 4.4, 1.75, "#4a5568")
    _flow_arrow(ax, RX, 2.85, 5.6, 1.75, "#4a5568")

    _card(ax, 5.0, 1.05, 7.0, 1.25,
          ["Two independent readouts select different subsets",
           "reported side by side, not merged into one score"], teal, tealf, fs_pt=8.0)

    save(fig, "fig2_workflow")


def load_docking():
    path = BASE / "docking_rerun" / "docking_results_clean.csv"
    data = {mol: {} for mol in MOLECULES}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if row["mol_id"] in data:
                data[row["mol_id"]][row["target"]] = float(row["affinity_kcal_mol"])
    return data


def load_calibration_band():
    path = BASE / "docking_rerun" / "calibration_reference_ligands.csv"
    refs = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            refs.append(float(row["ref_affinity_kcal_mol"]))
    return min(refs), max(refs)  # e.g. (-7.7, -5.5)


def load_s1(functional):
    """Return {mol: {'Toluene': (E, f), 'Water': (E, f)}} for state 1."""
    path = BASE / "solvatochromic_tddft" / f"solvatochromic_results_{functional}.csv"
    s1 = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if int(row["state"]) == 1 and row["mol_id"] in DIPOLE:
                s1.setdefault(row["mol_id"], {})[row["solvent"]] = (
                    float(row["energy_eV"]), float(row["oscillator_strength"]))
    return s1


# ----------------------------------------------------------------------------
# Fig 3 -- docking affinities: molecules on y, 4 target dots each, inhibitor band
# ----------------------------------------------------------------------------
def fig3_docking():
    data = load_docking()
    lo, hi = load_calibration_band()  # (-7.7, -5.5)

    order = MOLECULES[::-1]  # reactive at bottom, stable at top when plotted upward
    ypos = {m: i for i, m in enumerate(order)}

    fig, ax = plt.subplots(figsize=(DOUBLE * 0.72, 0.62 * DOUBLE * 0.72))

    # reference-inhibitor band
    ax.axvspan(lo, hi, color="#c9d6e5", alpha=0.45, zorder=0, lw=0)
    ax.text(hi, len(order) - 0.35, "  native-inhibitor range",
            va="center", ha="left", fontsize=7, color="#4a5568", style="italic")

    tmark = ["o", "s", "^", "D"]
    for m in order:
        y = ypos[m]
        vals = [data[m][t] for t in TARGETS]
        ax.plot([min(vals), max(vals)], [y, y], color=CLR[m], lw=2.2,
                alpha=0.35, solid_capstyle="round", zorder=2)
        for t, mk in zip(TARGETS, tmark):
            ax.scatter(data[m][t], y, marker=mk, s=34, color=CLR[m],
                       edgecolor="white", linewidth=0.6, zorder=3)

    # class separator between stable (top 4) and reactive (bottom 3)
    sep = len(REACTIVE) - 0.5
    ax.axhline(sep, color="#cbd2dc", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(0.015, 0.83, "stable donors", transform=ax.transAxes, fontsize=7.5,
            color=STABLE_INK, fontweight="bold", va="center", ha="left")
    ax.text(0.015, 0.16, "reactive", transform=ax.transAxes, fontsize=7.5,
            color=REACTIVE_INK, fontweight="bold", va="center", ha="left")

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    for tick, m in zip(ax.get_yticklabels(), order):
        tick.set_color(CLR[m]); tick.set_fontweight("bold")
    ax.set_xlabel("Docking affinity (kcal mol$^{-1}$)  —  stronger binding →")
    ax.invert_xaxis()  # more negative (stronger) to the right
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)

    tgt_handles = [Line2D([0], [0], marker=mk, color="#4a5568", ls="none",
                          ms=6, mfc="#4a5568", label=TARGET_LABELS[t])
                   for t, mk in zip(TARGETS, tmark)]
    ax.legend(handles=tgt_handles, title="Target", loc="upper left",
              bbox_to_anchor=(0.0, -0.16), ncol=4, columnspacing=1.2,
              handletextpad=0.3, title_fontsize=7.5)
    save(fig, "fig3_docking_affinities")


# ----------------------------------------------------------------------------
# Fig 5 -- solvatochromic dumbbell: toluene -> water, per molecule, 2 functionals
# ----------------------------------------------------------------------------
def fig5_solvatochromic():
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE * 0.82, 0.5 * DOUBLE * 0.82),
                             sharey=True)
    order = MOLECULES[::-1]
    for ax, functional, title in zip(
            axes, ["B3LYP", "CAM-B3LYP"], ["B3LYP/6-31G*", "CAM-B3LYP/6-31G*"]):
        s1 = load_s1(functional)
        for i, m in enumerate(order):
            tol = s1[m]["Toluene"][0]
            wat = s1[m]["Water"][0]
            ax.plot([tol, wat], [i, i], color=CLR[m], lw=2.4, alpha=0.55,
                    solid_capstyle="round", zorder=2)
            ax.scatter(tol, i, s=32, facecolor="white", edgecolor=CLR[m],
                       linewidth=1.4, zorder=3)
            ax.scatter(wat, i, s=40, color=CLR[m], edgecolor="white",
                       linewidth=0.6, zorder=4)
            d = wat - tol
            if abs(d) > 0.05:
                ax.annotate(f"{d:+.2f}", (max(tol, wat), i), xytext=(5, 0),
                            textcoords="offset points", va="center", ha="left",
                            fontsize=6.5, color=CLR[m], fontweight="bold")
        sep = len(REACTIVE) - 0.5
        ax.axhline(sep, color="#cbd2dc", lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.set_title(title, fontsize=8.5)
        ax.set_xlabel("S$_1$ energy (eV)")
        ax.grid(axis="x", color=GRID, lw=0.7)
        ax.set_axisbelow(True)

    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels(order)
    for tick, m in zip(axes[0].get_yticklabels(), order):
        tick.set_color(CLR[m]); tick.set_fontweight("bold")
    axes[0].set_ylim(-0.6, len(order) - 0.4)

    handles = [
        Line2D([0], [0], marker="o", ls="none", ms=6, mfc="white",
               mec="#4a5568", mew=1.4, label="Toluene (nonpolar)"),
        Line2D([0], [0], marker="o", ls="none", ms=6, mfc="#4a5568",
               mec="white", label="Water (polar)"),
    ]
    axes[1].legend(handles=handles, loc="upper right", handletextpad=0.3)
    fig.subplots_adjust(bottom=0.20)
    fig.text(0.5, 0.045,
             "line length = solvatochromic shift; stable donors (top) barely move",
             ha="center", fontsize=7, color="#4a5568", style="italic")
    save(fig, "fig5_solvatochromic_shift")


# ----------------------------------------------------------------------------
# Fig 6 -- synthesis: two panels vs dipole; marker area = S1 oscillator strength
# ----------------------------------------------------------------------------
def fig6_synthesis():
    docking = load_docking()
    best = {m: min(docking[m].values()) for m in MOLECULES}
    lo, hi = load_calibration_band()
    cam = load_s1("CAM-B3LYP")

    def shift_mag(m):
        return abs(cam[m]["Water"][0] - cam[m]["Toluene"][0])

    def osc(m):
        return cam[m]["Toluene"][1]

    # marker area scaled from oscillator strength (sqrt, floored so dark dots show)
    def msize(m):
        return 28 + 620 * np.sqrt(max(osc(m), 0.0))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(DOUBLE * 0.86, 0.46 * DOUBLE * 0.86))

    # ---- Panel A: dipole vs binding ----
    axL.axhspan(lo, hi, color="#c9d6e5", alpha=0.45, lw=0, zorder=0)
    for m in MOLECULES:
        axL.scatter(DIPOLE[m], best[m], s=msize(m), color=CLR[m],
                    edgecolor="white", linewidth=0.7, zorder=3, alpha=0.95)
        axL.annotate(m, (DIPOLE[m], best[m]), xytext=(0, 9),
                     textcoords="offset points", ha="center", fontsize=6.5,
                     color=CLR[m], fontweight="bold")
    axL.set_xlabel("Dipole moment (D)")
    axL.set_ylabel("Best docking affinity (kcal mol$^{-1}$)")
    axL.invert_yaxis()
    axL.margins(y=0.16)
    axL.text(0.98, 0.5 * (lo + hi), "native-inhibitor band", fontsize=6.8,
             color="#4a5568", va="center", ha="right", style="italic",
             transform=axL.get_yaxis_transform())
    axL.grid(color=GRID, lw=0.7); axL.set_axisbelow(True)
    axL.set_title("Molecular recognition", fontsize=8.5)

    # ---- Panel B: dipole vs |shift| ----
    for m in MOLECULES:
        axR.scatter(DIPOLE[m], shift_mag(m), s=msize(m), color=CLR[m],
                    edgecolor="white", linewidth=0.7, zorder=3, alpha=0.95)
        axR.annotate(m, (DIPOLE[m], shift_mag(m)), xytext=(0, 9),
                     textcoords="offset points", ha="center", fontsize=6.5,
                     color=CLR[m], fontweight="bold")
    axR.set_xlabel("Dipole moment (D)")
    axR.set_ylabel("|Solvatochromic shift| (eV, CAM-B3LYP)")
    axR.margins(y=0.16)
    axR.grid(color=GRID, lw=0.7); axR.set_axisbelow(True)
    axR.set_title("Optical environment-sensitivity", fontsize=8.5)

    # size legend (oscillator strength) — placed lower-right where panel B is empty
    for f, lab in [(0.02, "0.02"), (0.25, "0.25"), (0.48, "0.48")]:
        axR.scatter([], [], s=28 + 620 * np.sqrt(f), color="#8a94a6",
                    edgecolor="white", linewidth=0.7, label=lab)
    axR.legend(title="S$_1$ osc. strength $f$", loc="center right",
               labelspacing=1.3, borderpad=0.9, handletextpad=1.0,
               title_fontsize=7)

    # class legend on left panel
    cls_handles = [
        Line2D([0], [0], marker="o", ls="none", ms=6, mfc=STABLE_INK,
               mec="white", label="stable donors"),
        Line2D([0], [0], marker="o", ls="none", ms=6, mfc=REACTIVE_INK,
               mec="white", label="reactive structures"),
    ]
    axL.legend(handles=cls_handles, loc="lower left", handletextpad=0.3)
    save(fig, "fig6_dipole_switch_synthesis")


# ----------------------------------------------------------------------------
# Fig 7 -- high-quality RDKit structures on class-colored cards, 4 + 3 layout
# ----------------------------------------------------------------------------
def _render_mol(mol_id, px=520):
    mol = Chem.MolFromSmiles(SMILES[mol_id])
    Chem.rdDepictor.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(px, px)
    opt = d.drawOptions()
    opt.bondLineWidth = 2
    opt.padding = 0.08
    opt.clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText())).convert("RGBA")


# ----------------------------------------------------------------------------
# Fig 1 -- graphical abstract with the REAL molecules (RDKit), modern-minimal
# ----------------------------------------------------------------------------
def _mol_chip(ax, mol_id, x0, y0, w, h, ink):
    """Real RDKit structure on a small rounded class-colored chip (axes fraction)."""
    ax.add_patch(FancyBboxPatch((x0, y0), w, h,
                 boxstyle="round,pad=0,rounding_size=0.015",
                 transform=ax.transAxes, facecolor="white", edgecolor=ink,
                 linewidth=1.1, zorder=3))
    img = _render_mol(mol_id, px=340)
    pw, ph = w * 0.09, h * 0.10
    ax.imshow(img, extent=(x0 + pw, x0 + w - pw, y0 + ph, y0 + h - ph - h * 0.16),
              transform=ax.transAxes, zorder=4, aspect="auto")
    ax.text(x0 + w / 2, y0 + h - ph, mol_id, transform=ax.transAxes,
            ha="center", va="top", fontsize=6.2, fontweight="bold", color=ink, zorder=5)


def fig1_graphical_abstract():
    blue, orange, teal = STABLE_INK, REACTIVE_INK, "#009E73"
    bluef, orangef = "#eef5fb", "#fdf3ec"

    fig = plt.figure(figsize=(SINGLE * 2.0, SINGLE * 1.62))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.965, "Seven organic-semiconductor candidates, two computational readouts",
            ha="center", va="center", fontsize=8.8, fontweight="bold", color="#2d3748")

    # ---------- LEFT: predicted protein binding ----------
    ax.add_patch(FancyBboxPatch((0.03, 0.17), 0.45, 0.74,
                 boxstyle="round,pad=0,rounding_size=0.02",
                 facecolor=bluef, edgecolor=blue, linewidth=1.6, zorder=1))
    ax.text(0.255, 0.875, "Predicted protein binding", ha="center", va="center",
            fontsize=8.2, fontweight="bold", color=blue, zorder=6)
    for mid, xx in zip(STABLE, [0.055, 0.16, 0.265, 0.37]):
        _mol_chip(ax, mid, xx, 0.615, 0.098, 0.215, blue)
    ax.add_patch(FancyBboxPatch((0.06, 0.45), 0.39, 0.12,
                 boxstyle="round,pad=0,rounding_size=0.02",
                 facecolor="#c9d6e5", edgecolor="none", alpha=0.75, zorder=2))
    ax.text(0.445, 0.575, "native-inhibitor band", ha="right", va="bottom",
            fontsize=6.3, style="italic", color="#4a5568", zorder=6)
    for i in range(4):
        ax.scatter(0.11 + i * 0.085, 0.51, s=85, marker="o", color=blue,
                   edgecolor="white", linewidth=1.0, transform=ax.transAxes, zorder=5)
    ax.text(0.255, 0.31, "4 stable donors bind within\nthe known-inhibitor range",
            ha="center", va="center", fontsize=6.9, color=INKTXT, zorder=6)

    # ---------- RIGHT: optical environment-sensitivity ----------
    ax.add_patch(FancyBboxPatch((0.52, 0.17), 0.45, 0.74,
                 boxstyle="round,pad=0,rounding_size=0.02",
                 facecolor=orangef, edgecolor=orange, linewidth=1.6, zorder=1))
    ax.text(0.745, 0.875, "Optical environment-sensitivity", ha="center", va="center",
            fontsize=7.8, fontweight="bold", color=orange, zorder=6)
    for mid, xx in zip(REACTIVE, [0.565, 0.695, 0.825]):
        _mol_chip(ax, mid, xx, 0.615, 0.11, 0.215, orange)
    grad = np.linspace(0, 1, 256).reshape(1, -1)
    ax.imshow(grad, extent=(0.565, 0.925, 0.52, 0.57), transform=ax.transAxes,
              aspect="auto", cmap="YlOrBr", zorder=3, alpha=0.9)
    ax.text(0.565, 0.505, "toluene", ha="left", va="top", fontsize=6.1, color="#4a5568", zorder=6)
    ax.text(0.925, 0.505, "water", ha="right", va="top", fontsize=6.1, color="#4a5568", zorder=6)
    ax.annotate("", xy=(0.9, 0.42), xytext=(0.57, 0.42), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=orange, lw=2.2), zorder=5)
    ax.text(0.745, 0.31, "large shift only here —\nbut optically dark states",
            ha="center", va="center", fontsize=6.9, color=INKTXT, zorder=6)

    # ---------- bottom takeaway ----------
    ax.add_patch(FancyBboxPatch((0.14, 0.02), 0.72, 0.10,
                 boxstyle="round,pad=0,rounding_size=0.03",
                 facecolor="#e9f7f1", edgecolor=teal, linewidth=1.6, zorder=2))
    ax.text(0.5, 0.07, "The two readouts select different subsets — reported side by side, not merged",
            ha="center", va="center", fontsize=7.7, fontweight="bold", color=teal, zorder=3)

    # fork + merge arrows
    ax.annotate("", xy=(0.255, 0.915), xytext=(0.45, 0.945), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=blue, lw=1.8), zorder=6)
    ax.annotate("", xy=(0.745, 0.915), xytext=(0.55, 0.945), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=orange, lw=1.8), zorder=6)
    ax.annotate("", xy=(0.45, 0.12), xytext=(0.255, 0.17), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#4a5568", lw=1.6), zorder=6)
    ax.annotate("", xy=(0.55, 0.12), xytext=(0.745, 0.17), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#4a5568", lw=1.6), zorder=6)

    save(fig, "fig1_graphical_abstract_modern")


def fig7_structures():
    # 2 rows x 4 cols; row 0 = 4 stable donors, row 1 = 3 reactive (+1 blank)
    fig = plt.figure(figsize=(DOUBLE, DOUBLE * 0.56))
    gs = fig.add_gridspec(2, 4, hspace=0.28, wspace=0.10,
                          left=0.02, right=0.98, top=0.90, bottom=0.03)
    layout = STABLE + REACTIVE
    positions = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)]

    for mol_id, (r, c) in zip(layout, positions):
        ax = fig.add_subplot(gs[r, c])
        is_stable = mol_id in STABLE
        ink = STABLE_INK if is_stable else REACTIVE_INK
        face = "#f3f8fd" if is_stable else "#fdf5f0"
        # rounded card
        card = FancyBboxPatch((0.02, 0.02), 0.96, 0.96,
                              boxstyle="round,pad=0,rounding_size=0.04",
                              transform=ax.transAxes, facecolor=face,
                              edgecolor=ink, linewidth=1.3, zorder=0)
        ax.add_patch(card)
        img = _render_mol(mol_id)
        ax.imshow(img, extent=(0.06, 0.94, 0.06, 0.80), transform=ax.transAxes,
                  zorder=1, aspect="auto")
        ax.text(0.5, 0.90, mol_id, transform=ax.transAxes, ha="center",
                va="center", fontsize=10, fontweight="bold", color=ink, zorder=2)
        ax.text(0.5, 0.83, f"$\\mu$ = {DIPOLE[mol_id]:.2f} D",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color="#4a5568", zorder=2)
        ax.axis("off")

    # row labels
    fig.text(0.5, 0.955, "Stable donors", ha="center", fontsize=9.5,
             fontweight="bold", color=STABLE_INK)
    fig.text(0.5, 0.475, "Reactive structures", ha="center", fontsize=9.5,
             fontweight="bold", color=REACTIVE_INK)
    save(fig, "fig7_structures_dipole")


if __name__ == "__main__":
    FIGDIR.mkdir(exist_ok=True)
    fig1_graphical_abstract()
    fig2_workflow()
    fig3_docking()
    fig5_solvatochromic()
    fig6_synthesis()
    fig7_structures()
