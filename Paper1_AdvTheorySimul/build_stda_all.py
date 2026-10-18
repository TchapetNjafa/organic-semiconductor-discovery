"""Collect the lowest sTDA-xTB singlet and the lowest bright singlet of every PubChemQC shard molecule.

Reads NTO_and_related_study/data/xtb_results/<mol_id>/tda.dat (DATXY block: state, E [eV], oscillator strength f, ...).
S1 = lowest state; bright = lowest state with oscillator strength f >= F_BRIGHT. Molecules without a DATXY block are skipped.
usage: python -I build_stda_all.py <project_root> [F_BRIGHT]   -> outputs/extract/stda_all.csv
"""
import sys
from pathlib import Path

import pandas as pd

F_BRIGHT = 0.01


def read_tda(path, f_bright):
    lines = path.read_text().splitlines()
    try:
        i = next(k for k, l in enumerate(lines) if l.strip() == "DATXY")
    except StopIteration:
        return None
    states = []
    for l in lines[i + 1:]:
        p = l.split()
        if len(p) < 3:
            break
        states.append((float(p[1]), float(p[2])))
    if not states:
        return None
    bright = next(((e, f) for e, f in states if f >= f_bright), (float("nan"), float("nan")))
    return states[0][0], states[0][1], bright[0], bright[1]


if __name__ == "__main__":
    root = Path(sys.argv[1])
    fb = float(sys.argv[2]) if len(sys.argv) > 2 else F_BRIGHT
    rows = []
    for d in sorted((root / "NTO_and_related_study/data/xtb_results").iterdir(), key=lambda p: p.name):
        t = d / "tda.dat"
        if d.name.isdigit() and t.exists() and (r := read_tda(t, fb)) is not None:
            rows.append((int(d.name), *r))
    out = pd.DataFrame(rows, columns=["mol_id", "E_S1_eV", "f_S1", "E_bright_eV", "f_bright"])
    dest = Path(sys.argv[3]) if len(sys.argv) > 3 else root / "outputs/extract/stda_all.csv"
    out.to_csv(dest, index=False)
    print(len(out), "molecules ->", dest)
