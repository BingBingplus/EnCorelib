#!/usr/bin/env python3
"""Compare a table2.sh run against the printed Table 2, cell by cell.

    make table2                 # or: bash repro/table2.sh eg2 eg4 eg6
    python repro/compare.py

Reads repro/out/table2.txt (what this code produced) and
repro/table2_paper.txt (the printed table, transcribed) and reports, for
every cell present in both, whether B1, rho and #(C) agree.

Sampling a few examples is enough to see whether the build reproduces the
paper; `bash repro/table2.sh eg2 eg4 eg6` takes about a minute.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OURS = os.path.join(HERE, "out", "table2.txt")
PAPER = os.path.join(HERE, "table2_paper.txt")

B1 = r"\(.*?\) \+- \(.*?\)"
PAPER_ROW = re.compile(r"^(Eg\d)\s+(\S+)\s+Ours\(([\d.]+)\)\s+([\d.]+)\s+(" + B1 + r")\s+(\S+)\s+(\S+)\s+(\S+)")
OURS_ROW = re.compile(r"^(Eg\d)\s+(\S+)\s+([\d.]+)\s+([\d.]+)\s+(" + B1 + r")\s+(\S+)\s+(\S+)\s+(\S+)")


def read(path, pattern, eps_at, h_at):
    out = {}
    with open(path) as fh:
        for line in fh:
            m = pattern.match(line)
            if m:
                g = m.groups()
                key = (g[0], g[eps_at], "%g" % float(g[h_at]))
                out[key] = (g[4], g[5], g[7])       # B1, rho, #(C)
    return out


def same_rho(a, b):
    try:
        return abs(float(a) - float(b)) < 5e-3
    except ValueError:
        return a == b


def main():
    for p in (OURS, PAPER):
        if not os.path.exists(p):
            sys.exit("missing %s -- run `make table2` first" % p)

    ours = read(OURS, OURS_ROW, 2, 3)
    paper = read(PAPER, PAPER_ROW, 2, 3)

    shared = sorted(k for k in ours if k in paper)
    if not shared:
        sys.exit("no cells in common -- did table2.sh run with METHOD=2?")

    print("%-5s %-5s %-5s  %s" % ("case", "eps", "H", "printed Table 2  ->  this build"))
    print("-" * 96)

    agree = 0
    for key in shared:
        eg, eps, h = key
        pb, pr, pc = paper[key]
        ob, orr, oc = ours[key]
        notes = []
        if pb != ob:
            notes.append("B1 %s -> %s" % (pb, ob))
        if not same_rho(pr, orr):
            notes.append("rho %s -> %s" % (pr, orr))
        if pc != oc:
            notes.append("#(C) %s -> %s" % (pc, oc))
        if notes:
            print("%-5s %-5s %-5s  DIFF  %s" % (eg, eps, h, ";  ".join(notes)))
        else:
            print("%-5s %-5s %-5s  ok    %s   rho=%s  #(C)=%s" % (eg, eps, h, ob, orr, oc))
            agree += 1

    print()
    print("%d of %d compared cells reproduce the printed table exactly."
          % (agree, len(shared)))
    if agree != len(shared):
        print("Cells marked DIFF: see the note at the bottom of repro/table2_paper.txt,")
        print("which lists the entries of the printed table that are inconsistent with")
        print("their own rho column.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
