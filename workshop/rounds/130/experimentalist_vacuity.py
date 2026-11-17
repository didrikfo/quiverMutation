"""T7/H-010 round 130 (answer to referee): vacuity-aware counts for experimentalist_screen.
Imports experimentalist_screen.screen / cells unchanged. A cell is vacuous when the closure reaches
no rewrite at all (reached == 0). Prints, per pair shape and per family, clean / vacuous /
non-vacuous / lowering at the given depth; with --two the two-bystander run on the overlap-2 pairs;
with --iso the number of rewrites reached on each isolated pair at depths 1..3.
Usage: experimentalist_vacuity.py DEPTH [--two] [--iso]"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import defaultdict
from experimentalist_screen import screen, cells
from quivermutation import overlap as ov

SHAPES = [((8, a), (10, b)) for a in range(3, 8) for b in range(2, 8) if 8 + a < 10 + b]

def tab(pairs, bys, mmax, depth, two):
    per = defaultdict(lambda: [0, 0, 0])   # clean, vacuous, lowering
    for pair, by, pat in cells(pairs, bys, mmax, two):
        r = screen(pat, depth)
        if r is None:
            continue
        pov = pair[0][0] + pair[0][1] - pair[1][0]
        if r["base"] != pov or any(abs(s - p) == 1 for s, _ in by for p, _ in pair):
            continue
        c = per[pair]; c[0] += 1
        if r["reached"] == 0:
            c[1] += 1
        if r["min"] < r["base"]:
            c[2] += 1
    return per

def fam(pair):
    (s1, a), (s2, b) = pair
    ovl = a - 2
    return ("b=a-1" if b == a - 1 else "b>=a") + f" ov{ovl}"

if __name__ == "__main__":
    depth = int(sys.argv[1])
    if "--iso" in sys.argv:
        for pair in SHAPES:
            print("isolated", pair, "reached by depth 1..3:",
                  [screen(list(pair), d)["reached"] for d in (1, 2, 3)])
        sys.exit()
    two = "--two" in sys.argv
    pairs = [((8, 4), (10, 3)), ((8, 4), (10, 4))] if two else SHAPES
    per = tab(pairs, range(2, 18 if two else 21), 7 if two else 8, depth, two)
    print(f"## depth {depth}, {'two bystanders, overlap-2 pairs' if two else 'one bystander, 20 shapes'}")
    print("pair | pairOv | clean | vacuous | non-vacuous | lower")
    agg = defaultdict(lambda: [0, 0, 0]); byov = defaultdict(lambda: [0, 0, 0])
    for pair, (c, v, l) in sorted(per.items()):
        print(pair, "|", pair[0][1] - 2, "|", c, "|", v, "|", c - v, "|", l)
        for d, k in ((agg, fam(pair)), (byov, pair[0][1] - 2)):
            for i, x in enumerate((c, v, l)):
                d[k][i] += x
    print("family | clean | vacuous | non-vacuous | lower")
    for k, (c, v, l) in sorted(agg.items()):
        print(k, "|", c, "|", v, "|", c - v, "|", l)
    print("pair overlap | clean | vacuous | non-vacuous | lower")
    for k, (c, v, l) in sorted(byov.items()):
        print(k, "|", c, "|", v, "|", c - v, "|", l)
    t = [sum(x[i] for x in per.values()) for i in range(3)]
    print("total | clean", t[0], "| vacuous", t[1], "| non-vacuous", t[0] - t[1], "| lower", t[2])
