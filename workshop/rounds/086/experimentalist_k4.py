"""T7 / H-010 at k = 4 (and margin option), pair (8:3)(9:3) plus two bystanders. Extends E-230.
Usage: experimentalist_k4.py STEPS MARGIN CELLS [--jobs J] [--plan]
CELLS: semicolon list like "6,2:11,2;7,3:12,2" (b1:b2, each t,m)."""
import sys, os
sys.path.insert(0, "workshop/rounds/081")
from multiprocessing import Pool
import theorist_twoby as t
if __name__ == "__main__":
    steps, margin = int(sys.argv[1]), int(sys.argv[2])
    cs = [tuple(tuple(map(int, b.split(","))) for b in c.split(":")) for c in sys.argv[3].split(";")]
    jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 4
    print("steps", steps, "margin", margin, "cells", len(cs), flush=True)
    print("b1 b2 base reached min lowered secs", flush=True)
    with Pool(jobs) as p:
        for r in p.imap_unordered(t.run, [(a, b, steps, margin) for a, b in cs]): print(*r, flush=True)
