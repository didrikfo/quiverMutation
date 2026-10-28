"""T7 / H-010 at k = 5 (E-241 extended), pair (8:3)(9:3) plus two bystanders. Adapted from 086/experimentalist_k4.py.
Usage: experimentalist_k5.py STEPS MARGIN CELLS [--jobs J] [--plan]
CELLS: semicolon list like "6,2:11,2;7,3:12,2" (b1:b2, each t,m).
--plan prints each cell with the k = 4 seconds recorded in E-241 and a x4.5 per step estimate; runs nothing."""
import sys
sys.path.insert(0, "workshop/rounds/081")
from multiprocessing import Pool
import theorist_twoby as t
K4 = {"11,2:12,2": 175, "6,2:7,2": 311, "7,2:11,2": 335, "6,2:11,2": 505, "7,2:12,2": 507,
      "7,3:12,2": 355, "10,3:5,2": 533, "6,4:12,2": 511}
if __name__ == "__main__":
    steps, margin = int(sys.argv[1]), int(sys.argv[2])
    raw = sys.argv[3].split(";")
    cs = [tuple(tuple(map(int, b.split(","))) for b in c.split(":")) for c in raw]
    jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 4
    print("steps", steps, "margin", margin, "cells", len(cs), flush=True)
    if "--plan" in sys.argv:
        for c in raw:
            k4 = K4.get(c)
            print(c, "k4 secs", k4, "est k%d secs" % steps, round(k4 * 4.5 ** (steps - 4)) if k4 else "?")
        sys.exit()
    print("b1 b2 base reached min lowered secs", flush=True)
    with Pool(jobs) as p:
        for r in p.imap_unordered(t.run, [(a, b, steps, margin) for a, b in cs]): print(*r, flush=True)
