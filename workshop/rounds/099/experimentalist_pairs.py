"""T7 / H-010, other pair shapes (not (8:3)(9:3)), two bystanders, k steps. Extends E-230/E-241.
Usage: experimentalist_pairs.py STEPS MARGIN PAIR CELLS [--jobs J] [--plan]
PAIR like "8,4:10,4" (start,len each); CELLS like "6,2:11,2;7,2:12,2" (bystanders, each start,len)."""
import sys, time
from multiprocessing import Pool
from quivermutation import lnaMoves as lm, overlap as ov

def run(args):
    pair, b1, b2, steps, margin = args
    pat = list(pair) + [b1, b2]
    lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
    off = max(0, 7 - lo); pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
    L = hi + margin + 7
    rl = lm.embedPattern(L, pat, 0)
    if rl is None: return (pair, b1, b2, None)
    t0 = time.time()
    reached = lm.localMutationSequences(L, rl, lo, hi, steps, margin)
    base = ov.maxOverlap(rl)
    ovs = [ov.maxOverlap([int(c) for c in n]) for n in reached]
    return (pair, b1, b2, "base", base, "reached", len(reached), "min", min(ovs) if ovs else None,
            "lowered", sum(o < base for o in ovs), "secs", round(time.time() - t0, 1))

def parse(c): return tuple(tuple(map(int, b.split(","))) for b in c.split(":"))
if __name__ == "__main__":
    steps, margin = int(sys.argv[1]), int(sys.argv[2])
    pair = parse(sys.argv[3]); cs = [parse(c) for c in sys.argv[4].split(";")]
    jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 4
    pa = set(a for s, k in pair for a in range(s, s + k))
    print("pair", pair, "steps", steps, "margin", margin, "cells", len(cs), flush=True)
    for a, b in cs:
        sh = [len(set(range(s, s + k)) & pa) for s, k in (a, b)]
        print(" cell", a, b, "shared arrows with pair", sh, "bystanders share", len(set(range(a[0], a[0]+a[1])) & set(range(b[0], b[0]+b[1]))))
    if "--plan" in sys.argv: sys.exit()
    with Pool(jobs) as p:
        for r in p.imap_unordered(run, [(pair, a, b, steps, margin) for a, b in cs]): print(*r, flush=True)
