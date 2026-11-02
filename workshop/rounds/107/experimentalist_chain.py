"""T7 / H-010 round 107: lowering controls for gap-2 pair shapes, and k>=5 cells.
Usage: experimentalist_chain.py STEPS MARGIN "s,l:s,l:s,l;..." [--jobs J]
Each cell is a list of relations (start,len) embedded as given (shifted so min start >= 7), all
mutated within MARGIN of the whole pattern. Prints base overlap, #LNAs reached, min overlap, #lowered,
and the lowest-overlap LNA (relation lengths). Wraps lnaMoves.localMutationSequences (as 099 pairs)."""
import sys, time
from multiprocessing import Pool
from quivermutation import lnaMoves as lm, overlap as ov

def run(a):
    pat, steps, margin = a
    pat = list(pat)
    lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
    off = max(0, 7 - lo); pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
    L = hi + margin + 7
    rl = lm.embedPattern(L, pat, 0)
    if rl is None: return (a[0], None)
    t0 = time.time()
    reached = lm.localMutationSequences(L, rl, lo, hi, steps, margin)
    base = ov.maxOverlap(rl)
    ovs = {n: ov.maxOverlap([int(c) for c in n]) for n in reached}
    best = min(ovs, key=ovs.get) if ovs else None
    return (a[0], "base", base, "reached", len(reached), "min", ovs[best] if ovs else None,
            "lowered", sum(o < base for o in ovs.values()), "secs", round(time.time()-t0,1))

if __name__ == "__main__":
    steps, margin = int(sys.argv[1]), int(sys.argv[2])
    cells = [tuple(tuple(map(int, r.split(","))) for r in c.split(":")) for c in sys.argv[3].split(";")]
    jobs = int(sys.argv[sys.argv.index("--jobs")+1]) if "--jobs" in sys.argv else 4
    print("steps", steps, "margin", margin, "cells", len(cells), flush=True)
    with Pool(jobs) as p:
        for r in p.imap_unordered(run, [(c, steps, margin) for c in cells]): print(*r, flush=True)
