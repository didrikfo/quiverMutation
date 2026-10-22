"""H-010, two bystanders, k = 3 cell.  Pair (8:3)(9:3) (arrows 8..11) plus TWO bystanders
(t1:m1),(t2:m2), each sharing <= 1 arrow with the pair (so neither alone is a "run of three"),
m in {2,3}.  Enumerates all admissible mutation sequences of <= STEPS mutations at vertices within
MARGIN of the pattern (lnaMoves.localMutationSequences) and reports base max overlap, the least
max overlap reached, and how many reached LNAs are lower.
Usage: theorist_twoby.py STEPS MARGIN TLO THI [--plan] [--jobs J]   (t1<t2 over [TLO,THI])"""
import sys, time, itertools
from multiprocessing import Pool
from quivermutation import lnaMoves as lm, overlap as ov

PAIR = [(8, 3), (9, 3)]
PAIRARR = set(range(8, 12))

def cells(tlo, thi):
    out = []
    for m1, m2 in itertools.product((2, 3), repeat=2):
        for t1 in range(tlo, thi + 1):
            for t2 in range(t1 + 1, thi + 1):
                if t1 in (8, 9) or t2 in (8, 9): continue
                a1, a2 = set(range(t1, t1 + m1)), set(range(t2, t2 + m2))
                if len(a1 & PAIRARR) > 1 or len(a2 & PAIRARR) > 1: continue
                out.append(((t1, m1), (t2, m2)))
    return out

def run(args):
    b1, b2, steps, margin = args
    pat = PAIR + [b1, b2]
    lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
    off = max(0, 7 - lo)
    pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
    L = hi + margin + 7
    rl = lm.embedPattern(L, pat, 0)
    if rl is None: return (b1, b2, None)
    t0 = time.time()
    reached = lm.localMutationSequences(L, rl, lo, hi, steps, margin)
    base = ov.maxOverlap(rl)
    ovs = [ov.maxOverlap([int(c) for c in n]) for n in reached]
    return (b1, b2, base, len(reached), min(ovs) if ovs else None, sum(o < base for o in ovs), round(time.time() - t0, 1))

if __name__ == "__main__":
    steps, margin, tlo, thi = map(int, sys.argv[1:5])
    jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 4
    cs = cells(tlo, thi)
    print("cells", len(cs))
    if "--plan" in sys.argv: sys.exit()
    print("b1 b2 base reached min lowered secs", flush=True)
    with Pool(jobs) as p:
        for r in p.imap_unordered(run, [(a, b, steps, margin) for a, b in cs]):
            print(*r, flush=True)
