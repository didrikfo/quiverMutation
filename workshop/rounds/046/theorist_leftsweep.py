"""Pair (8:3)(9:3) plus ONE bystander (t:m), t=2..16 (left and right, overlaps allowed),
m=2..4, STEPS mutations, MARGIN.  Pattern shifted so >=6 empty vertices lie left.
Usage: theorist_leftsweep.py STEPS MARGIN [tmin tmax]"""
import sys, time
from quivermutation import lnaMoves as lm, overlap as ov
steps, margin = int(sys.argv[1]), int(sys.argv[2])
tmin, tmax = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (2, 16)
print("t m shared_with_pair  base reached min lowered")
for m in (2, 3, 4):
    for t in range(tmin, tmax + 1):
        pat = [(8, 3), (9, 3), (t, m)]
        if t in (8, 9): continue
        lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
        off = max(0, 7 - lo)
        pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
        L = hi + margin + 7
        rl = lm.embedPattern(L, pat, 0)
        if rl is None:
            print(t, m, "inadmissible"); continue
        arrows = set(range(t, t + m)); sh = len(arrows & set(range(8, 12)))
        t0 = time.time()
        reached = lm.localMutationSequences(L, rl, lo, hi, steps, margin)
        base = ov.maxOverlap(rl)
        ovs = [ov.maxOverlap([int(c) for c in n]) for n in reached]
        print(t, m, sh, f"base={base}", len(reached), min(ovs) if ovs else '-', sum(o < base for o in ovs), f"({time.time()-t0:.0f}s)", flush=True)
