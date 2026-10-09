"""H-010 locality scan: an isolated pair (1:3)(2:3) planted in the interior of a
line, plus ONE bystander relation of m arrows to its right at gap g (g = number of
arrows between the pair's last arrow and the bystander's first; g < 0 = overlap).
Enumerates every admissible mutation sequence of <= steps mutations at vertices
within `margin` of the pair+bystander, and reports the least maximum overlap of an
LNA reached.  Usage: theorist_bystander.py STEPS MARGIN."""
import sys, time
from quivermutation import lnaMoves as lm, overlap as ov

steps, margin = int(sys.argv[1]), int(sys.argv[2])
P = 8                        # pair starts at vertex P: relations (P:3),(P+1:3); arrows P..P+3
last = P + 3
print("steps", steps, "margin", margin)
print("m  g  start  bystander-arrows  nReached  minMaxOverlap  lowered(count)")
for m in (2, 3, 4):
    for g in range(-2, 6):
        t = last + 1 + g          # bystander first arrow
        if t <= P + 1: continue
        end = t + m - 1
        L = end + margin + 6
        pattern = [(P, 3), (P + 1, 3), (t, m)]
        rl = lm.embedPattern(L, pattern, 0)
        if rl is None:
            continue
        t0 = time.time()
        reached = lm.localMutationSequences(L, rl, P, end, steps, margin)
        base = ov.maxOverlap(rl)
        ovs = [ov.maxOverlap([int(c) for c in name]) for name in reached]
        low = sum(1 for o in ovs if o < base)
        print(f"{m}  {g:2d}  {t:3d}  base={base}  {len(reached):6d}  {min(ovs) if ovs else '-'}  {low}  ({time.time()-t0:.0f}s)", flush=True)
