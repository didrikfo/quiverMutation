"""E-152 cap-cut k = 3 rows.  Pair (8:3)(9:3) plus ONE bystander (t:m), STEPS mutations, MARGIN 2.
Right side: bystander first arrow t = 12 + g (g = arrows between pair end (arrow 11) and t; same grid as
theorist_bystander.py).  Left side: t = 8 - m - g ... (pattern shifted right so >= 6 empty vertices lie left, as in
theorist_leftsweep.py).  Usage: capcut.py SIDE(right|left) STEPS M G_LO G_HI [--plan]
Output per row: base overlap, number of LNAs reached, min overlap reached, number lowered."""
import sys, time
from quivermutation import lnaMoves as lm, overlap as ov
side, steps, m, glo, ghi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
plan = '--plan' in sys.argv
margin = 2
for g in range(glo, ghi + 1):
    if side == 'right':
        t = 12 + g
    else:
        t = 8 - g - m          # bystander arrows t..t+m-1, g empty arrows between it and arrow 8
    pat = [(8, 3), (9, 3), (t, m)]
    lo = min(s for s, _ in pat); hi = max(s + k - 1 for s, k in pat)
    off = max(0, 7 - lo)
    pat = [(s + off, k) for s, k in pat]; lo += off; hi += off
    L = hi + margin + (7 if side == 'left' else 6)
    rl = lm.embedPattern(L, pat, 0)
    if rl is None:
        print(side, m, g, "inadmissible", flush=True); continue
    arrows = set(range(t + off, t + off + m)); sh = len(arrows & set(range(8 + off, 12 + off)))
    if plan:
        print(side, m, g, "pattern", pat, "L", L, "shared arrows", sh); continue
    t0 = time.time()
    reached = lm.localMutationSequences(L, rl, lo, hi, steps, margin)
    base = ov.maxOverlap(rl)
    ovs = [ov.maxOverlap([int(c) for c in nme]) for nme in reached]
    print(f"{side} k={steps} m={m} g={g} shared={sh} base={base} reached={len(reached)} min={min(ovs) if ovs else '-'} lowered={sum(o < base for o in ovs)} ({time.time()-t0:.0f}s)", flush=True)
