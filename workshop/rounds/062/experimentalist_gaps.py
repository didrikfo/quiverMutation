"""Round 062 T4: floating VERIFIED_MOVES of width w at the lengths E-165 skipped.
usage: experimentalist_gaps.py L WIDTH [START STEP] [--plan]
(theorist_rulelen.py filters L >= w+5, which excludes w=7 L=10 and w=8 L=11,12; this does not.)"""
import sys, time
from quivermutation import lnaMoves as lm
a = [x for x in sys.argv[1:] if not x.startswith('--')]
L, w = int(a[0]), int(a[1]); start = int(a[2]) if len(a) > 2 else 0; step = int(a[3]) if len(a) > 3 else 1
rules = [r for r in lm.VERIFIED_MOVES if lm.anchorOf(r) is None and r[0] == w][start::step]
print("length", L, "width", w, "rules", len(rules), "slice %d::%d" % (start, step), flush=True)
if '--plan' in sys.argv: sys.exit()
t0 = time.time(); tot = 0; bad = 0
for i, r in enumerate(rules):
    c, f = lm.verifyMove(r, [L]); tot += c; bad += len(f)
    if f: print("FAIL", r, f[:2], flush=True)
    print(i, tot, bad, "%.0fs" % (time.time()-t0), flush=True)
print("DONE length %d width %d slice %d::%d rules %d confirmed %d failures %d (%.0fs)" % (L, w, start, step, len(rules), tot, bad, time.time()-t0))
