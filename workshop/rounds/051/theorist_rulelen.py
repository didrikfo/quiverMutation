"""T4/H1: floating rules of width w in WIDTHS re-verified at exactly one length L (>= w+5).
usage: theorist_rulelen.py L WIDTHS(comma) [START STEP]   -- rules sliced [START::STEP] of the chosen set
  --plan prints counts only.  Prints per-rule failures and a total."""
import sys, time
from collections import Counter
from quivermutation import lnaMoves as lm
L = int(sys.argv[1]); widths = [int(x) for x in sys.argv[2].split(',')]
start = int(sys.argv[3]) if len(sys.argv) > 3 else 0
step = int(sys.argv[4]) if len(sys.argv) > 4 else 1
rules = [r for r in lm.VERIFIED_MOVES if lm.anchorOf(r) is None and r[0] in widths and L >= r[0] + 5]
rules = rules[start::step]
print("length", L, "widths", widths, "rules", len(rules), Counter(r[0] for r in rules), flush=True)
t0 = time.time(); tot = 0; bad = 0
for i, r in enumerate(rules):
    c, f = lm.verifyMove(r, [L])
    tot += c; bad += len(f)
    if f: print("FAIL", r, f[:2], flush=True)
    print(i, tot, bad, "%.0fs" % (time.time()-t0), flush=True)
print("DONE length %d widths %s slice %d::%d rules %d confirmed %d failures %d (%.0fs)" % (L, widths, start, step, len(rules), tot, bad, time.time()-t0))
