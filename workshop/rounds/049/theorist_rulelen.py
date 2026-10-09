"""T4: is each floating rule of lnaMoves.VERIFIED_MOVES length-independent?
Re-verify rules at lengths beyond those the test suite uses (width+1..width+4).
usage: theorist_rulelen.py MAXLEN [MAXWIDTH]  -- checks, for each floating rule of
window width w <= MAXWIDTH, the lengths w+5..MAXLEN (all admissible LNAs).
"""
import sys, time
from collections import Counter
from quivermutation import lnaMoves as lm
maxlen = int(sys.argv[1]); maxw = int(sys.argv[2]) if len(sys.argv) > 2 else 4
rules = [r for r in lm.VERIFIED_MOVES if lm.anchorOf(r) is None]
print("floating rules", len(rules), "widths", sorted(Counter(r[0] for r in rules).items()))
t0 = time.time(); tot = 0; bad = 0
for r in rules:
    w = r[0]
    if w > maxw: continue
    lens = range(w + 5, maxlen + 1)
    if not lens: continue
    c, f = lm.verifyMove(r, lens)
    tot += c; bad += len(f)
    if f: print("FAIL", r, f[:2])
print("checked rules width<=%d at lengths w+5..%d: confirmed %d failures %d (%.0fs)" % (maxw, maxlen, tot, bad, time.time()-t0))
