"""T6 (round 022, theorist): test the 'peeling' formula for the minimum cord depth of a doubly blocked LNA:
depth = 1 + min over relations of >= 3 arrows of min(a, b), a = number of consecutive 2-arrow relations starting at x0-1, x0-2, ...,
b = number of consecutive 2-arrow relations starting at x_m-1, x_m, ... (x0 = start, x_m = x0 + m).
usage: theorist_peel.py FILE...   (lines 'k digits pred p found f depth d' for LNAs with depth >= 2; repository root)."""
import sys
def formula(d):
    r = [0] + list(map(int, d)) + [0] * 5
    best = None
    for s in range(1, len(d) + 1):
        m = r[s]
        if m < 3: continue
        a = 0
        while s - 1 - a >= 1 and r[s - 1 - a] == 2: a += 1
        b = 0
        while r[s + m - 1 + b] == 2: b += 1
        v = 1 + min(a, b); best = v if best is None else min(best, v)
    return best
tot = bad = 0
for f in sys.argv[1:]:
    for line in open(f):
        t = line.split()
        if len(t) < 8 or t[2] != 'pred' or t[5] != '1': continue
        tot += 1
        if formula(t[1]) != int(t[7]):
            bad += 1; print('mismatch', len(t[1]) + 2, t[1], 'data', t[7], 'formula', formula(t[1]))
print('LNAs', tot, 'mismatches', bad)
