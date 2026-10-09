"""T6 (round 022, theorist): fit the depth-1 data of theorist_cordcrit.py (L = 1 files) with the 'blocked ends' rule:
a relation of m >= 3 arrows from x0 to x0+m is LEFT-blocked if a relation starts at x0-1 and RIGHT-blocked if a relation ends at x0+m+1;
variant 'interior': x0 (resp. x0+m) is an interior vertex of some relation. rule: a depth-1 cord exists iff some relation of >= 3 arrows is not blocked on both sides.
usage: theorist_blocked.py FILE...   (files of lines 'k digits pred p found f depth d'; repository root)."""
import sys
def rels(d):
    return [(i + 1, a) for i, a in enumerate(map(int, d)) if a]
def rule(d, variant):
    R = rels(d); starts = {s: a for s, a in R}; ends = {s + a for s, a in R}
    for s, a in R:
        if a < 3: continue
        if variant == 'sym':
            L = (s + 1) in ends; Rb = (s + a - 1) in starts
        elif variant == 'interior':
            L = any(s2 < s < s2 + a2 for s2, a2 in R); Rb = any(s2 < s + a < s2 + a2 for s2, a2 in R)
        else:
            L = (s - 1) in starts
            Rb = (s + a + 1) in ends if variant == 'end' else (s + a - 1) in starts
        if not (L and Rb): return True
    return False
for f in sys.argv[1:]:
    for variant in ('end', 'start', 'interior', 'sym'):
        bad = []; tot = 0
        for line in open(f):
            t = line.split()
            if len(t) < 8 or t[2] != 'pred': continue
            tot += 1; d1 = (t[7] == '1')
            if rule(t[1], variant) != d1: bad.append(t[1])
        print(f, variant, 'LNAs', tot, 'mismatches', len(bad), bad[:12])
