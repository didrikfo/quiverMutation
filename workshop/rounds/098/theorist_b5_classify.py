"""T4 round 098: does the first step alone (E-187 right / E-205 left) decide class A vs B for the 24 floating width-5 rules?
Prediction P: rule is B iff at some applicable match the row-model first step is 'notLNA'. Compare with the full replay class (lengths 6..8).
Also: for B rules, the index of the failing step and, for the 'notLNA' middle algebra, its arrow count (engine)."""
import sys, collections
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == 5]
agree = collections.Counter()
for d in flo:
    full = collections.Counter(); first = collections.Counter(); failstep = collections.Counter(); arrows = collections.Counter()
    for L in range(6, 9):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if not lm.matchesAt(L, r, d, ws): continue
                o = replay(L, r, d, ws)
                full['A' if o == 'agree' else ('noapply' if o == 'noapply' else 'B')] += 1
                if o == 'noapply': continue
                if o != 'agree': failstep[o[0]] += 1
                st, _ = step(L, r, lm.applyAt(L, r, d, ws)[1][0]); first[st] += 1
    cls = 'A' if full['A'] and not full['B'] else ('B' if full['B'] and not full['A'] else 'mixed')
    pred = 'B' if first['notLNA'] else 'A'
    agree[(cls, pred)] += 1
    print(d[1:], 'class', cls, 'first-step', dict(first), 'failing step idx', dict(failstep), 'noapply', full['noapply'])
print('(replay class, first-step prediction):', dict(agree))
