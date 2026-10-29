"""T4 round 098: extend theorist_b5_classify to all 414 floating rules, lengths w+1..w+2 only.
Per rule: replay class (A all agree / B none agree / C mixed / D no match) vs first-step-only prediction (B iff some match has first step notLNA).
Also the step index at which B rules fail (0 = first step)."""
import sys, collections
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None]
tab = collections.Counter(); late = collections.Counter(); bad = []
for d in flo:
    full = collections.Counter(); first = collections.Counter(); fs = collections.Counter()
    for L in range(d[0] + 1, d[0] + 3):
        for r in rowsOf(L):
            for ws in lm.windowStartsFor(L, d):
                if not lm.matchesAt(L, r, d, ws): continue
                o = replay(L, r, d, ws)
                if o == 'noapply': continue
                full['A' if o == 'agree' else 'B'] += 1
                if o != 'agree': fs[o[0]] += 1
                first[step(L, r, lm.applyAt(L, r, d, ws)[1][0])[0]] += 1
    cls = 'D' if not full else ('A' if not full['B'] else ('B' if not full['A'] else 'C'))
    pred = 'B' if first['notLNA'] else ('A' if first else 'D')
    tab[(d[0], cls, pred)] += 1
    if cls != pred: bad.append((d, dict(full), dict(first), dict(fs)))
    if cls == 'B': late[(d[0], 'fail at step 0' if set(fs) == {0} else 'later step involved ' + str(sorted(fs)))] += 1
for k in sorted(tab): print('width', k[0], 'replay class', k[1], 'first-step pred', k[2], tab[k])
for k in sorted(late): print(k, late[k])
print('mismatching rules', len(bad))
for b in bad[:12]: print(b)
