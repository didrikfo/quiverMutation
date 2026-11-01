"""Skeptic r106 response: out-of-width test of the post-hoc rule found on width 5: notLNA first step has (arrows, rels) delta = (1, nadd-1),
nadd = number of relation entries in the rule's ADDED window pair list. Pre-registered here before running on widths 6 and 7."""
import sys, collections
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
W = [int(a) for a in sys.argv[1].split(',')]; LMAX = int(sys.argv[2])
for w in W:
    flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == w]
    res = collections.Counter(); rules = collections.defaultdict(set); bad = collections.Counter()
    for i, d in enumerate(flo):
        for L in range(w + 1, LMAX + 1):
            for r in rowsOf(L):
                for ws in lm.windowStartsFor(L, d):
                    if not lm.matchesAt(L, r, d, ws): continue
                    s0 = lm.applyAt(L, r, d, ws)[1][0]
                    st, _ = step(L, r, s0)
                    if st != 'notLNA': res['other:' + st] += 1; continue
                    alg = nakayama.LinearNakayamaAlgebra(L, r)
                    m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s0])
                    y = (len(list(m.arrows())) - len(list(alg.arrows())), len(m.rels) - len(alg.rels))
                    pred = (1, len(d[2]) - 1)
                    res['hit' if y == pred else 'miss'] += 1; res['arrow+1' if y[0] == 1 else 'arrowOTHER'] += 1; rules[i].add(y)
                    if y != pred: bad[(d[1], d[2], y, pred)] += 1
    print('width', w, 'lengths', w + 1, '..', LMAX, 'rules', len(flo), 'notLNA rules', len(rules), dict(res))
    for k, v in list(bad.items())[:8]: print('  miss', k, v)
