"""Skeptic r118, item 14 (T4/H-020 H3/H4): is the first-step relation delta of a notLNA class-B rule a function of the RULE alone
(same value on every matching row, every length, every window position)?  Pre-registered: YES = for every rule one delta pair over
all matches at all lengths run; NO = some rule with two deltas.  Then the held-out transfer: delta from lengths <= Lmax-1 predicts Lmax.
Also: which rule features (len(added), len(d[1]), width) predict the pair, tried as functions of rule data only.
usage: skeptic_ruledelta.py W1,W2 LMAX"""
import sys, collections, json
sys.path.insert(0, 'workshop/rounds/095')
from theorist_rowmodel import *
W = [int(a) for a in sys.argv[1].split(',')]; LMAX = int(sys.argv[2])
for w in W:
    flo = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == w]
    perrule = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    for i, d in enumerate(flo):
        for L in range(w + 1, LMAX + 1):
            for r in rowsOf(L):
                for ws in lm.windowStartsFor(L, d):
                    if not lm.matchesAt(L, r, d, ws): continue
                    s0 = lm.applyAt(L, r, d, ws)[1][0]
                    st, _ = step(L, r, s0)
                    if st != 'notLNA': continue
                    alg = nakayama.LinearNakayamaAlgebra(L, r)
                    m = lm._quiet(mutation.quiverMutationAtVertices, lm._copy(alg), [s0])
                    y = (len(list(m.arrows())) - len(list(alg.arrows())), len(m.rels) - len(alg.rels))
                    perrule[i][L][y] += 1
    nconst = nmulti = 0; trans_ok = trans_n = 0; rows = []
    for i, byL in perrule.items():
        allv = set(y for c in byL.values() for y in c)
        if len(allv) == 1: nconst += 1
        else: nmulti += 1; rows.append((i, {L: dict(c) for L, c in byL.items()}))
        Ls = sorted(byL)
        for L in Ls[1:]:
            hist = set(y for K in Ls if K < L for y in byL[K])
            if len(hist) == 1:
                trans_n += 1; trans_ok += set(byL[L]) <= hist
    print('width', w, 'L<=', LMAX, 'notLNA rules', len(perrule), 'one delta over all matches', nconst, 'several deltas', nmulti,
          'length-transfer (single hist delta -> next L)', trans_ok, '/', trans_n)
    dist = collections.Counter(next(iter(set(y for c in byL.values() for y in c))) for byL in perrule.values() if len(set(y for c in byL.values() for y in c)) == 1)
    nm = collections.Counter(sum(sum(c.values()) for c in byL.values()) for byL in perrule.values())
    print('  delta distribution over rules', dict(dist), '; matches per rule', dict(sorted(nm.items())))
    # control: coarse keys that must fail to be functions
    for name, kf in (('len(added)', lambda d: len(d[2])), ('len(d[1])', lambda d: len(d[1])), ('width only', lambda d: 0)):
        grp = collections.defaultdict(set)
        for i, byL in perrule.items():
            for c in byL.values(): grp[kf(flo[i])] |= set(c)
        print('  key', name, ': groups', len(grp), 'functional groups', sum(len(v) == 1 for v in grp.values()))
    json.dump([[list(flo[i][:3]) , sorted(set(y for c in byL.values() for y in c))[0]] for i, byL in perrule.items()], open('workshop/rounds/118/skeptic_ruledelta_w%d.json' % w, 'w'))
    for i, dd in rows[:6]: print('  multi', flo[i][1:3], dd)
