"""Round 022 (experimentalist): control for E-099's long-sided square.
Guarded BFS as rounds/021/scholar_step7_entries.py (class C, from LNAs and duals), but at EVERY gate-admitted step (parent, v)
record tiltingPlus and the two shape tests of E-099 at v (strict A5 = hasA5; long square = hasLongSquare, both copied from
rounds/021/scholar_step7_entries.py), de-duplicated on (canonical key of parent, v). Also: v has exactly one out-arrow.
Usage: experimentalist_shapectl.py n [--class I] [--budget-sec S] [--maxexp N]"""
import argparse, sys, time
from collections import Counter
sys.path.insert(0, '.'); _argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
p = argparse.ArgumentParser(); p.add_argument('n', type=int)
p.add_argument('--class', type=int, dest='cls', default=0); p.add_argument('--budget-sec', type=float, default=0, dest='budget')
p.add_argument('--maxexp', type=int, default=0); a = p.parse_args()
def hasA5(quiver, v):
    outs = ap.arrowsOutOf(quiver, v)
    if len(outs) != 1: return False
    ins = [x for x in set(t for (t, h, _) in ap.arrowsInto(quiver, v))]
    for b in ins:
        for c in ins:
            if b < c:
                srcb = {t for (t, h, _) in ap.arrowsInto(quiver, b)}; srcc = {t for (t, h, _) in ap.arrowsInto(quiver, c)}
                if srcb & srcc: return True
    return False
def hasLongSquare(alg, v):
    outs = ap.arrowsOutOf(alg.quiver, v)
    if len(outs) != 1: return False
    e = outs[0][1]
    for rel in alg.rels:
        if len(rel) >= 2 and all(len(q) >= 3 and q[-1] == e and q[-2] == v for q in rel) and len({q[0] for q in rel}) == 1 and len({q[-3] for q in rel}) == len(rel):
            return True
    return False
seen_pv = {}; steps = Counter()
def check(alg, v, base):
    rels = procedure.relationsFrom(alg)
    if not mutation.mutationIsPossibleAtVertex(alg, v): return None
    t = bool(tiltingPlus(alg.quiver, rels, v))
    raw = mutation.quiverMutationAtVertex(alg, v)
    if any(ap.isIllegalRelation(raw.quiver, r) for r in procedure.relationsFrom(raw)): return None
    ch = reduction.reducePathAlgebra(raw)
    guard = search._coxeterKeyOrNone(ch) == base
    pk = fingerprint.canonicalKey(alg) or str(alg.rels)
    steps[t] += 1
    if (pk, v) not in seen_pv:
        outs = ap.arrowsOutOf(alg.quiver, v)
        seen_pv[(pk, v)] = (t, len(outs) == 1, hasA5(alg.quiver, v), hasLongSquare(alg, v))
    return ch, guard
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(a.n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
t0 = time.time(); base = order[a.cls]; seen = set(); frontier = []; nexp = 0
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
stop = False
while frontier and not stop:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if (a.budget and time.time() - t0 > a.budget) or (a.maxexp and nexp >= a.maxexp): stop = True; break
        nexp += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            r = check(alg, v, base)
            if r and r[1] and r[0] is not None: add(r[0])
print('n', a.n, 'class', a.cls, 'expansions', nexp, 'algebras', len(seen), '%.0fs' % (time.time() - t0), 'steps tilting/non', steps[True], steps[False])
for t in (True, False):
    rows = [x for x in seen_pv.values() if x[0] == t]; N = len(rows)
    one = [x for x in rows if x[1]]
    print('tilting' if t else 'REJECTING', 'distinct (parent,v):', N, '| v has one out-arrow:', len(one),
          '| strict A5:', sum(x[2] for x in rows), '| long square:', sum(x[3] for x in rows),
          '| long square among one-out:', sum(x[3] for x in one), '| A5 but not long:', sum(x[2] and not x[3] for x in rows))
