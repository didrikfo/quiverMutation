"""Round 025 (scholar): the n = 8 class-0 gate-admitted rejects with out-degree >= 2 (found by the r023 referee, 42 steps).
Guarded BFS as rounds/023/scholar_longsquare.py. At each gate-admitted step with J != 0 (sum dim ker g_i > 0) and out-degree >= 2: record distinct (parent,v)
by canonical key; classify: 'D' = for every arrow out of v some relation (>= 2 paths) has all paths ending ..,v,e_beta (commutes into each arrow);
'D-part' = into some arrows but not all; 'none' = no such relation on alg.rels (G-like: presentation differs from generators). Then run the actual rewrite
(quiverMutationAtVertex + reducePathAlgebra) and test Cartan(child) == R C R^T (rplus as rounds/018). Prints the first reject in full.
Usage: scholar_n8rejects.py n [--class I] [--budget-sec S]"""
import argparse, sys, time
from collections import Counter
sys.path.insert(0, '.'); _argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
p = argparse.ArgumentParser(); p.add_argument('n', type=int)
p.add_argument('--class', type=int, dest='cls', default=0); p.add_argument('--budget-sec', type=float, default=0, dest='budget'); a = p.parse_args()

def kerdim(alg, v, rels):
    quiver = alg.quiver; outs = ap.arrowsOutOf(quiver, v); tot = 0; mono = False
    for i in quiver.nodes:
        if i == v: continue
        P = ap.allPathsBetween(quiver, i, v)
        if not P: continue
        dimV = len(P) - len(ap.idealBasis(quiver, rels, i, v)); rows = []
        for q in P:
            row = {}
            for b in outs:
                for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(quiver, rels, i, b[1])).items(): row[(b, kk)] = x
            rows.append(row)
            if not row and not ap.isInIdeal(quiver, rels, ap.combination([q])): mono = True
        tot += dimV - (rank(rows) if rows else 0)
    return tot, mono
def intoArrow(alg, v, b):
    e = b[1]
    return any(len(rel) >= 2 and all(len(q) >= 3 and q[-1] == e and q[-2] == v for q in rel) for rel in alg.rels)
def classify(alg, v):
    outs = ap.arrowsOutOf(alg.quiver, v); k = sum(intoArrow(alg, v, b) for b in outs)
    return 'D' if k == len(outs) else ('D-part%d/%d' % (k, len(outs)) if k else 'none')
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(a.n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
t0 = time.time(); base = order[a.cls]; seen = set(); frontier = []; tab = Counter(); rej = {}; steps = Counter()
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
while frontier:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if a.budget and time.time() - t0 > a.budget: frontier.clear(); break
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(alg.quiver, v)
            kd, mono = kerdim(alg, v, rels)
            if kd and len(outs) >= 2:
                steps['outdeg>=2 J!=0 steps'] += 1
                key = (fingerprint.canonicalKey(alg) or str(alg.rels), v)
                if key not in rej: rej[key] = (alg, v, kd, classify(alg, v))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('n', a.n, 'class', a.cls, 'algebras', len(seen), '%.0fs' % (time.time() - t0), dict(steps), 'distinct (parent,v):', len(rej))
shown = False
for key, (alg, v, kd, cl) in rej.items():
    rels = procedure.relationsFrom(alg); verts = sorted(alg.vertices()); tp = bool(tiltingPlus(alg.quiver, rels, v))
    raw = mutation.quiverMutationAtVertex(alg, v); ok = not any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw))
    cong = 'illegal'
    if ok:
        ch = reduction.reducePathAlgebra(raw)
        if len(ch.vertices()) != len(verts): cong = 'dim'
        else: cong = bool((rplus(alg, v, verts).dot(cartan(alg)).dot(rplus(alg, v, verts).T) == cartan(ch)).all())
    tab[(cl, 'outdeg', len(ap.arrowsOutOf(alg.quiver, v)), 'dimJ', kd, 'tiltingPlus', tp, 'cartanCong', cong)] += 1
    if not shown:
        shown = True; print('FIRST REJECT: v =', v, 'arrows', sorted(set((x, y) for x, y, _ in alg.quiver.edges(keys=True))) if hasattr(alg.quiver, 'edges') else alg.quiver)
        print(' rels =', alg.rels); print(' class', cl, 'dimJ', kd, 'tiltingPlus', tp, 'cartanCong', cong)
for k, x in sorted(tab.items(), key=str): print(k, x)
