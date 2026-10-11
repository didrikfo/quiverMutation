"""Round 027 (experimentalist): test conjecture W (E-109) on a fresh guarded BFS walk, parallel-arrow rows included.
Usage: experimentalist_w.py n class budget_sec [max_exp]   (walk as rounds/023/scholar_longsquare.py; counts are cap lower bounds)
Per gate-admitted (algebra, v): K = J != 0 (kerdim), W for out-degree 2 (relations of relationsFrom, arrows keyed so parallel pairs are fine).
Tally key: (outdeg capped 3, parallel-out (two out-arrows with same head), W, K). Mismatches are printed."""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); _a = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec')); sys.argv = _a
n, cls, budget = int(_a[1]), int(_a[2]), float(_a[3]); maxexp = int(_a[4]) if len(_a) > 4 else 0
def W(alg, v):
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v)
    for R in rels:
        if len(R) < 2: continue
        for b1 in outs:
            if not all(len(q) >= 2 and q[-1] == b1 and q[-2][1] == v for q in R): continue
            b2 = [b for b in outs if b != b1][0]
            x = ap.combination({q[:-1]: c for q, c in R.items()})
            if ap.isInIdeal(Q, rels, x): continue
            if ap.isInIdeal(Q, rels, ap.combination({q[:-1] + (b2,): c for q, c in R.items()})): return True
    return False
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; tab = Counter(); mism = []; nalg = 0; stop = False
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
while frontier and not stop:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget or (maxexp and nalg >= maxexp): stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            outs = ap.arrowsOutOf(alg.quiver, v); od = len(outs)
            par = len({b[1] for b in outs}) < od
            kd, _ = kerdim(alg, v, procedure.relationsFrom(alg)); K = bool(kd)
            w = W(alg, v) if od == 2 else None
            tab[(min(od, 3), par, w, K)] += 1
            if od == 2 and w != K and len(mism) < 8: mism.append((v, sorted(alg.quiver.edges(keys=True)), alg.rels, kd))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('n', n, 'class', cls, 'classsize', len(classes[base]), 'expanded', nalg, 'seen', len(seen), 'stopped', stop, '%.0fs' % (time.time() - t0))
print('key = (outdeg capped 3, parallel-out, W (None if outdeg!=2), J!=0): rows')
for k, x in sorted(tab.items(), key=str): print(k, x)
for m in mism: print('MISMATCH', m)
