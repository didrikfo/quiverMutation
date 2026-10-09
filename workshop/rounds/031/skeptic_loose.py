"""Round 031 (skeptic): classify loose-shape (L) gate-admitted out-degree-2 rows by which W clause fails.
Usage: skeptic_loose.py n class max_exp
For each L row (relation R=p1*b1 - p2*b1 through v, x=p1-p2 not in I, monomial M ending (r, b2) with r ending at v),
record pattern (p1*b2 in I?, p2*b2 in I?, (p1-p2)*b2 in I = W-clause, K=J!=0, M-suffix-of p1?, suffix-of p2?, len p1, len p2, len M)."""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); _a = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec')); sys.argv = _a
n, cls, maxexp = int(_a[1]), int(_a[2]), int(_a[3])
def rows(alg, v):
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v); out = []
    for b1 in outs:
        b2 = [b for b in outs if b != b1][0]
        for R in rels:
            if len(R) < 2 or not all(len(q) >= 2 and q[-1] == b1 and q[-2][1] == v for q in R): continue
            x = ap.combination({q[:-1]: c for q, c in R.items()})
            if ap.isInIdeal(Q, rels, x): continue
            Ms = [M for M in rels if len(M) == 1 and all(len(q) >= 2 and q[-1] == b2 and q[-2][1] == v for q in M)]
            if not Ms: continue
            ps = [q[:-1] for q in R]
            ind = tuple(ap.isInIdeal(Q, rels, ap.combination([p + (b2,)])) for p in ps)
            W = ap.isInIdeal(Q, rels, ap.combination({q[:-1] + (b2,): c for q, c in R.items()}))
            Mq = [list(M)[0] for M in Ms]
            suf = tuple(tuple(any(p[len(p)-(len(m)-1):] == m[:-1] for m in Mq) for p in ps) for _ in [0])[0]
            out.append((ind, W, suf, tuple(len(p) for p in ps), tuple(len(m) for m in Mq), str(R), str(Mq)))
    return out
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; seen = set(); frontier = []; tab = Counter(); ex = {}; nalg = 0
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
while frontier and nalg < maxexp:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if nalg >= maxexp: break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            if len(ap.arrowsOutOf(alg.quiver, v)) == 2:
                kd, _ = kerdim(alg, v, procedure.relationsFrom(alg)); K = bool(kd)
                for r in rows(alg, v):
                    key = (K,) + r[:5]; tab[key] += 1; ex.setdefault(key, (r[5], r[6], str(alg.quiver.edges) if hasattr(alg.quiver,'edges') else ''))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('n', n, 'class', cls, 'expanded', nalg)
print('key = (K, (p1*b2 in I, p2*b2 in I), W-clause(diff in I), (M suffix of p1, of p2), (len p1, len p2), len M)')
for k, c in sorted(tab.items(), key=str): print(c, k, '\n    e.g. R=', ex[k][0], ' M=', ex[k][1])
