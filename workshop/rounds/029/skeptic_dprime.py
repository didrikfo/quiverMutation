"""Round 029 (skeptic): gate admission of D' rejects, analogues in other n/classes, and a selectivity null.
Usage: skeptic_dprime.py n class_lo class_hi max_exp [budget_sec_per_class]
Per gate-admitted (algebra, v) of the round-023 guarded walk (same as experimentalist_w.py), tally
 key = (outdeg capped 3, L, W, K): K = J != 0 (kerdim); W = E-107 rule (out-degree 2 only);
 L = loose D' shape: some two-term relation with both terms ending b1 through v and x=p1-p2 not in I (commute into one out-arrow),
     and some monomial relation whose last arrow is the other out-arrow b2 with penultimate vertex v (zero relation into b2).
For K rows: witness pattern of the gate: over nonzero paths p into v, which out-arrows have p*arrow nonzero in A."""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); _a = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec')); sys.argv = _a
n, lo, hi, maxexp = int(_a[1]), int(_a[2]), int(_a[3]), int(_a[4]); bud = float(_a[5]) if len(_a) > 5 else 1e9
def shapes(alg, v):
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v); W = Lw = False
    for b1 in outs:
        b2 = [b for b in outs if b != b1][0]
        for R in rels:
            if len(R) < 2 or not all(len(q) >= 2 and q[-1] == b1 and q[-2][1] == v for q in R): continue
            x = ap.combination({q[:-1]: c for q, c in R.items()})
            if ap.isInIdeal(Q, rels, x): continue
            if any(len(M) == 1 and all(len(q) >= 2 and q[-1] == b2 and q[-2][1] == v for q in M) for M in rels): Lw = True
            if ap.isInIdeal(Q, rels, ap.combination({q[:-1] + (b2,): c for q, c in R.items()})): W = True
    return Lw, W
def witness(alg, v):
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v); c = Counter()
    for i in Q.nodes:
        if i == v: continue
        for p in ap.allPathsBetween(Q, i, v):
            if ap.isInIdeal(Q, rels, ap.combination([p])): continue
            w = tuple(sorted(j for j, b in enumerate(outs) if not ap.isInIdeal(Q, rels, ap.combination([p + (b,)]))))
            c[w] += 1
    return c
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
print('n', n, 'classes', len(order), 'sizes', [len(classes[k]) for k in order])
for cls in range(lo, min(hi + 1, len(order))):
    base = order[cls]; t0 = time.time(); seen = set(); frontier = []; tab = Counter(); wit = Counter(); nalg = 0; stop = False; rej = 0
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
            if (maxexp and nalg >= maxexp) or time.time() - t0 > bud: stop = True; break
            nalg += 1
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                od = len(ap.arrowsOutOf(alg.quiver, v))
                kd, _ = kerdim(alg, v, procedure.relationsFrom(alg)); K = bool(kd)
                Lw, W = shapes(alg, v) if od == 2 else (None, None)
                tab[(min(od, 3), Lw, W, K)] += 1
                if od == 2 and K:
                    rej += 1
                    for w, c in witness(alg, v).items(): wit[w] += c
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if search._coxeterKeyOrNone(ch) == base: add(ch)
    print('CLASS', cls, 'size', len(classes[base]), 'expanded', nalg, 'seen', len(seen), 'stopped', stop, '%.0fs' % (time.time() - t0), flush=True)
    print('  key=(outdeg,L,W,K):', {k: x for k, x in sorted(tab.items(), key=str)}, flush=True)
    print('  witness patterns over nonzero paths into v at out-2 rejecting rows (arrow indices with p*arrow != 0):', dict(wit), flush=True)
