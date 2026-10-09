"""Round 027 (scholar): is J != 0 explained by a PAIR of paths with equal signature (W2), and is W (E-109) the b2-zero special case?
J = {x in e_i A e_v : x b = 0 in A for all arrows b out of v}.  sig(p) = (normal form of p*b in e_i A e_{t(b)})_b.
W2: two paths p1,p2 : i -> v, p1 - p2 != 0 in A, sig(p1) = sig(p2).  Then x = p1 - p2 is in J.  (W of E-109 = W2 with b2-component zero.)
Usage: scholar_pairtest.py --hand | n budget_sec class"""
import sys, time
from collections import Counter
sys.path.insert(0, '.'); MY = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/023/scholar_longsquare.py').read()
pre = src.split("if a.hand:")[0].replace("a = p.parse_args()", "a = p.parse_args([])")
exec(compile(pre, 'ls', 'exec')); sys.argv = MY

def analyse(alg, v):
    """(kerdim, kinds of W2 pairs found, list of per-i data).  kinds: tuple over out-arrows of 'z'(both zero)/'n'(equal nonzero)."""
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v); kinds = set()
    kd, mono = kerdim(alg, v, rels)
    for i in Q.nodes:
        if i == v: continue
        P = ap.allPathsBetween(Q, i, v)
        if not P: continue
        Bv = ap.idealBasis(Q, rels, i, v); NF = {}
        for q in P:
            nf = ap.reduceAgainstPivots(ap.combination([q]), Bv)
            comps = [ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(Q, rels, i, b[1])) for b in outs]
            vec = {(j, k): c for j, d in enumerate(comps) for k, c in d.items()}
            if not nf or not vec: continue
            lead = vec[min(vec, key=str)]          # scale so that sig(p)/lead is canonical: p1 - lam*p2 in J iff same canonical sig
            nsig = tuple(sorted(((str(k), c / lead) for k, c in vec.items())))
            NF.setdefault((nsig, tuple(bool(d) for d in comps)), []).append(frozenset((k, c / lead) for k, c in nf.items()))
        for (nsig, nz), L in NF.items():
            if len(set(L)) >= 2: kinds.add(''.join('n' if z else 'z' for z in nz))
    return kd, mono, kinds

def hand():
    sq = [(1,2),(1,3),(2,4),(3,4)]
    cases = {
     'D cancels: commutes into both b':   (sq+[(4,5),(4,6)], [[[1,2,4,5],[1,3,4,5]],[[1,2,4,6],[1,3,4,6]]], 4),
     'W-type: commute b1, kill both b2':  (sq+[(4,5),(4,6)], [[[1,2,4,5],[1,3,4,5]],[[1,2,4,6]],[[1,3,4,6]]], 4),
     'G three paths (two relations)':     (sq+[(4,5),(1,7),(7,5)], [[[1,2,4,5],[1,7,5]],[[1,3,4,5],[1,7,5]]], 4),
     # chain: 3-term kernel element, no two classes with equal signature: gamma1,gamma2 agree under b1, gamma2,gamma3 agree under b2
     'H chain (3-term kernel)':           ([(1,2),(1,3),(1,7),(2,4),(3,4),(7,4),(4,5),(4,6)],
        [[[1,2,4,5],[1,3,4,5]], [[1,3,4,6],[1,7,4,6]], [[1,2,4,6]], [[1,7,4,5]]], 4),
    }
    for name, (arr, rels, v) in cases.items():
        nodes = sorted({x for e in arr for x in e}); A = build(arr, rels, nodes)
        kd, mono, kinds = analyse(A, v); gate = mutation.mutationIsPossibleAtVertex(A, v)
        r = procedure.relationsFrom(A)
        try:
            procedure.mutateAtVertex(A.quiver, r, v, checkCartan=True); cart = 'congruent'
        except procedure.CartanCongruenceError: cart = 'Cartan FAILS'
        print('%-36s gate=%s kerdim=%d single-path-in-J=%s W2-kinds=%s rewrite: %s' % (name, gate, kd, mono, sorted(kinds), cart))

if MY[1] == '--hand': hand(); sys.exit()
n, budget, cls = int(MY[1]), float(MY[2]), int(MY[3])
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
base = order[cls]; t0 = time.time(); seen = set(); frontier = []; tab = Counter(); ex = {}; nalg = 0
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
        if time.time() - t0 > budget: frontier.clear(); break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            kd, mono, kinds = analyse(alg, v); od = min(len(ap.arrowsOutOf(alg.quiver, v)), 3)
            key = ('outdeg', od, 'J!=0' if kd else 'J=0', 'W2' if kinds else 'noW2', ','.join(sorted(kinds)))
            tab[key] += 1
            if kd and not kinds: ex.setdefault('J!=0 no pair', []).append((v, alg.rels, sorted(alg.quiver.edges(keys=True))))
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == base: add(ch)
print('n', n, 'class', cls, 'algebras', nalg, '%.0fs' % (time.time() - t0))
for k in sorted(tab, key=str): print(k, tab[k])
for m in ex.get('J!=0 no pair', [])[:5]: print('J!=0 but no equal-signature pair:', m)
import pickle; pickle.dump(ex, open('/tmp/ex.pkl','wb'))
