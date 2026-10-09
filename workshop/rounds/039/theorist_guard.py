"""Round 039 (theorist; adapted from rounds/038/experimentalist_review_cartan.py): for H != 0 steps record guard pass, det C_A, det C_B, and test the determinant lemma
Original: Round 038 referee (experimentalist): test C_B = r C_A r^T + H on walk algebras (n, class, budget s).
H_{v i} = dim J_i (perI), all else 0.  B = reducePathAlgebra(quiverMutationAtVertex) as on the walk.
Usage: experimentalist_review_cartan.py n cls budget_sec"""
import sys, time
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
n, cls, budget = int(_a[1]), int(_a[2]), float(_a[3])
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for s in classes[base]: add(s)
from collections import Counter
tab = Counter(); tab2 = Counter(); bad = []; vsetbad = 0
while frontier and not stop:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget: stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        Q = alg.quiver; V = sorted(Q.nodes)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: vsetbad += 1; continue
            CA = np.array(invariants.cartanMatrix(alg), dtype=int); CB = np.array(invariants.cartanMatrix(ch), dtype=int)
            k = V.index(v); r = np.eye(n, dtype=int)
            for j, w in enumerate(V):
                r[k, j] = -(j == k) + sum(1 for (a, b) in Q.edges() if a == v and b == w)
            # edges() on a MultiDiGraph yields one per parallel arrow
            H = np.zeros((n, n), dtype=int)
            for i, jj in perI(alg, v).items(): H[k, V.index(i)] = jj
            Hn = H.any()
            if Hn:
                Cp = r @ CA @ r.T; jv = H[k].copy(); ev = np.zeros(n, dtype=int); ev[k] = 1
                passed = search._coxeterKeyOrNone(ch) == base
                dA, dB = round(np.linalg.det(CA)), round(np.linalg.det(CB))
                t = float(jv @ np.linalg.inv(Cp) @ ev)          # det CB / det C' = 1 + t
                # determinant lemma at x = 2, 3, 5:  det(xCB+CB^T)/det(xC'+C'^T) = (1+x a)(1+e) - x b c
                lem = True
                for x in (2.0, 3.0, 5.0):
                    N = np.linalg.inv(x * Cp + Cp.T); a_ = jv @ N @ ev; e_ = ev @ N @ jv; b_ = jv @ N @ jv; c_ = ev @ N @ ev
                    lhs = np.linalg.det(x * CB + CB.T) / np.linalg.det(x * Cp + Cp.T)
                    lem &= abs(lhs - ((1 + x * a_) * (1 + e_) - x * b_ * c_)) < 1e-6
                tab2[(passed, dA, dB, int(jv.sum()), lem)] += 1
            if search._coxeterKeyOrNone(ch) == base: add(ch)
tab2_print = True
print('expanded', nalg, 'vertex-set-mismatch', vsetbad, 'secs %.0f' % (time.time()-t0))
print('(passed, detCA, detCB, sum dimJ, lemma ok): count'); [print('  ', k, c) for k, c in sorted(tab2.items(), key=str)]

