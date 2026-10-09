"""Round 043 (skeptic): random acyclic parents (m vertices) -- every gate-admitted J != 0 step, tally H1/H2 shape, lowest term of
D(x) = det(xC_B+C_B^T)-det(xC_A+C_A^T), orbit exponent s (e_i = F^s e_w). Random generator as rounds/041/experimentalist_control.py (re-typed);
helpers detpoly/shape/orbitrel are shared with skeptic_c2.py.   Usage: skeptic_rand.py m seed budget_sec [pEdge]"""
import sys, time, random, itertools
from collections import Counter
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
h = open('workshop/rounds/043/skeptic_c2.py').read(); exec(h[h.index('def detpoly'):h.index('classes = {}')])
m, seed, budget = int(_a[1]), int(_a[2]), float(_a[3]); pE = float(_a[4]) if len(_a) > 4 else 0.35
rng = random.Random(seed); t0 = time.time(); N = list(range(1, m + 1))
tab = Counter(); orb = Counter(); cross = Counter(); cand = []; tried = 0; built = 0; nj = 0; dist = set()
def paths(arr, s, t):
    out = []
    def go(p):
        if p[-1] == t and len(p) > 1: out.append(tuple(p))
        for (a, b) in arr:
            if a == p[-1] and b <= t: go(p + [b])
    go([s]); return out
while time.time() - t0 < budget:
    tried += 1
    arr = [(a, b) for a in N for b in N if a < b and rng.random() < pE]
    rels = []
    for s in N:
        for t in N:
            if s >= t: continue
            P = [p for p in paths(arr, s, t) if len(p) >= 3]
            if len(P) >= 2 and rng.random() < 0.7:
                k = rng.choice([2, len(P)]) if len(P) > 2 else 2
                rels.append([list(p) for p in rng.sample(P, k)])
    for (a, b) in arr:
        for (b2, c) in arr:
            if b2 == b and rng.random() < 0.15: rels.append([[a, b, c]])
    try:
        A = build(arr, rels, N)
        if any(ap.isIllegalRelation(A.quiver, rr) for rr in procedure.relationsFrom(A)): continue
    except Exception: continue
    built += 1
    for v in N:
        try:
            if not mutation.mutationIsPossibleAtVertex(A, v): continue
            J = perI(A, v)
            if not J: continue
            raw = mutation.quiverMutationAtVertex(A, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != N: continue
            CA = np.array(invariants.cartanMatrix(A), dtype=int); CB = np.array(invariants.cartanMatrix(ch), dtype=int)
        except Exception: continue
        vi = N.index(v); dk = (CA.tobytes(), CB.tobytes(), vi)
        if dk in dist: continue
        dist.add(dk); nj += 1
        ok, w, orient = shape(CA, CB, vi)
        dA, dB = detpoly(CA), detpoly(CB); D = [b - a for a, b in zip(dA, dB)]
        low = next((k for k, c in enumerate(D) if c), None); lowc = (low, D[low] if low is not None else 0)
        tab[(ok, lowc)] += 1
        if ok:
            i = N.index(list(J)[0]); s = orbitrel(CA, vi, w, i, orient); orb[(len(J), s)] += 1; cross[(strict(CA, CB, vi, w, N.index(list(J)[0]), orient), s is not None, lowc)] += 1
            if lowc != (2, 1) and s is not None and cross and strict(CA, CB, vi, w, i, orient): cand.append((D, 's=', s, arr, rels, v, J))
print('m', m, 'seed', seed, 'tried', tried, 'built', built, 'distinct J!=0', nj, '%.0fs' % (time.time() - t0))
print('(H1H2, (low deg, coeff)): count'); [print('  ', k, c) for k, c in sorted(tab.items(), key=str)]
print('H1H2 (|supp J|, s): count'); [print('  ', k, c) for k, c in sorted(orb.items(), key=str)]
print('H1H2 cross (|supp J|, orbit relation found, (low deg, coeff)): count'); [print('  ', k, c) for k, c in sorted(cross.items(), key=str)]
print('CANDIDATES', len(cand)); [print('  ', c) for c in cand[:6]]
