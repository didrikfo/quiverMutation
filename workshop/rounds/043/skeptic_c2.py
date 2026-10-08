"""Round 043 (skeptic): search gate-admitted J != 0 steps for one satisfying H1, H2 (C_B = [[Z,0],[e_w^T,1]], Z = C_A off v)
whose D(x) = det(xC_B+C_B^T) - det(xC_A+C_A^T) has x^2 coefficient != 1 (c_2 != 0 candidate), parents NOT restricted to the
class walk: guard-off BFS from seed LNAs (and duals) of length n, canonicalKey dedup.
Own matrix code (sympy-free integer polynomials via numpy interpolation); library used only for algebra mutation/gate/reduction/Cartan.
Usage: skeptic_c2.py n seed_lo seed_hi budget_sec maxdepth   (seeds = sorted-by-key class index range; run from repo root)"""
import sys, time, itertools
from collections import Counter
import numpy as np
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec')); _a = ARGV
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
n, lo, hi, budget, maxdepth = int(_a[1]), int(_a[2]), int(_a[3]), float(_a[4]), int(_a[5]); tguard = len(_a) > 6 and _a[6] == 'tguard'; guard = len(_a) > 6 and _a[6] in ('guard', 'tguard')

def detpoly(C):
    """coefficients (low->high) of det(x C + C^T), exact via integer evaluation at 2n+1 points + Lagrange (fractions)."""
    from fractions import Fraction
    m = C.shape[0]; xs = list(range(-m, m + 1)); ys = []
    for x in xs:
        M = x * C + C.T
        ys.append(int(round(np.linalg.det(M.astype(float)))))
    # solve Vandermonde exactly
    import sympy
    X = sympy.symbols('X'); P = sympy.interpolate(list(zip(xs, ys)), X)
    p = sympy.Poly(P, X); co = [int(p.coeff_monomial(X**k)) for k in range(m + 1)]
    return co

def shape(CA, CB, v):
    """H1/H2 test in either transpose convention. Return (ok, w, orient)."""
    m = CA.shape[0]; idx = [k for k in range(m) if k != v]
    for orient in (0, 1):
        A = CA if orient == 0 else CA.T; B = CB if orient == 0 else CB.T
        if not np.array_equal(A[np.ix_(idx, idx)], B[np.ix_(idx, idx)]): continue
        col = B[:, v]; row = B[v, :]
        if col[v] != 1 or any(col[k] for k in idx): continue
        r = [row[k] for k in idx]
        if sorted(r) == [0] * (len(r) - 1) + [1]:
            return True, idx[r.index(1)], orient
    return False, None, None

def strict(CA, CB, v, w, i, orient):
    """C' := C_B - E_vi has polynomial det(xC'+C'^T) == that of C_A, row v of C' = e_w - e_i, w != i."""
    if w == i: return False
    A = CA if orient == 0 else CA.T; B = CB if orient == 0 else CB.T
    Cp = B.copy(); Cp[v, i] -= 1
    return detpoly(Cp) == detpoly(A)

def orbitrel(CA, v, w, i, orient):
    A = CA if orient == 0 else CA.T
    idx = [k for k in range(CA.shape[0]) if k != v]
    Z = A[np.ix_(idx, idx)].astype(float); F = Z @ np.linalg.inv(Z.T)   # F = Z Z^{-T}
    ew = np.zeros(len(idx)); ew[idx.index(w)] = 1; ei = np.zeros(len(idx)); ei[idx.index(i)] = 1
    for sgn in (1, -1):
        Fm = F if sgn == 1 else np.linalg.inv(F); x = ew.copy()
        for s in range(0, 40):
            if np.allclose(x, ei, atol=1e-6): return sgn * s
            x = Fm @ x
    return None

classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
print('classes', len(order), 'sizes', [len(classes[k]) for k in order][:40], 'using', lo, hi, flush=True)
t0 = time.time(); seen = set(); frontier = []; nalg = 0; stop = False; level = 0
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
for c in order[lo:hi]:
    for s in classes[c]: add(s)
tilt = Counter(); zero = []; cross = Counter(); tab = Counter(); cand = []; nj = 0; ngate = 0; distinct = set(); lowdeg = Counter(); orb = Counter()
while frontier and not stop and level <= maxdepth:
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if time.time() - t0 > budget: stop = True; break
        nalg += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        V = sorted(alg.quiver.nodes)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            ngate += 1
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if sorted(ch.vertices()) != V: continue
            tpv = (not tguard) or bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v))
            J = perI(alg, v)
            if J:
                nj += 1
                CA = np.array(invariants.cartanMatrix(alg), dtype=int); CB = np.array(invariants.cartanMatrix(ch), dtype=int)
                vi = V.index(v); ok, w, orient = shape(CA, CB, vi)
                dk = (CA.tobytes(), CB.tobytes(), vi)
                if dk in distinct:
                    if tpv and ((not guard) or search._coxeterKeyOrNone(ch) == order[lo]): add(ch)
                    continue
                distinct.add(dk)
                dA, dB = detpoly(CA), detpoly(CB); D = [b - a for a, b in zip(dA, dB)]
                low = next((k for k, c in enumerate(D) if c), None)
                lowc = (low, D[low] if low is not None else 0)
                tp = bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v)); tilt[(tp, low is None)] += 1
                tab[(ok, lowc)] += 1
                if low is None and len(zero) < 4:
                    zero.append(dict(v=v, J=J, pe=sorted(alg.quiver.edges(keys=True)), pr=repr(procedure.relationsFrom(alg)),
                                     ce=sorted(ch.quiver.edges(keys=True)), cr=repr(procedure.relationsFrom(ch)), pkey=search._coxeterKeyOrNone(alg), ckey=search._coxeterKeyOrNone(ch), ok=ok, tp=tp))
                    if len(zero) >= 4: print('ZERO', zero[-1], flush=True)
                if ok:
                    i = [V.index(x) for x in J][0]
                    s = orbitrel(CA, vi, w, i, orient); orb[(len(J), s)] += 1; st = strict(CA, CB, vi, w, i, orient); cross[(st, s is not None, lowc)] += 1
                    if lowc != (2, 1) or s is None: cand.append((D, s, sorted(alg.quiver.edges(keys=True)), V, v, J))
            if tpv and ((not guard) or search._coxeterKeyOrNone(ch) == order[lo]): add(ch)
    level += 1
print('n', n, 'seeds', lo, hi, 'expanded', nalg, 'levels', level, 'stopped', stop, 'gate', ngate, 'J!=0 steps', nj, 'distinct', len(distinct), '%.0fs' % (time.time() - t0))
print('(H1H2 shape, (lowest degree, coeff) of D): count'); [print('  ', k, c) for k, c in sorted(tab.items(), key=str)]
print('H1/H2 steps: (|supp J|, orbit exponent s): count'); [print('  ', k, c) for k, c in sorted(orb.items(), key=str)]
print('H1H2 cross (strict, orbit relation found, (low deg, coeff)): count'); [print('  ', k, c) for k, c in sorted(cross.items(), key=str)]
print('(tiltingPlus, D == 0): count', dict(tilt))
import pickle; pickle.dump(zero, open('workshop/rounds/043/skeptic_zero_n%dc%d.pkl' % (n, lo), 'wb'))
print('CANDIDATES (H1H2 with D low term != x^2 or no orbit relation):', len(cand))
for c in cand[:5]: print('  ', c)
