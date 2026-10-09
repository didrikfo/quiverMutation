"""H-015 second-invariant audit of every step the Coxeter guard admits.

For every distinct algebra reached (BFS to --depth from every LNA of length n
and its relation dual) and every vertex the gate (procedure.isMutable) admits:
  gate   : the gate admits (necessary condition only)
  guard  : mutated-and-reduced algebra has the same Coxeter key
  tilt   : Ladkani 1001.4765 Prop 2.3(c) (iff) on the parent, linear combinations
  cong   : child Cartan == r+ C r+^T exactly (Prop 3.6), a Z-congruence
Reports counts of the cross-tabulation guard x tilt x cong.
"""
import argparse, sys, time, copy
from fractions import Fraction
from collections import Counter
sys.path.insert(0, '.')
import numpy as np
import networkx as nx
from quivermutation import nakayama as nk, mutation, procedure, reduction, invariants, fingerprint
from quivermutation import arrowPaths as ap, pathAlgebra
from quivermutation import search

def rank(rows):
    rows = [dict(r) for r in rows if r]
    rk = 0; piv = {}
    for r in rows:
        r = {k: Fraction(v) for k, v in r.items()}
        while r:
            h = min(r)
            if h in piv:
                f = r[h]; pr = piv[h]
                r = {k: r.get(k, 0) - f * pr.get(k, 0) for k in set(r) | set(pr)}
                r = {k: v for k, v in r.items() if v != 0}
            else:
                f = r[h]; piv[h] = {k: v / f for k, v in r.items()}; rk += 1; break
    return rk

def tiltingPlus(quiver, rels, k):
    """Prop 2.3(c): p |-> (p beta)_beta injective on paths i~>k mod I, all i != k."""
    outs = ap.arrowsOutOf(quiver, k)
    cache = {}
    def piv(i, j):
        if (i, j) not in cache:
            cache[(i, j)] = ap.idealBasis(quiver, rels, i, j)
        return cache[(i, j)]
    for i in quiver.nodes:
        if i == k: continue
        P = ap.allPathsBetween(quiver, i, k)
        if not P: continue
        dimV = len(P) - len(piv(i, k))
        rows = []
        for p in P:
            row = {}
            for b in outs:
                res = ap.reduceAgainstPivots(ap.combination([p + (b,)]), piv(i, b[1]))
                for kk, v in res.items():
                    row[(b, kk)] = v
            rows.append(row)
        if rank(rows) < dimV:
            return False
    return True

def cartan(alg):
    M = invariants.cartanMatrix(alg, exact=True)
    return np.array(M.tolist(), dtype=object)

def rplus(alg, k, verts):
    n = len(verts); idx = {v: i for i, v in enumerate(verts)}
    R = np.identity(n, dtype=object)
    for j in verts:
        R[idx[k], idx[j]] = -(1 if j == k else 0) + sum(1 for (t, h, _) in ap.arrowsOutOf(alg.quiver, k) if h == j)
    return R

def main():
    a = argparse.ArgumentParser()
    a.add_argument('n', type=int); a.add_argument('--depth', type=int, default=3)
    a.add_argument('--limit', type=int, default=0); a.add_argument('--plan', action='store_true')
    a.add_argument('--unguarded', action='store_true'); a.add_argument('--show', type=int, default=5)
    # overnight.py passes --budget-hours to every job; stop there, print what was
    # tallied, and exit 2 ("out of budget") so it is not restarted.
    a.add_argument('--budget-hours', type=float, default=0, dest='budgetHours')
    a = a.parse_args()
    starts = []
    for lna in nk.LinearNakayamaAlgebra.allOfLength(a.n):
        starts.append(lna); starts.append(pathAlgebra.dualPathAlgebra(lna))
    seen = {}; frontier = []
    def add(alg):
        key = fingerprint.canonicalKey(alg)
        if key is not None:
            if key in seen: return False
            seen[key] = 1
        frontier.append(alg); return True
    for s in starts: add(s)
    tab = Counter(); bad = []; t0 = time.time(); expanded = 0
    outOfBudget = False
    for d in range(a.depth):
        cur, frontier[:] = frontier[:], []
        for alg in cur:
            if a.budgetHours and time.time() - t0 > a.budgetHours * 3600:
                outOfBudget = True; break
            expanded += 1
            baseKey = search._coxeterKeyOrNone(alg)
            if baseKey is None: continue
            if list(nx.simple_cycles(alg.quiver)): continue
            verts = sorted(alg.vertices())
            for v in verts:
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                child = mutation.quiverMutationAtVertex(alg, v)
                illegal = any(ap.isIllegalRelation(child.quiver, r) for r in procedure.relationsFrom(child))
                if illegal: tab['illegal-relation'] += 1; continue
                child = reduction.reducePathAlgebra(child)
                ck = search._coxeterKeyOrNone(child)
                if ck is None: tab['cyclic-child'] += 1; continue
                guard = (ck == baseKey)
                rels = procedure.relationsFrom(alg)
                tilt = tiltingPlus(alg.quiver, rels, v)
                R = rplus(alg, v, verts)
                cong = bool((R.dot(cartan(alg)).dot(R.T) == cartan(child)).all())
                tab[('guard' if guard else 'noguard', 'tilt' if tilt else 'NOTtilt', 'cong' if cong else 'NOTcong')] += 1
                if ((guard and not (tilt and cong)) or (tilt and not guard) or (cong != tilt)) and len(bad) < 50:
                    bad.append((d, alg.rels, v))
                if guard or a.unguarded: add(child)   # walk only where the guard walks
        print('depth', d + 1, 'expanded so far', expanded, 'distinct next', len(frontier), '%.0fs' % (time.time() - t0), flush=True)
        if a.plan and d >= 0: break
        if outOfBudget:
            print('stopped on the budget during depth', d + 1, '-- counts are partial', flush=True); break
    for k, v in sorted(tab.items(), key=str): print(k, v)
    print('suspicious', len(bad))
    for b in bad[:a.show]: print(b)
    if outOfBudget: sys.exit(2)
main()
