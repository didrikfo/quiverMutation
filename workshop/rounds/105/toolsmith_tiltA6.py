"""Round 105 (toolsmith), H-027 / T10: independent tilting test applied to every one-step mutation of B_6 (and controls).
Test = workshop/rounds/050/skeptic_tilt.py (chain-map Hom in K^b(proj), no tiltingPlus / perI / gate / E-187 rule for the verdict).
Per (algebra a, vertex k): T = (+_{i!=k} P_i) + T_k, T_k = (P_k -> +_{k->h} P_h).  Verdicts:
  hom   : Hom(T,T[-1]) = 0 and Hom(T,T[1]) = 0
  count : T has |V| indecomposable summands (by construction), generation by construction (P_k = cone-shift of T_k over P_h, h != k, acyclic)
  det   : det H == det Cartan(a)  (K0: a necessary condition for generation, independent of the Hom vanishing)
  tilt  := hom and det.
Then: is H (Cartan of End(T)) equal to Cartan(A_6) up to simultaneous permutation, or its transpose?  (necessary for End(T) ~ A_6.)
Left mutation = same test on a.opposite().  Usage: toolsmith_tiltA6.py [n]   (n = 6 default; n = 5 gives B_5 vs A_5)"""
import sys, itertools
sys.path.insert(0, '.')
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import numpy as np
from quivermutation import nakayama as nk, invariants, procedure, mutation, reduction
from quivermutation.procedure import mutateAtVertex

def LNA(n, row): return nk.LinearNakayamaAlgebra(n, list(row))
def cart(a): return np.array(invariants.cartanMatrix(a, exact=True).tolist(), dtype=int)

def test(a, k, S=-1):
    V, r = stepTest(a, k, S)
    H = np.array(r[0], dtype=int)
    hom = (not np.array(r[-1]).any()) and (not np.array(r[1]).any())
    detok = round(abs(np.linalg.det(H))) == round(abs(np.linalg.det(cart(a))))
    return H, hom, detok, int(np.abs(np.array(r[-1])).sum()), int(np.abs(np.array(r[1])).sum())

def permEq(X, Y):
    n = len(X)
    if sorted(X.flatten()) != sorted(Y.flatten()): return False
    for p in itertools.permutations(range(n)):
        if (X[np.ix_(p, p)] == Y).all(): return True
    return False
def matchesA(H, CA): return permEq(H, CA) or permEq(H.T, CA)

def repoChild(a, k, left=False):
    """context only: the repo's (reduced) child Cartan or the refusal."""
    b = a.opposite() if left else a
    try:
        return 'gate %s' % ('admits' if procedure.isMutable(b.quiver, procedure.relationsFrom(b), k) else 'REFUSES')
    except Exception as e:
        return 'gate error (%s)' % type(e).__name__

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6   # n = 6 (default) or 5: B_n = lone 3-arrow end relation, A_n = rad^(n-1)
    A6 = LNA(n, (n - 1,) + (0,) * (n - 3)); B6 = LNA(n, (0,) * (n - 4) + (3, 0)); CA = cart(A6)
    print('n = %d  A row %s  B row %s  det Cartan A %d B %d' % (n, (n - 1,) + (0,) * (n - 3), (0,) * (n - 4) + (3, 0), round(np.linalg.det(CA)), round(np.linalg.det(cart(B6)))))
    print('== CONTROLS ==')
    # known tilting: A_3 (no relations) and LNA steps; test output vs the repo's Cartan of the child
    for nm, a in (('A3 hereditary', LNA(3, (0,))), ('LNA 4 (2,0)', LNA(4, (2, 0))), ('B6', B6)):
        for k in range(1, a.length):
            H, hom, detok, m1, p1 = test(a, k)
            print(nm, 'k', k, 'hom0', hom, 'det', detok, '(Hom[-1],Hom[1] mass %d,%d)' % (m1, p1), repoChild(a, k))
    # known non-tilting 1: the planted error, T_k with P_k -> sum P_h replaced by the wrong direction is not built here; instead use
    # shift-mixed T: summand T_k placed in degrees (S,S+1) but the others in degree 1 (planted wrong shift)
    print('-- planted error: A3 hereditary, P_1 stays in degree 0 but P_2 moved to degree 1 (T = P_1 + P_2[-1] + P_3: not tilting)')
    a = LNA(3, (0,)); A = Alg(a)
    T = {1: ({0: [1]}, {}), 2: ({1: [2]}, {}), 3: ({0: [3]}, {})}
    bad = {(i, j, m): homdim(A, T[i], T[j], m) for i in T for j in T for m in (-1, 1)}
    print('   nonzero Hom(T,T[+-1]) entries:', sum(1 for v in bad.values() if v))
    print('-- planted error 2: T = P_1 + P_1 + P_3 (repeated summand, count wrong / det 0)')
    H = np.array([[homdim(A, ({0:[x]},{}), ({0:[y]},{}), 0) for y in (1,1,3)] for x in (1,1,3)])
    print('   det H =', round(np.linalg.det(H)), ' vs det Cartan', round(np.linalg.det(cart(a))))
    print('== B6 -> A6: every vertex, right and left ==')
    out = []
    for left in (False, True):
        a = B6.opposite() if left else B6
        for k in sorted(a.quiver.nodes):
            if not list(a.quiver.successors(k)): print('left' if left else 'right', 'k', k, 'no out-arrows (T_k = P_k[?]) skipped'); continue
            H, hom, detok, m1, p1 = test(a, k)
            m = matchesA(H, CA)
            print('left ' if left else 'right', 'k', k, 'hom0', hom, 'det', detok, 'tilting', hom and detok, 'mass', (m1, p1),
                  'Cartan(End T) == Cartan(A6) up to perm/transp:', m, '|', repoChild(B6, k, left))
            print('     H rows', H.tolist())

    print('== AGREEMENT with the gate over all LNAs, right steps, n = 4, 5, 6 (power: both verdicts occur) ==')
    from collections import Counter
    for nn in (4, 5, 6):
        c = Counter()
        for lna in nk.LinearNakayamaAlgebra.allOfLength(nn):
            for k in range(1, nn):
                H, hom, detok, m1, p1 = test(lna, k)
                g = procedure.isMutable(lna.quiver, procedure.relationsFrom(lna), k)
                c[(g, hom and detok)] += 1
        print('n', nn, '(gate admits, indep tilting) : count', dict(c))
