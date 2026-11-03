"""Round 111 toolsmith, T10: for the key-kept J != 0, D = 0 steps of n = 7 key class 1 (E-201), decide Z-congruence of the Cartan matrices
P C_A P^T = C_B, P in GL_7(Z) (necessary for derived equivalence; E-269 style), independently of tiltingPlus/Hom.
Input: the checkpoint /tmp/skeptic_c1_13_ck1.pkl written by workshop/rounds/069/skeptic_c1_13.py (run to its end first).
Method: invariants (signature, det, SNF of S = C + C^T, charpoly of Phi = -C^-T C), then exact search if S is definite:
rows of P are vectors p_i with p_i S_A p_j^T = S_B[i,j] (diag 2: roots); backtrack, then test P C_A P^T == C_B and det P = +-1.
Controls: (YES) C vs P0 C P0^T for random unimodular P0; (NO) C_A vs a C_B' with one off-diagonal entry changed.
usage (repo root): .venv/bin/python workshop/rounds/111/toolsmith_c1_13_cong.py"""
import sys, pickle, itertools, random
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
import numpy as np, sympy as sp
from cypari import pari
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
from quivermutation import fingerprint
D0 = None
def idet(M): return int(sp.Matrix(M.tolist()).det())
def inv_(C):
    S = C + C.T; M = sp.Matrix(C.tolist()); Phi = -(M.T).inv() * M
    return (idet(S), tuple(np.linalg.eigvalsh(S.astype(float)).round(6) > 0),
            tuple(sorted(int(x) for x in pari.matsnf(pari.matrix(7, 7, [int(v) for v in S.flatten()])))),
            tuple(sp.Poly(Phi.charpoly().as_expr()).all_coeffs()))
def roots(S, b):
    """all integer x with x^T S x = b (S positive definite)"""
    n = S.shape[0]; Sg = pari.matrix(n, n, [int(v) for v in S.flatten()])
    r = pari.qfminim(Sg, b, 10 ** 7, 2)
    M = r[2]; k = int(pari.matsize(M)[1]); out = []
    for j in range(k):
        v = np.array([int(M[i, j]) for i in range(n)], dtype=object); out += [v, -v]
    return out
def congruent(CA, CB, limit=10 ** 7):
    """decide exists P in GL_n(Z): P CA P^T = CB; S definite required. returns (True, P) / (False, None)"""
    n = CA.shape[0]; SA, SB = CA + CA.T, CB + CB.T
    assert (np.linalg.eigvalsh(SA.astype(float)) > 1e-9).all(), 'S not definite'
    cand = roots(SA, 2); cand = [v for v in cand]
    CAo = CA.astype(object)
    # P rows p_i, need p_i CA p_j^T == CB[i,j] for all i,j (ordered, not just symmetric part)
    def rec(rows):
        i = len(rows)
        if i == n:
            P = np.array(rows, dtype=object)
            if abs(round(float(np.linalg.det(P.astype(float))))) == 1: return P
            return None
        for v in cand:
            if v @ CAo @ v != CB[i, i]: continue
            ok = all(rows[j] @ CAo @ v == CB[j, i] and v @ CAo @ rows[j] == CB[i, j] for j in range(i))
            if ok:
                r = rec(rows + [v])
                if r is not None: return r
        return None
    P = rec([])
    return (P is not None), P
if __name__ == '__main__':
    ck = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb')); seen, queue, byCart, byKey, nexp, snaps = ck
    print('checkpoint nexp', nexp, 'queue', len(queue), 'distinct Cartan', len(byCart))
    rows = sorted(byCart.values(), key=lambda r: r['nexp'])
    print('J != 0 steps', len(rows))
    sel = []
    for r in rows:
        dA = tuple(int(x) for x in r['CA'].flatten());
        sel.append(r)
    out = []
    for i, r in enumerate(rows):
        CA, CB = r['CA'], r['CB']
        iA, iB = inv_(CA), inv_(CB)
        d0 = bool(r.get('D0', False)) if 'D0' in r else None
        res = dict(i=i, nexp=r['nexp'], kept=r['kept'], tp=r['tp'], v=r['v'], D0=d0, same_inv=(iA == iB), det=(iA[0], iB[0]), posdef=all(iA[1]))
        if res['same_inv'] and res['posdef']:
            ok, P = congruent(CA, CB); res['congruent'] = ok
        else: res['congruent'] = None
        out.append(res)
        if r['kept'] and not r['tp']: print(res, flush=True)
    pickle.dump(out, open('/tmp/toolsmith_c1_13_cong.pkl', 'wb'))
