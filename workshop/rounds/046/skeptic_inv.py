"""Round 046 (skeptic), T10 (i): congruence invariants of Cartan matrices (C_B = P C_A P^T, P in GL_n(Z) is necessary for derived equivalence of finite
global dimension algebras).  Each is invariant under ANY such P, not only the step's own map, so inequality proves different derived class, equality proves nothing.
invariants(C): Coxeter charpoly; Smith forms of x C + C^T for x in -3..3 (x=1 symmetrised, x=-1 skew part); Smith forms of f(Phi), f^2(Phi) for irreducible
factors f of the charpoly, Phi = -C^{-T} C; signature of C + C^T; histogram of q(v) = v^T C v mod m, m = 2..7 on (Z/m)^n.
Usage: skeptic_inv.py in.pkl [power]     (power: instead test separation of LNAs within key class)"""
import sys, pickle, itertools
import numpy as np, sympy
from sympy import Matrix, ZZ, symbols, Poly, factor_list
from sympy.matrices.normalforms import smith_normal_form
X = symbols('X')
def snf(M):
    S = smith_normal_form(M, domain=ZZ); n = min(S.shape)
    return tuple(sorted(abs(int(S[i, i])) for i in range(n)))
def qhist(C, m):
    n = C.shape[0]; Cn = np.array(C.tolist(), dtype=np.int64)
    # enumerate (Z/m)^n in chunks
    grids = np.indices((m,) * n).reshape(n, -1).T.astype(np.int64)
    q = np.einsum('ij,jk,ik->i', grids, Cn, grids) % m
    return tuple(np.bincount(q, minlength=m).tolist())
def invariants(Cl, mmax=7, deep=True):
    C = Matrix(Cl); n = C.shape[0]; I = sympy.eye(n)
    Phi = -(C.T.inv()) * C
    cp = Poly(Phi.charpoly(X).as_expr(), X)
    inv = {'charpoly': tuple(cp.all_coeffs())}
    for x in range(-3, 4): inv['snf(%dC+C^T)' % x] = snf(x * C + C.T)
    for f, e in factor_list(cp.as_expr(), X)[1]:
        fp = Poly(f, X); M = sympy.zeros(n, n)
        for k, c in enumerate(reversed(fp.all_coeffs())): M += c * Phi ** k
        inv['snf f(Phi) f=%s' % (tuple(fp.all_coeffs()),)] = snf(M)
        inv['snf f^2(Phi) f=%s' % (tuple(fp.all_coeffs()),)] = snf(M * M)
    ev = np.linalg.eigvalsh(np.array((C + C.T).tolist(), dtype=float))
    inv['sig'] = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()), int((abs(ev) <= 1e-9).sum()))
    if deep:
        for m in range(2, mmax + 1): inv['q mod %d' % m] = qhist(C, m)
    return inv
def diff(a, b):
    return [k for k in a if a[k] != b.get(k)]
if __name__ == '__main__':
    recs = pickle.load(open(sys.argv[1], 'rb'))
    for kind in ('ctrl', 'fail'):
        rs = [r for r in recs if r['kind'] == kind]; cnt = 0; names = {}
        for r in rs:
            a, b = invariants(r['CA']), invariants(r['CB']); d = diff(a, b)
            if d: cnt += 1
            for k in d: names[k] = names.get(k, 0) + 1
            r['diff'] = d
        print(sys.argv[1], kind, 'steps', len(rs), 'invariant differs in', cnt, names)
    pickle.dump(recs, open(sys.argv[1] + '.inv', 'wb'))
