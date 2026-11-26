"""Round 139 (scholar): find an explicit isometry of (non-symmetric) Euler forms
K0(Lambda(n,n-1)) ~ K0(kD_n), i.e. P in GL_n(Z) with P^T E_D P = E_L, E = C^{-1}.
With Happel 1987 Thm 5.12 as quoted in Hubery-Krause arXiv:1310.1907v2 Thm 7.4
(basic, simply connected A over alg. closed k: A ~der kQ, Q Dynkin, iff the
lattices are isometric), this is a third route to Lambda(n,n-1) ~ kD_n.
Method: backtracking for v_1..v_n in the root set {x : x^T E_L x = 1} (finite: the
symmetric form is PD, E-352) with v_i^T E_L v_j = E_D[i][j]; P^-1 = [v_1..v_n]
must be unimodular. Same Cartan convention as lnaCartanMatrix (M[t][s] = 1 iff a
path s -> t). Controls (must find NO isometry to D_n): A_n-type N_n(2) (rad^2 = 0, ~der kA_n), E_n = N_n(n-2)
at n = 6, 7, and N_n(n-1) against kA_n.
usage (from repo root): .venv/bin/python -I workshop/rounds/139/scholar_isometry.py [NMAX]"""
import itertools
import sys
sys.path.insert(0, '.')
import sympy
from quivermutation.coxeterTables import lnaCartanMatrix


def rels(n, r):
    return [r] * (n - r) + [0] * (r - 2)


def hereditary(n, kind):
    """Cartan matrix of kQ: 'D' arrows 0->2, 1->2, 2->3, ..., n-2->n-1; 'A' linear."""
    arrows = [(0, 2), (1, 2)] + [(i, i + 1) for i in range(2, n - 1)] if kind == 'D' \
        else [(i, i + 1) for i in range(n - 1)]
    M = sympy.eye(n)
    reach = {i: {i} for i in range(n)}
    for _ in range(n):
        for s, t in arrows:
            for i in range(n):
                if s in reach[i]:
                    reach[i].add(t)
    for s in range(n):
        for t in reach[s]:
            M[t, s] = 1
    return M


def roots(E, n, bound):
    S = (E + E.T) / 2
    out = []
    for x in itertools.product(range(-bound, bound + 1), repeat=n):
        v = sympy.Matrix(x)
        if any(x) and (v.T * S * v)[0] == 1:
            out.append(v)
    return out


def isometry(EL, ED, R):
    n = ED.shape[0]
    pick = []

    def rec(i):
        if i == n:
            V = sympy.Matrix.hstack(*pick)
            return V if abs(V.det()) == 1 else None
        for v in R:
            if all((pick[j].T * EL * v)[0] == ED[j, i] and (v.T * EL * pick[j])[0] == ED[i, j]
                   for j in range(i)) and (v.T * EL * v)[0] == ED[i, i]:
                pick.append(v)
                got = rec(i + 1)
                if got is not None:
                    return got
                pick.pop()
        return None
    return rec(0)


def test(CL, CH, bound=2):
    EL, ED = CL.inv(), CH.inv()
    n = CL.shape[0]
    R = roots(EL, n, bound)
    V = isometry(EL, ED, R)
    if V is None:
        return None, len(R)
    assert V.T * EL * V == ED
    return V, len(R)


nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
for n in range(4, nmax + 1):
    CL = sympy.Matrix(lnaCartanMatrix(n, rels(n, n - 1)))
    V, nr = test(CL, hereditary(n, 'D'))
    print(f'N_{n}({n-1}) vs kD_{n}: roots {nr} (D_n has {2*n*(n-1)}); isometry {"FOUND" if V is not None else "NONE"}')
    if V is not None:
        print('   P^-1 columns =', [list(V[:, i]) for i in range(n)])
print('controls (must be NONE):')
for n, r, kind in [(5, 2, 'D'), (6, 2, 'D'), (6, 4, 'D'), (7, 5, 'D'), (5, 4, 'A'), (6, 5, 'A')]:
    CL = sympy.Matrix(lnaCartanMatrix(n, rels(n, r)))
    V, nr = test(CL, hereditary(n, kind))
    print(f'  N_{n}({r}) vs k{kind}_{n}: roots {nr}; isometry {"FOUND (control failed)" if V is not None else "NONE, OK"}')
