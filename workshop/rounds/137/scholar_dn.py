"""Re-derive 'Lambda(n, n-1) = kA_n/rad^(n-1) is derived equivalent to kD_n' from
Bruestle 2001 Thm 1.1 + Thm 1.2 (Table I, corank 0 row): a tree algebra with
positive definite Euler form is derived-tame (1.1) and then derived equivalent to
k Delta, Delta Dynkin, determined by n, corank and Dynkin type (1.2); read from the
PDF, research/literature/bruestle-derived-tame-tree-algebras.md l. 32-41.
General n (R137 answer to referee, argument by maverick): C^-1 + C^-T is the
signed n-cycle matrix 2I - (path adjacency) + (e_1 e_n^T + e_n e_1^T), and
2q(x) = (x_1 + x_n)^2 + sum_{i<n} (x_i - x_{i+1})^2 is PD; det = 4 = det of D_n's
Cartan matrix only (A_n: n+1, E_6/7/8: 3/2/1). This script checks per n:
(a) the signed-cycle identity and det = 4 exactly, (b) PD by leading minors,
(c) the Coxeter polynomial is D_n's and not A_n's or E_n's.
usage: .venv/bin/python -I workshop/rounds/137/scholar_dn.py [NMAX]   (from repo root)
Controls: N_9(3) (affine E8~, form only semidefinite -> NOT PD), N_12(3) (not PH,
NOT PD), N_7(3) (E_7: PD, polynomial E_7 not D_7), N_8(3) (E_8)."""
import sys
sys.path.insert(0, '.')
import networkx as nx
import sympy
from quivermutation.coxeterTables import lnaCartanMatrix
from quivermutation.invariants import coxeterCoefficients
from quivermutation.treeSearch import treeCoxeterKey


def rels(n, r):
    return [r] * (n - r) + [0] * (r - 2)


def star(a, b, c):
    g = nx.Graph(); g.add_node(0); k = 1
    for arm in (a, b, c):
        prev = 0
        for _ in range(arm - 1):
            g.add_edge(prev, k); prev = k; k += 1
    return g


def path(n):
    g = nx.path_graph(n)
    return g


def report(n, r):
    C = sympy.Matrix(lnaCartanMatrix(n, rels(n, r)))
    E = C.inv()
    S = (E + E.T) / 2
    if r == n - 1:
        cyc = sympy.zeros(n, n)
        for i in range(n):
            cyc[i, i] = 2
        for i in range(n - 1):
            cyc[i, i + 1] = cyc[i + 1, i] = -1
        cyc[0, n - 1] = cyc[n - 1, 0] = 1
        assert E + E.T == cyc and cyc.det() == 4, ('signed-cycle identity fails', n)
    pd = all(S[:k, :k].det() > 0 for k in range(1, n + 1))
    key = tuple(coxeterCoefficients(lnaCartanMatrix(n, rels(n, r))))
    names = {'A': path(n)}
    if n >= 4:
        names['D'] = star(2, 2, n - 2)
    if n in (6, 7, 8):
        names['E'] = star(2, 3, n - 3)
    match = [nm for nm, g in names.items() if tuple(treeCoxeterKey(g)) == key]
    return pd, match


nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
bad = 0
for n in range(4, nmax + 1):
    pd, match = report(n, n - 1)
    ok = pd and match == ['D']
    bad += not ok
    print(f'Lambda({n},{n-1}): signed-cycle form, det 4; PD={pd} Coxeter matches {match}  {"OK" if ok else "FAIL"}')
print('controls:')
for n, r, want in [(9, 3, (False, [])), (12, 3, (False, [])), (7, 3, (True, ['E'])), (8, 3, (True, ['E']))]:
    got = report(n, r)
    print(f'  N_{n}({r}): PD={got[0]} matches {got[1]}  expected {want}  {"OK" if got == want else "FAIL"}')
print('failures among Lambda(n,n-1):', bad)
