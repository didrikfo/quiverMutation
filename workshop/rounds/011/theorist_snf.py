"""Is the Z-conjugacy class of the Coxeter matrix (not just its polynomial) what separates the orbits P, Q?
For rows of a line algebra: Cartan C (C_ij = 1 iff path i->j nonzero), Coxeter Phi = -C^T C^{-1}-type matrix;
invariants: Smith normal forms of f(Phi) for f in (Phi-1, Phi+1, Phi^2+Phi+1, Phi^2+1, ...).
  .venv/bin/python workshop/rounds/011/theorist_snf.py N word [word..]    (one line per offset of each word)
"""
import sys
sys.path.insert(0, '.')
import batch
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

def cartan(row):
    n = len(row) + 2
    reach = [n] * (n + 2)
    for i in range(n, 0, -1):
        reach[i] = min(reach[i + 1] if i < n else n, n)
        if i <= len(row) and row[i - 1]:
            reach[i] = min(reach[i], i + row[i - 1] - 1)
    # reach(i) = min over a >= i of (a + r_a - 1); recompute cleanly
    reach = []
    for i in range(1, n + 1):
        m = n
        for a in range(i, n - 1):
            if row[a - 1]: m = min(m, a + row[a - 1] - 1)
        reach.append(m)
    C = sp.zeros(n, n)
    for i in range(1, n + 1):
        for j in range(i, reach[i - 1] + 1): C[i - 1, j - 1] = 1
    return C

def invariants(row):
    C = cartan(row)
    Phi = -C.T * C.inv()
    I = sp.eye(C.shape[0])
    out = []
    for name, M in (('Phi-1', Phi - I), ('Phi+1', Phi + I), ('P2+P+1', Phi**2 + Phi + I), ('P2+1', Phi**2 + I),
                    ('P2-P+1', Phi**2 - Phi + I)):
        d = smith_normal_form(M, domain=sp.ZZ)
        diag = sorted(abs(d[i, i]) for i in range(d.shape[0]))
        out.append((name, tuple(x for x in diag if x != 1)))
    return tuple(out)

if __name__ == '__main__':
    n = int(sys.argv[1])
    for w in sys.argv[2:]:
        for o in range(n):
            row = batch._rowFor(n, w, o)
            if row is None: continue
            print(w, o, invariants(tuple(row)), flush=True)
