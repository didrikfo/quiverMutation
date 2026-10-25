"""Round 087 scholar: the phi(m) misprint of arXiv:1310.1557 (Thm 2(c), l. 227 of the v1 LaTeX): the companion matrix of Phi_m has
multiplicative order m, not phi(m).  Also the order of a block-diagonal sum is lcm of the m, not lcm of phi(m).
usage (repo root): .venv/bin/python workshop/rounds/087/scholar_phi_misprint.py"""
import math, sympy as sp, numpy as np
T = sp.symbols('T')
def order(M):
    I = np.eye(len(M), dtype=object); X = M.copy()
    for k in range(1, 400):
        if (X == I).all(): return k
        X = X.dot(M)
bad = []
for m in range(2, 31):
    c = sp.Poly(sp.cyclotomic_poly(m, T), T).all_coeffs(); d = len(c) - 1
    M = np.zeros((d, d), dtype=object)
    for i in range(d - 1): M[i + 1, i] = 1
    for i in range(d): M[i, d - 1] = -c[d - i]
    o = order(M)
    if o != m: print("ORDER != m at", m, o)
    if o != sp.totient(m): bad.append((m, int(sp.totient(m)), o))
print("m, phi(m), order of companion matrix, where order != phi(m):", bad)
for ms in [(5, 6), (3, 4, 5), (2, 3, 7), (5, 6, 7, 8)]:
    print("M =", ms, "lcm m =", math.lcm(*ms), " lcm phi(m) =", math.lcm(*[int(sp.totient(m)) for m in ms]))
