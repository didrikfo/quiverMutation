"""Round 131 scholar: check the radical-power LNAs that Lenzing-Meltzer-Ruan (arXiv:2112.15587v2, Thm
classification-of-Fuchsian-nakayama-cat) identify with Fuchsian singularity categories against Lenzing-de la Pena
(arXiv:math/0611532v1, Table 2): Coxeter polynomial factorisation and the order of Phi.
N_n(r) = kA_n / rad^r; Cartan C[i][j] = 1 iff 0 <= j - i < r; Phi = -C^{-T} C.  Exact arithmetic (sympy).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/131/scholar_lmr.py"""
import sympy as sp

x = sp.symbols('x')

def phi(n, r):
    C = sp.Matrix(n, n, lambda i, j: 1 if 0 <= j - i < r else 0)
    return -C.inv().T * C

def order(M, cap=400):
    n = M.shape[0]; I = sp.eye(n); P = M
    for k in range(1, cap + 1):
        if P == I: return k
        P = P * M
    return None

def cyclo_factor(poly):
    """factorisation as {m: e} for Phi_m^e, or the raw factor list if a non-cyclotomic factor occurs"""
    out = {}
    for f, e in sp.factor_list(sp.Poly(poly, x))[1]:
        hit = None
        for m in range(1, 400):
            if sp.Poly(sp.cyclotomic_poly(m, x), x) == f: hit = m; break
        if hit is None: return ('non-cyclotomic', f.as_expr(), e)
        out[hit] = out.get(hit, 0) + e
    return dict(sorted(out.items()))

# (n, r, LMR type, LdlP Table-2 factorisation of f_A, LdlP period, LdlP mark)
CASES = [
    (11, 4, '<2,5,4] = (2,4,5)', {2: 1, 6: 1, 30: 1}, 30, 'hypersurface'),
    (11, 5, '<2,5,4] = (2,4,5)', {2: 1, 6: 1, 30: 1}, 30, 'hypersurface'),
    (12, 3, '<2,3,7]', {42: 1}, 42, 'hypersurface'),
    (12, 4, '<2,5,5]', {5: 1, 20: 1}, 20, 'hypersurface'),
    (12, 5, '<2,5,5]', {5: 1, 20: 1}, 20, 'hypersurface'),
    (13, 3, '<2,3,8]', {2: 1, 10: 1, 30: 1}, 30, 'hypersurface'),
    (13, 4, '<2,5,6]', {2: 1, 8: 1, 16: 1}, 16, 'hypersurface'),
    (13, 5, '<2,5,6]', {2: 1, 8: 1, 16: 1}, 16, 'hypersurface'),
    (13, 6, '<2,4,7]', {2: 1, 9: 1, 18: 1}, 18, 'hypersurface'),
    (14, 3, '<2,3,9]', {3: 1, 12: 1, 24: 1}, 24, 'hypersurface'),
    (14, 4, '<2,5,7]', {11: 1, 12: 1}, 132, 'complete intersection, 4 generators'),
    (14, 6, '<2,5,7]', {11: 1, 12: 1}, 132, 'complete intersection, 4 generators'),
    (14, 5, '<2,6,6]', {2: 2, 3: 1, 6: 1, 10: 1, 12: 1}, 60, 'complete intersection, 4 generators'),
    (15, 3, '<2,3,10]', {2: 1, 16: 1, 18: 1}, 72, 'complete intersection, 4 generators'),
    (15, 8, '<2,3,10]', {2: 1, 16: 1, 18: 1}, 72, 'complete intersection, 4 generators'),
    (16, 8, '<2,3,11]', None, None, 'complete intersection; LdlP: roots not all on S^1'),
]

if __name__ == '__main__':
    bad = 0
    for n, r, typ, fac, per, mark in CASES:
        M = phi(n, r)
        cp = M.charpoly(x).as_expr()
        f = cyclo_factor(cp)
        o = order(M)
        relstr = str(r) * (n - r) + '0' * (r - 2)
        ok = (f == fac and o == per) if fac is not None else (o is None)
        bad += not ok
        print(f"N_{n}({r}) {relstr:>16}  {typ:<20} f={f}  order={o}  LdlP: {fac} period {per}  [{mark}]  {'OK' if ok else 'MISMATCH'}")
    print('mismatches', bad)
