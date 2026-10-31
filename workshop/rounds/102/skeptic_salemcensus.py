"""Round 102 skeptic: which Salem numbers lam < 1.5 occur as the Salem factor of char(Phi) for LNAs at n <= 12 (one rep per lnaKeys group),
and are they consistent with published tables (Lehmer poly deg 10; t^10-t^7-t^5-t^3+1 = 1.230391 [Calegari? no: Kellerhals-Kolpakov-type, arXiv 2008.10965 Prop. 1];
2409.11159 table rows 2d=12: 1.302268805, 1.315914432).
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/102/skeptic_salemcensus.py"""
import numpy as np, sympy as sp, collections
from quivermutation import coxeterTables as ct
x = sp.symbols('x')
pub = {"Lehmer deg10 (2409.11159 P0; 2008.10965)": sp.Poly(x**10+x**9-x**7-x**6-x**5-x**4-x**3+x+1, x),
       "t^10-t^7-t^5-t^3+1 (2008.10965 Prop.1)": sp.Poly(x**10-x**7-x**5-x**3+1, x),
       "2409 row 12 / 1.302268805": sp.Poly(x**12-x**11-x**7+x**6-x**5-x+1, x),
       "2409 row 12 / 1.315914432": sp.Poly(x**12-x**8-x**7-x**6-x**5-x**4+1, x)}
for k, P in pub.items():
    r = np.roots([float(c) for c in P.all_coeffs()]); print(k, "largest real root %.9f" % max(r.real[abs(r.imag) < 1e-9]), "reciprocal:", P.all_coeffs() == P.all_coeffs()[::-1])
found = collections.defaultdict(set)
for n in range(4, 13):
    keys = ct.lnaKeys(n); seen = set()
    for rl, kk in keys.items():
        if kk in seen: continue
        seen.add(kk)
        C = sp.Matrix(np.array(ct.lnaCartanMatrix(n, rl)).astype(int).tolist())
        if C.det() != 1: continue
        Ph = -(C.T).inv() * C
        for f, m in sp.factor_list(Ph.charpoly(x).as_expr())[1]:
            P = sp.Poly(f, x); r = np.roots([float(c) for c in P.all_coeffs()])
            big = [z.real for z in r if abs(z.imag) < 1e-9 and z.real > 1.0001]
            if len(big) == 1 and P.degree() >= 4 and np.all(np.abs(np.delete(r, np.argmax(r.real)) ) <= 1 + 1e-7) and big[0] < 1.5:
                found[(round(big[0], 6), tuple(P.all_coeffs()))].add(n)
    print("n", n, "groups", len(seen), flush=True)
for (lam, co), ns in sorted(found.items()): print("lam %.6f deg %d  n in %s  poly %s" % (lam, len(co) - 1, sorted(ns), co))
