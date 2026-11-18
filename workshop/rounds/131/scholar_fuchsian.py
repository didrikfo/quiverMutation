"""Round 131 scholar: Coxeter polynomials of extended canonical algebras (genus-0 Fuchsian singularity categories)
against LNAs.
(1) phi_T(p) = prod v_{a_i} * [ (1 - 2x - 2x^2 + x^3) + x^2 * sum v_{a_i - 1} / v_{a_i} ]  (g = 0), read off from
    arXiv:1909.10362v1 Prop 3.8 (p_R times phi_X; the printed bracket there is checked here, not trusted);
    controlled on the LMR identifications N_n(r) ~ <a,b,c] (2112.15587v2 Thm 6.3) and on LdlP Table 2 (math/0611532v1).
(2) every LNA at n = 8..MAXN whose Coxeter polynomial equals phi_T of a Table-2 weight sequence (periodic, rho = 1).
Coxeter polynomial equality is necessary for derived equivalence, not sufficient: (2) lists candidates only.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/131/scholar_fuchsian.py [MAXN [MINN]] [--ci-only]"""
import sys, os, itertools, collections
import numpy as np
import sympy as sp
sys.path.insert(0, os.getcwd())
from quivermutation.nakayama import allRelationLengths

x = sp.symbols('x')
v = lambda a: sum(x**i for i in range(a))

def phiT(p):
    expr = sp.prod([v(a) for a in p]) * ((1 - 2*x - 2*x**2 + x**3) + x**2 * sum(v(a - 1) / v(a) for a in p))
    return sp.Poly(sp.cancel(expr), x)

def cyclo(poly):
    out = {}
    for f, e in sp.factor_list(poly)[1]:
        m = next((m for m in range(1, 300) if sp.Poly(sp.cyclotomic_poly(m, x), x) == f), None)
        if m is None: return None
        out[m] = out.get(m, 0) + e
    return dict(sorted(out.items()))

# Table 2 of math/0611532v1 (tex l. 1500-1585): weight sequence -> (printed factorisation, printed period, mark)
# mark: 'A' Arnold hypersurface, 'H' other hypersurface, '' complete intersection, not hypersurface
T2 = {(2,3,7): ({42:1}, 42, 'A'), (2,3,8): ({2:1,10:1,30:1}, 30, 'A'), (2,3,9): ({3:1,12:1,24:1}, 24, 'A'),
      (2,3,10): ({2:1,16:1,18:1}, 72, ''), (2,4,5): ({2:1,6:1,30:1}, 30, 'A'), (2,4,6): ({2:2,22:1}, 22, 'A'),
      (2,4,7): ({2:1,9:1,18:1}, 18, 'A'), (2,4,8): ({2:2,4:1,12:1,14:1}, 84, ''), (2,5,5): ({5:1,20:1}, 20, 'A'),
      (2,5,6): ({2:1,8:1,16:1}, 16, 'A'), (2,5,7): ({11:1,12:1}, 132, ''), (2,6,6): ({2:2,3:1,6:1,10:1,12:1}, 60, ''),
      (3,3,4): ({3:1,24:1}, 24, 'A'), (3,3,5): ({2:1,3:1,6:1,18:1}, 18, 'A'), (3,3,6): ({3:2,15:1}, 15, 'A'),
      (3,3,7): ({2:1,3:1,4:1,10:1,12:1}, 60, ''), (3,4,4): ({2:1,4:1,16:1}, 16, 'A'), (3,4,5): ({13:1}, 13, 'A'),
      (3,4,6): ({2:1,3:1,9:1,10:1}, 90, ''), (3,5,5): ({2:1,5:1,8:1,10:1}, 40, ''),
      (4,4,4): ({2:2,4:2,6:1,12:1}, 12, 'A'), (4,4,5): ({2:1,4:1,8:1,9:1}, 72, ''),
      (2,2,2,3): ({2:2,18:1}, 18, 'H'), (2,2,2,4): ({2:2,14:1}, 14, 'H'), (2,2,2,5): ({2:2,3:1,6:1,12:1}, 12, 'H'),
      (2,2,2,6): ({2:2,8:1,10:1}, 40, ''), (2,2,3,3): ({2:1,3:1,4:1,12:1}, 12, 'H'), (2,2,3,4): ({2:2,5:1,10:1}, 10, 'H'),
      (2,2,3,5): ({2:1,7:1,8:1}, 56, ''), (2,2,4,4): ({2:2,4:1,6:1,8:1}, 24, ''), (2,3,3,3): ({3:2,9:1}, 9, 'H'),
      (2,3,3,4): ({2:1,3:1,6:1,7:1}, 42, ''), (3,3,3,3): ({2:1,3:3,6:2}, None, 'not periodic'),
      (2,2,2,2,2): ({2:4,10:1}, 10, 'H'), (2,2,2,2,3): ({2:3,4:1,8:1}, 8, 'H'),
      (2,2,2,2,4): ({2:2,3:1,6:2}, None, 'not periodic'), (2,2,2,3,3): ({2:2,3:1,5:1,6:1}, 30, ''),
      (2,2,2,2,2,2): ({2:5,4:1,6:1}, 12, '')}

def lna_cartan(rl):
    n = len(rl) + 2
    e = []
    for i in range(1, n + 1):
        m = n
        for j in range(i, n + 1):
            if j - 1 < len(rl) and rl[j - 1] > 0: m = min(m, j + rl[j - 1] - 1)
        e.append(m)
    return np.array([[1 if i <= j <= e[i - 1] else 0 for j in range(1, n + 1)] for i in range(1, n + 1)], dtype=float)

def coxpoly(C):
    Phi = -np.linalg.solve(C.T, C)
    return tuple(int(round(c)) for c in np.poly(Phi))       # high degree first

def exact_phi(C):
    """Phi = -C^{-T} C over Z (C upper unitriangular, so C^T is lower unitriangular: forward substitution)"""
    n = len(C); A = [[int(C[j][i]) for j in range(n)] for i in range(n)]      # A = C^T
    B = [[int(C[i][j]) for j in range(n)] for i in range(n)]
    X = [[0] * n for _ in range(n)]                                           # solve A X = B
    for i in range(n):
        for j in range(n):
            X[i][j] = B[i][j] - sum(A[i][k] * X[k][j] for k in range(i))
    return [[-X[i][j] for j in range(n)] for i in range(n)]

def mul(P, Q):
    n = len(P); Qt = list(zip(*Q))
    return [[sum(a * b for a, b in zip(P[i], Qt[j])) for j in range(n)] for i in range(n)]

def exact_charpoly(M):
    """Faddeev-LeVerrier over Z, high degree first"""
    n = len(M); c = [1]; Mk = [[0] * n for _ in range(n)]; I = [[int(i == j) for j in range(n)] for i in range(n)]
    for k in range(1, n + 1):
        Mk = mul(M, [[Mk[i][j] + c[-1] * I[i][j] for j in range(n)] for i in range(n)])
        tr = sum(Mk[i][i] for i in range(n)); assert tr % k == 0
        c.append(-tr // k)
    return tuple(c)

def lna_order(Phi, cap=400):
    n = len(Phi); I = [[int(i == j) for j in range(n)] for i in range(n)]; M = Phi
    for k in range(1, cap + 1):
        if M == I: return k
        M = mul(M, Phi)
    return None

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    maxn = int(args[0]) if args else 12
    minn = int(args[1]) if len(args) > 1 else 8
    ci_only = '--ci-only' in sys.argv          # targets restricted to complete intersections that are not hypersurfaces
    print('(1a) phi_T against LdlP Table 2')
    target = {}
    for p, (fac, per, mark) in T2.items():
        P = phiT(p); f = cyclo(P); n = P.degree(); rank = sum(a - 1 for a in p) + 3
        flag = 'OK' if f == fac else 'TABLE DIFFERS'
        print(f"  {p}: rank {rank}, deg {n}, computed {f}, printed {fac} period {per} [{mark or 'CI, not hypersurface'}] {flag}")
        if not ci_only or mark == '': target[tuple(int(c) for c in P.all_coeffs())] = (p, mark)
    print('(1b) control: phi_T(<a,b,c]) = Coxeter polynomial of N_n(r) for LMR Thm 6.3 identifications')
    for (n, r, p) in [(11,4,(2,4,5)), (12,3,(2,3,7)), (12,5,(2,5,5)), (13,3,(2,3,8)), (13,5,(2,5,6)), (13,6,(2,4,7)),
                      (14,3,(2,3,9)), (14,4,(2,5,7)), (14,5,(2,6,6)), (15,3,(2,3,10)), (15,9,(2,3,10)), (15,10,(2,3,10)),
                      (16,8,(2,3,11)), (16,10,(2,3,11)), (17,10,(2,3,11))]:   # last row: planted wrong pair (rank 16 != 17)
        rl = [r] * (n - r) + [0] * (r - 2)
        ok = coxpoly(lna_cartan(rl)) == tuple(int(c) for c in phiT(p).all_coeffs())
        print(f"  N_{n}({r}) vs {p}: {'equal' if ok else 'DIFFERENT'}")
    print(f'(2) LNAs n = {minn}..{maxn} (ci_only={ci_only}) whose Coxeter polynomial is phi_T of a Table-2 sequence')
    for n in range(minn, maxn + 1):
        hits = collections.defaultdict(list); fp = 0
        for rl in allRelationLengths(n):
            C = lna_cartan(rl); cp = coxpoly(C)
            if cp in target:                                  # float filter, then exact confirmation
                Phi = exact_phi(C)
                if exact_charpoly(Phi) == cp: hits[target[cp][0]].append((rl, Phi))
                else: fp += 1
        for p, L in hits.items():
            mark = target[tuple(int(c) for c in phiT(p).all_coeffs())][1]
            orders = collections.Counter(lna_order(Phi) for rl, Phi in L)
            show = [''.join(str(d) if d < 10 else f'({d})' for d in rl) for rl, _ in L[:4]]
            print(f"  n={n} {p} [{mark or 'CI, not hypersurface'}]: {len(L)} LNAs, Phi-orders {dict(orders)}, first {show}", flush=True)
        print(f"  n={n}: float-filter false positives {fp}", flush=True)
