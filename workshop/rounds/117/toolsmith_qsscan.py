"""Round 117 (toolsmith): backlog 1 / T10. For every LNA at n = 8, 9 compute (signature of S=C+C^T, det S, q_S histogram of x^T S x mod 2 on S^-1 Z^n / Z^n, sign-exact Fractions);
group by Coxeter key; report key groups whose members differ in that Cartan-level invariant (=> derived-inequivalent, same key).  Also the F-010 triple.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/117/toolsmith_qsscan.py [n ...]"""
import sys, collections
from fractions import Fraction
import numpy as np, sympy as sp
from quivermutation import coxeterTables as ct
def cart(N, rl): return sp.Matrix(np.array(ct.lnaCartanMatrix(N, rl)).astype(int).tolist())
def inv(N, rl):
    C = cart(N, rl); S = C + C.T; n = N
    ev = np.linalg.eigvalsh(np.array(S.tolist(), dtype=float)); sig = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()), int((abs(ev) <= 1e-9).sum()))
    d = int(S.det())
    if d == 0: return sig, 0, None
    Si = S.inv(); cols = [[Fraction(int(Si[i, j].p), int(Si[i, j].q)) for i in range(n)] for j in range(n)]
    key = lambda v: tuple(x % 1 for x in v)
    seen = {key([Fraction(0)] * n)}; fr = list(seen)
    while fr:
        nf = []
        for g in fr:
            for c in cols:
                h = key([a + b for a, b in zip(g, c)])
                if h not in seen: seen.add(h); nf.append(h)
        fr = nf
    Sl = [[int(S[i, j]) for j in range(n)] for i in range(n)]
    h = collections.Counter()
    for y in seen: h[sum(y[i] * Sl[i][j] * y[j] for i in range(n) for j in range(n)) % 2] += 1
    return sig, abs(d), tuple(sorted(h.items()))
if __name__ == '__main__':
    for N in [int(a) for a in sys.argv[1:]] or [8, 9]:
        idx = ct.lnaKeyIndex(N); st = ct.lnaStatus(N); tot = sep = 0
        for key, mem in idx.items():
            if len(mem) < 2: continue
            tot += 1
            vals = collections.defaultdict(list)
            for rl in mem: vals[inv(N, rl)].append("".join(map(str, rl)))
            print(f"n={N} key group of {len(mem)} LNAs: {len(vals)} Cartan-level class(es)", "SEPARATED" if len(vals) > 1 else "", flush=True)
            if len(vals) > 1:
                sep += 1
                for v, m in vals.items(): print("    sig,det,qhist", v[0], v[1], v[2], "members", len(m), m[:3])
        print(f"n={N}: {tot} multi-LNA key groups, {sep} separated by (sig, det, q_S hist)")
