"""Round 111 toolsmith, T10: discriminant quadratic form q_S (histogram of x^T S x mod 2 over x in S^-1 Z^n / Z^n, exact fractions) and signature of
S = C + C^T for parent and child of every key-kept J != 0 step of n = 7 key class 1, from the checkpoint of workshop/rounds/069/skeptic_c1_13.py
(/tmp/skeptic_c1_13_ck1.pkl).  q_S and signature differ => no P in GL_n(Z) with P C_A P^T = C_B => child not derived equivalent to parent (E-272 necessity).
Control NO: parent vs. child of a DIFFERENT-det planted matrix is not available here, so the control is the E-292 discriminant (a pair that differs) -- see submission.
usage (repo root): .venv/bin/python workshop/rounds/111/toolsmith_c1_13_qs.py"""
import sys, pickle, collections
from fractions import Fraction
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
import numpy as np, sympy as sp
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
def qhist(S):
    n = S.shape[0]; M = sp.Matrix(S.tolist()); Si = M.inv()
    Sinv = [[Fraction(int(Si[i, j].p), int(Si[i, j].q)) for j in range(n)] for i in range(n)]
    cols = [tuple(Sinv[i][j] % 1 for i in range(n)) for j in range(n)]
    add = lambda a, b: tuple((x + y) % 1 for x, y in zip(a, b))
    elems = {tuple([Fraction(0)] * n)}
    for c in cols:
        if c in elems: continue
        cur = c; new = set(elems)
        while cur not in elems:
            new |= {add(e, cur) for e in elems}; cur = add(cur, c)
        elems = new
    Sf = [[Fraction(int(v)) for v in r] for r in S]
    h = collections.Counter()
    for x in elems:
        # x = S^-1 v mod Z^n ; q(x) = x^T S x mod 2 (well defined mod 2 for even S)
        Sx = [sum(Sf[i][j] * x[j] for j in range(n)) for i in range(n)]
        q = sum(x[i] * Sx[i] for i in range(n)) % 2
        h[q] += 1
    return tuple(sorted(h.items()))
def sig(S):
    w = np.linalg.eigvalsh(S.astype(float)); return (int((w > 1e-9).sum()), int((w < -1e-9).sum()), int((abs(w) <= 1e-9).sum()))
if __name__ == '__main__':
    seen, queue, byCart, byKey, nexp, snaps = pickle.load(open('/tmp/skeptic_c1_13_ck1.pkl', 'rb'))
    rows = sorted(byCart.values(), key=lambda r: r['nexp'])
    cache = {}
    def info(C):
        k = C.tobytes()
        if k not in cache:
            S = C + C.T; cache[k] = (sig(S), int(sp.Matrix(S.tolist()).det()), qhist(S) if abs(int(sp.Matrix(S.tolist()).det())) > 0 else None)
        return cache[k]
    tab = collections.Counter(); out = []
    for i, r in enumerate(rows):
        a, b = info(r['CA']), info(r['CB']); same = (a == b)
        out.append((i, r['nexp'], r['kept'], r['tp'], same, a[0], a[1], a[2], b[2]))
        tab[(r['kept'], r['tp'], same)] += 1
    print('steps', len(rows), 'table (kept, tiltingPlus, same sig/det/qS):', dict(tab))
    bad = [o for o in out if o[2] and not o[3] and not o[4]]
    print('key-kept, tiltingPlus-false steps with DIFFERENT sig/det/q_S:', len(bad)); [print(o) for o in bad[:20]]
    first = [o for o in out if o[2] and not o[3]][:13]
    print('distinct (sig, det, q_S) over key-kept failing parents:', {(o[5], o[6], o[7]) for o in out if o[2] and not o[3]})
