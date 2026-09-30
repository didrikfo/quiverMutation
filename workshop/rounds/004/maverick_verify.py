"""Sharp test of H-017: do below-diagonal (relations < cords) candidates with the class's Coxeter polynomial
AND Smith profile actually reach an LNA of that class by mutation search?  Whole cells (rels-cords = -2,-1),
up to K candidates per (poly, cords, rels) cell, search depth D.  usage: ... N DEPTH K [maxdiag]"""
import sys, collections, time
sys.argv_ = sys.argv
import families as fm
import numpy as np, sympy
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import nakayama as nk, invariants as inv
def cart(m):
    C = np.eye(n, dtype=int); succ = collections.defaultdict(list)
    for (a, b), bit in zip(m['edges'], m['orientation']):
        (t, h) = (a, b) if bit else (b, a); succ[t].append(h)
    rels = [tuple(p) for p in m['relations']]
    def bad(p): return any(p[i:i+len(r)] == r for r in rels for i in range(len(p)-len(r)+1))
    def walk(p):
        for h in succ[p[-1]]:
            q = p + (h,)
            if not bad(q): C[p[0], h] += 1; walk(q)
    for v in range(n): walk((v,))
    return C
def snf(C):
    S = smith_normal_form(sympy.Matrix((C + C.T).tolist()), domain=sympy.ZZ)
    return tuple(sorted(abs(int(S[i, i])) for i in range(n)))
from quivermutation import quipuRelations as qr, coxeterTables as ct
n, depth, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
maxdiag = int(sys.argv[4]) if len(sys.argv) > 4 else -1
res = qr.search(n, minArrows=2, statuses=(ct.NOT_QUIPU, ct.UNPLACED), keepPerKey=None)
seen = collections.Counter()
lnasnf = {}

for key, group in res['examples'].items():
    for m in group:
        cords = sum(1 for x in m['parameters'][1] if x > 0); r = len(m['relations'])
        if r - cords > maxdiag: continue
        for nm, _st in m['lnas']:
            if nm not in lnasnf: lnasnf[nm] = snf(np.array(inv.cartanMatrix(nk.LinearNakayamaAlgebra(n, [int(c) for c in nm])), dtype=int))
        if snf(cart(m)) not in {lnasnf[nm] for nm, _ in m['lnas']}: continue   # Smith-form filter (same survivors as the full F-047 profile)
        cell = (key[:3], cords, r)
        if seen[cell] >= K: continue
        seen[cell] += 1
        t = time.time(); reached = fm.verify(m, n, depth)
        print("cords", cords, "rels", r, m['quipu'], "arrows", fm.arrowsOf(m), "rels", fm.relationsAsVertices(m),
              "-> reached", sorted(reached), "%.0fs" % (time.time() - t), flush=True)
