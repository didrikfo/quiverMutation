"""Apply the F-047 profile (Smith form of g(Phi) per irreducible factor g of the Coxeter polynomial)
plus the Smith form of C+C^T to every polynomial-match candidate at order N; tabulate survivors by (cords, rels).
usage: python workshop/rounds/004/maverick_profile.py N"""
import sys, collections
import numpy as np, sympy
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, invariants as inv
sys.path.insert(0, "workshop/rounds/004")
import importlib.util
spec = importlib.util.spec_from_file_location("me", "workshop/rounds/004/maverick_euler.py")
n = int(sys.argv[1])
T = sympy.symbols('T')
def cart_match(m):
    C = np.eye(n, dtype=int); succ = collections.defaultdict(list)
    for (a, b), bit in zip(m['edges'], m['orientation']):
        (t, h) = (a, b) if bit else (b, a); succ[t].append(h)
    rels = [tuple(p) for p in m['relations']]
    def bad(p): return any(p[i:i+len(r)] == r for r in rels for i in range(len(p)-len(r)+1))
    def walk(p):
        for h in succ[p[-1]]:
            q = p + (h,)
            if bad(q): continue
            C[p[0], h] += 1; walk(q)
    for v in range(n): walk((v,))
    return C
cache = {}
def profile(C):
    k = C.tobytes()
    if k in cache: return cache[k]
    M = sympy.Matrix(C.tolist())
    Phi = -(M.inv().T) * M
    cp = sympy.Poly(Phi.charpoly(T).as_expr(), T)
    out = [tuple(sorted(abs(int(s[i, i])) for i in range(n))) for s in
           [smith_normal_form(sympy.Matrix(C.tolist()) + sympy.Matrix(C.T.tolist()), domain=sympy.ZZ)]]
    for g, e in sympy.factor_list(cp.as_expr())[1]:
        gm = sympy.Poly(g, T).as_expr()
        G = sympy.zeros(n, n)
        for (d,), c in sympy.Poly(g, T).terms():
            G += c * Phi**d
        s = smith_normal_form(G, domain=sympy.ZZ)
        out.append((str(g), e, tuple(sorted(abs(int(s[i, i])) for i in range(n)))))
    cache[k] = tuple(out); return cache[k]
status = ct.lnaStatus(n)
res = qr.search(n, minArrows=2, statuses=(ct.NOT_QUIPU, ct.UNPLACED), keepPerKey=None)
lp = {}
for rl in status:
    if status[rl] != ct.QUIPU:
        lp[ct.className(rl)] = profile(np.array(inv.cartanMatrix(nk.LinearNakayamaAlgebra(n, list(rl))), dtype=int))
print("LNA profiles distinct:", len(set(lp.values())), {k: hash(v) % 1000 for k, v in lp.items()})
tab = collections.defaultdict(collections.Counter)
for key, group in res['examples'].items():
    want = {lp[nm] for nm, st in group[0]['lnas']}
    for m in group:
        cords = sum(1 for x in m['parameters'][1] if x > 0)
        tab[(key, profile(cart_match(m)) in want)][(cords, len(m['relations']))] += 1
for (key, ok), c in sorted(tab.items(), key=str):
    print("poly", key, "profile matches:", ok, "n=", sum(c.values()))
    print("   ", {f"{a},{b}": v for (a, b), v in sorted(c.items())})
