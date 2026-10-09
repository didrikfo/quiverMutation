"""Do Cartan-matrix (Euler form) invariants cut the polynomial-match candidates by (cords, relations)?
For each quipu-with-relations algebra with the polynomial of an LNA outside a quipu class:
signature and Smith form of the symmetrised Euler form C+C^T (derived invariants, Ladkani math/0610685 3.15).
usage: python workshop/rounds/004/maverick_euler.py N"""
import sys, collections
import numpy as np, sympy
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk, invariants as inv

n = int(sys.argv[1])
def cartan_from_match(m):
    C = np.eye(n, dtype=int)
    succ = collections.defaultdict(list)
    for (a, b), bit in zip(m['edges'], m['orientation']):
        (t, h) = (a, b) if bit else (b, a); succ[t].append(h)
    rels = [tuple(p) for p in m['relations']]
    def has_rel(path):
        return any(path[i:i+len(r)] == r for r in rels for i in range(len(path)-len(r)+1))
    def walk(path):
        for h in succ[path[-1]]:
            q = path + (h,)
            if has_rel(q): continue
            C[path[0], h] += 1; walk(q)
    for v in range(n): walk((v,))
    return C
def sym_inv(C):
    G = C + C.T
    ev = np.linalg.eigvalsh(G.astype(float))
    sig = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()), int((abs(ev) <= 1e-9).sum()))
    snf = smith_normal_form(sympy.Matrix(G.tolist()), domain=sympy.ZZ)
    return sig, tuple(sorted(abs(int(snf[i, i])) for i in range(n)))
status = ct.lnaStatus(n)
res = qr.search(n, minArrows=2, statuses=(ct.NOT_QUIPU, ct.UNPLACED), keepPerKey=None)
# LNA invariants
lna = {}
for rl in status:
    if status[rl] != ct.QUIPU:
        A = nk.LinearNakayamaAlgebra(n, list(rl)); C = np.array(inv.cartanMatrix(A), dtype=int)
        lna[rl] = sym_inv(C)
print("LNA sym invariants:", {ct.className(r): v for r, v in lna.items()})
tab = collections.defaultdict(collections.Counter)
for key, group in res['examples'].items():
    targets = {sym_inv(np.array(inv.cartanMatrix(nk.LinearNakayamaAlgebra(n, list(rl))), dtype=int))
               for (name, st) in group[0]['lnas'] for rl in [tuple(map(int, name))] if False} 
    lnas = [nm for nm, st in group[0]['lnas']]
    want = {lna[tuple(int(c) for c in nm)] for nm in lnas}
    for m in group:
        cords = sum(1 for x in m['parameters'][1] if x > 0)
        iv = sym_inv(cartan_from_match(m))
        tab[(key, iv in want)][(cords, len(m['relations']))] += 1
for (key, ok), c in tab.items():
    print("poly", key, "symmetric-form invariant matches LNA:", ok, "n=", sum(c.values()))
    print("   (cords,rels):diag counts", {f"{a},{b}": v for (a, b), v in sorted(c.items())})
