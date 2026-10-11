"""Round 053 (skeptic), T10 agenda 1: replay the 3 printed E-160 paths (rounds/050/toolsmith.md, 'Printed paths') edge by edge with the
Hom_K(T,T[m]) test of rounds/050/skeptic_tilt.py (same machinery as skeptic_replay13.py / E-161 / E-163); the gate, perI, tiltingPlus are printed for context only.
Usage (repo root): skeptic_replay3.py in.pkl cls childIdx ldx 'child moves' 'LNA/dual moves'
  in.pkl = rounds/050/toolsmith_collect.py 7 cls 20000 in.pkl 100.   childIdx indexes the 'fail' records, ldx indexes classes[base] (LNAs+duals).
Per edge: tilt := Hom(T,T[-1]) = 0 and Hom(T,T[1]) = 0; cartan := H = Cartan(c)^T (labelled).  Also: end keys of the two sides equal and not None;
every intermediate node's canonicalKey tag (None shown), since a None key is where the toolsmith's matcher is blind."""
import sys, re, pickle
import numpy as np
A_ = sys.argv; PK, CLS, IDX, LDX, CS, TS = A_[1], int(A_[2]), int(A_[3]), int(A_[4]), A_[5].split(), A_[6].split()
sys.argv = ['x', PK, str(CLS), '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import hashlib
from quivermutation import invariants
def tag(alg):
    k = fingerprint.canonicalKey(alg)
    return 'None' if k is None else hashlib.md5(repr(k).encode()).hexdigest()[:6]
def step(x, kind, v):
    a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
    V = sorted(x.quiver.nodes)
    raw = mutation.quiverMutationAtVertex(a, v)
    c = reduction.reducePathAlgebra(raw)
    if sorted(c.vertices()) != V: return None
    y = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
    return y, a, c, V
def verdict(x, kind, v):
    r = step(x, kind, v)
    if r is None: return None
    y, a, c, V = r
    Vs, res = stepTest(a, v, -1)
    zm1 = int(np.abs(np.array(res[-1])).sum()); zp1 = int(np.abs(np.array(res[1])).sum())
    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
    H = np.array(res[0], dtype=int)
    cart = bool((H == C.T).all())
    ctx = dict(J=bool(perI(a, v)), tp=bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v)))
    return y, (zm1 == 0 and zp1 == 0), cart, ctx, zm1, zp1
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
x0 = recs[IDX]['childObj']; print('child', IDX, 'key', tag(x0), 'LNA/dual #%d' % LDX, tag(classes[base][LDX]))
ends = []; allok = True; ne = 0
for side, x, mv in [('child', x0, CS), ('lna', classes[base][LDX], TS)]:
    for s_ in mv:
        kind, v = s_[0], int(s_[1:]); vd = verdict(x, kind, v)
        if vd is None: print(side, s_, 'NOSTEP'); allok = False; break
        y, tilt, cart, ctx, a1, a2 = vd; ne += 1
        allok &= tilt and cart
        print(side, s_, 'tilt', tilt, 'cartan', cart, 'hm1', a1, 'hp1', a2, ctx, 'to', tag(y), flush=True); x = y
    ends.append(fingerprint.canonicalKey(x))
meet = ends[0] == ends[1] and ends[0] is not None
print('RESULT edges', ne, 'all-accept', allok, 'meet', meet)
