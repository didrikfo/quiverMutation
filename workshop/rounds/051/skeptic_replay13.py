"""Round 050 (skeptic), T10 (i): replay the printed E-157 paths (rounds/049/toolsmith_paths_logs.txt) and test EVERY edge with the independent
tilting test of skeptic_tilt.py (no tiltingPlus / perI / gate used for the verdict; they are printed for context).
Usage (repo root): skeptic_replay.py in.pkl cls   where in.pkl = rounds/049/toolsmith_collect.py output for class cls (7 1 / 7 2).
Per edge x -> y (move F v: step at v of x; move R v: step at v of dual(x), carried back): a = x or dual(x), c = reduced mutation of a at v.
Verdict 'tilt' := Hom(T,T[-1]) = 0 and Hom(T,T[1]) = 0 for T = T_v(a) of skeptic_tilt.py; 'cartan' := H = Cartan(c)^T (labelled, c = child before dualizing).
Also checks the paths: child-side end key == LNA-side end key."""
import sys, re, pickle
import numpy as np
A_ = sys.argv; PK, CLS = '/tmp/tsm/c1.pkl', 1
sys.argv = ['x', PK, str(CLS), '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import hashlib
from quivermutation import invariants
def tag(alg): return hashlib.md5(repr(fingerprint.canonicalKey(alg)).encode()).hexdigest()[:6]
def step(x, kind, v):
    """returns (y, a, c_a, V) or None"""
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
    zero = not np.array(res[-1]).any() and not np.array(res[1]).any()
    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
    H = np.array(res[0], dtype=int)
    cart = bool((H == C.T).all())
    ctx = dict(J=bool(perI(a, v)), tp=bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v)))
    return y, zero, cart, ctx, int(np.abs(np.array(res[-1])).sum())

recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
CS='F1 F3 F1 F5 F4 R7 R1'.split(); TS='F7 R2 R1 F2 F7 R3'.split()
x0 = recs[14]['childObj']; print('child key', tag(x0))
ends=[]
for x, mv in [(x0, CS), (classes[base][9], TS)]:
    for s_ in mv:
        kind, v = s_[0], int(s_[1:]); vd = verdict(x, kind, v)
        if vd is None: print(s_, 'NOSTEP'); break
        y, zero, cart, ctx, hm1 = vd
        print(s_, 'tilt', zero, 'cartan', cart, 'hm1', hm1, ctx, flush=True); x = y
    ends.append(fingerprint.canonicalKey(x))
print('meet', ends[0]==ends[1], ends[0] is not None)
