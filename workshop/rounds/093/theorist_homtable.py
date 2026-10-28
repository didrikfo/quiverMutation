"""Round 093 theorist, T10: Hom(T,T[-1]), Hom(T,T[1]) (dim totals) at S = -1, 0 along the recorded paths to the key-kept
J != 0 steps E-233 (root 3060000, F1 R9 R8 R8 | R8), E-239 (root 3304000, F2 R8 R9 R7 R7 | R7),
E-254 (root 2400230, F1 F1 F2 F3 R9 R7 R9 F3 | F3). Every prefix step is a J = 0 step of the E-212 move set (tiltingPlus true),
so the table lists J = 0 steps (controls) and the one J != 0 step per path.
Run from repo root:  .venv/bin/python workshop/rounds/093/theorist_homtable.py"""
import sys, time
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import numpy as np
from quivermutation import invariants

def child(alg, kind, v):
    a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
    raw = mutation.quiverMutationAtVertex(a, v)
    c = reduction.reducePathAlgebra(raw)
    return c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)

def row(x, kind, v):
    a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
    tp = tiltingPlus(a.quiver, procedure.relationsFrom(a), v); pi = perI(a, v)
    out = []
    for S in (-1, 0):
        Vs, res = stepTest(a, v, S)
        out.append((int(np.array(res[-1]).sum()), int(np.array(res[1]).sum())))
    return tp, pi, out

CASES = [('E-233', '3060000', 'F1 R9 R8 R8', 'R8'), ('E-239', '3304000', 'F2 R8 R9 R7 R7', 'R7'),
         ('E-254', '2400230', 'F1 F1 F2 F3 R9 R7 R9 F3', 'F3')]
print('case step kind v | tiltingPlus perI | S=-1: (sum Hom[-1], sum Hom[1]) | S=0: (..)')
for name, root, path, last in CASES:
    x = nk.LinearNakayamaAlgebra.fromClassName(root)
    for st in path.split() + [last]:
        kind, v = st[0], int(st[1:])
        t0 = time.time()
        tp, pi, out = row(x, kind, v)
        tag = 'J!=0 (final)' if st == last and (st is path.split()[-1] or True) and False else ''
        print(name, st, '|', tp, pi, '| S=-1', out[0], '| S=0', out[1], '| %.0fs' % (time.time() - t0), flush=True)
        if st != last or True:
            try: x = child(x, kind, v)
            except Exception as e: print('child failed', e); break
    print('---')
