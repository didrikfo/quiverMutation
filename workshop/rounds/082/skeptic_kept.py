"""Round 082 (skeptic), T10: find the key-kept perI-flagged (J != 0) step of E-229 in the 3060000 ball and Hom-test it.
usage: skeptic_kept.py find ROOT DEPTH OUT.pkl   (BFS as experimentalist_jkey.py with parent pointers; saves kept steps + J=0 path from root)
       skeptic_kept.py test IN.pkl               (Hom(T,T[m]) m=-1,0,1, End vs Cartan, S=0 and S=-1, at the kept step)
Run from the repository root.  Machinery (moves, perI, stepTest) is the repo's / round 049 / round 050 code, loaded as in experimentalist_jkey.py."""
import sys, time, pickle, collections
A_ = sys.argv; MODE = A_[1]
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
import numpy as np
from quivermutation import invariants
base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000'))

def child(alg, kind, v):
    a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
    raw = mutation.quiverMutationAtVertex(a, v)
    c = reduction.reducePathAlgebra(raw)
    return c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)

def jkept(alg):
    out = []; V = sorted(alg.quiver.nodes)
    for kind in 'FR':
        a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(a, v) or not perI(a, v): continue
            try:
                raw = mutation.quiverMutationAtVertex(a, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                c = reduction.reducePathAlgebra(raw)
                if sorted(c.vertices()) != V: continue
                if kind == 'R': c = pathAlgebra.dualPathAlgebra(c)
                if search._coxeterKeyOrNone(c) == base: out.append((kind, v))
            except Exception: pass
    return out

def rels(a): return [list(r) for r in procedure.relationsFrom(a)]

if MODE == 'find':
    ROOT, DEPTH, OUT = A_[2], int(A_[3]), A_[4]
    root = nk.LinearNakayamaAlgebra.fromClassName(ROOT)
    seen = {fingerprint.canonicalKey(root): (root, None, None, 0)}; front = [root]; found = []
    T0 = time.time()
    for d in range(DEPTH):
        nxt = []
        for x in front:
            kx = fingerprint.canonicalKey(x)
            for kind, v in jkept(x):
                print('KEPT parent depth', d, kind, v, 'parent rels', rels(x), flush=True)
                path = []; k = kx
                while seen[k][1] is not None: path.append(seen[k][2]); k = seen[k][1]
                found.append(dict(parent=x, kind=kind, v=v, depth=d, path=path[::-1], rootalg=root))
            for kind, v, c in moves(x):
                k = fingerprint.canonicalKey(c)
                if k is not None and k not in seen: seen[k] = (c, kx, (kind, v), d + 1); nxt.append(c)
        print('level', d, 'ball', len(seen), '%.0fs' % (time.time() - T0), flush=True)
        front = nxt
    pickle.dump(found, open(OUT, 'wb'))
elif MODE == 'test':
    found = pickle.load(open(A_[2], 'rb'))
    for f in found:
        x = f['parent']; kind, v = f['kind'], f['v']
        print('parent depth', f['depth'], 'path from root', f['path'])
        print('parent relations', rels(x), 'quiver arrows', sorted(x.quiver.edges()))
        a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
        print('step', kind, v, 'perI', perI(a, v), 'tiltingPlus', tiltingPlus(a.quiver, procedure.relationsFrom(a), v))
        c = child(x, kind, v); print('child relations', rels(c))
        ac = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
        C = np.array(invariants.cartanMatrix(ac, exact=True).tolist(), dtype=int)
        for S in (-1, 0):
            Vs, res = stepTest(a, v, S)
            h = {m: np.array(res[m], dtype=int) for m in (-1, 0, 1)}
            print('S', S, 'Hom(T,T[-1]) total', h[-1].sum(), 'entries', [(Vs[i], Vs[j], int(h[-1][i][j])) for i in range(len(Vs)) for j in range(len(Vs)) if h[-1][i][j]],
                  '| Hom(T,T[1]) total', h[1].sum(), [(Vs[i], Vs[j], int(h[1][i][j])) for i in range(len(Vs)) for j in range(len(Vs)) if h[1][i][j]])
            print('   End(T) dims == Cartan(child)^T', bool((h[0] == C.T).all()), ' == Cartan(child)', bool((h[0] == C).all()))
            if not (h[0] == C.T).all(): print('   H0\n', h[0], '\n   C\n', C)
        print('Coxeter key parent == child', search._coxeterKeyOrNone(x) == search._coxeterKeyOrNone(c))
        print('Cartan det parent/child', round(np.linalg.det(np.array(invariants.cartanMatrix(a, exact=True).tolist(), dtype=float))), round(np.linalg.det(C.astype(float))))
# Meets (after `skeptic_join.py side ...` for C = pkl:KEPT:child, R = 3060000, P = pkl:KEPT:parent, Q = 3304000, depth 5, ckpt files ckC/ckR/ckP/ckQ):
#   skeptic_join.py meet ckC ckR   -> total 7, child path R1 R3 || root path F1 F7 F6 F5 F4, Hom-verified, iso meet
#   skeptic_join.py meet ckC ckP   -> total 6 (parent ball completed to depth 4), F8 F8 || F8 R7 R6 R5
#   skeptic_join.py meet ckC ckQ   -> no meet (0 WL-bucket coincidences)
