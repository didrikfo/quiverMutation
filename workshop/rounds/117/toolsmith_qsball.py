"""Round 117 (toolsmith): backlog 1 / T10. NO-capable Cartan-level instrument on REAL key-kept steps at n = 9.
Roots: the two q_S classes of the F-010 key group (3060000 vs 3304000; also 2400230). BFS over the repo's key-kept J = 0 tilting moves (`moves`) to DEPTH, and
for every ball node the key-kept perI J != 0 children (`jcands` + `child`).  For each parent->child step compare the invariant (sig, |det S|, q_S histogram) of
the integer Cartan matrices.  A step whose invariant differs would be a key-kept step that is NOT a derived equivalence (NO).  Also the planted control: the
cross-class pairs (parent from one ball, child from the other ball) -- the invariant must say NO on them.
usage (repo root): DEPTH=4 timeout 10m .venv/bin/python workshop/rounds/117/toolsmith_qsball.py ROOT [ROOT2 ...]"""
import sys, os, time, collections
ROOTS = sys.argv[1:] or ['3060000', '3304000']; DEPTH = int(os.environ.get('DEPTH', '4'))
sys.path.insert(0, '.'); sys.argv = ['x', 'none', '2', '5', '6', '400000', 'paths']; KEYROOT='3060000'
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName('3060000'))
sys.path.insert(0, 'workshop/rounds/117')
from fractions import Fraction
import numpy as np, sympy as sp
from quivermutation import invariants
import toolsmith_qsscan as Q
_sj = open('workshop/rounds/082/skeptic_join.py').read(); exec(_sj[_sj.index('def child('):_sj.index('def struct(')])
def jcands(alg):
    V = sorted(alg.quiver.nodes); out = []
    for kind in 'FR':
        a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(a, v) or not perI(a, v): continue
            try:
                raw = mutation.quiverMutationAtVertex(a, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): out.append((kind, v, 'illegal')); continue
                c = reduction.reducePathAlgebra(raw)
                if sorted(c.vertices()) != V: out.append((kind, v, 'verts')); continue
                if kind == 'R': c = pathAlgebra.dualPathAlgebra(c)
                k = search._coxeterKeyOrNone(c)
                out.append((kind, v, 'none' if k is None else ('kept' if k == base else 'moved')))
            except Exception as e: out.append((kind, v, 'err ' + type(e).__name__))
    return out

_cache = {}
def qinv(alg):
    C = np.array(invariants.integerCartanMatrix(alg)).astype(int); k = (C.tobytes(), C.shape)
    if k in _cache: return _cache[k]
    S = sp.Matrix((C + C.T).tolist()); n = C.shape[0]
    ev = np.linalg.eigvalsh(np.array(S.tolist(), dtype=float)); sg = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
    d = int(S.det()); Si = S.inv(); cols = [[Fraction(int(Si[i, j].p), int(Si[i, j].q)) for i in range(n)] for j in range(n)]
    kf = lambda v: tuple(x % 1 for x in v); seen = {kf([Fraction(0)] * n)}; fr = list(seen)
    while fr:
        nf = []
        for g in fr:
            for c in cols:
                h = kf([a + b for a, b in zip(g, c)])
                if h not in seen: seen.add(h); nf.append(h)
        fr = nf
    Sl = [[int(S[i, j]) for j in range(n)] for i in range(n)]
    hh = collections.Counter(sum(y[i] * Sl[i][j] * y[j] for i in range(n) for j in range(n)) % 2 for y in seen)
    _cache[k] = (sg, abs(d), tuple(sorted(hh.items()))); return _cache[k]
def ballrun(root):
    t0 = time.time(); r = nk.LinearNakayamaAlgebra.fromClassName(root)
    seen = {fingerprint.canonicalKey(r): r}; front = [r]; steps = []; tags = []
    for d in range(DEPTH):
        nxt = []
        for x in front:
            qx = qinv(x)
            for kind, v, c in moves(x):
                steps.append(('J0', qx, qinv(c)))
                k = fingerprint.canonicalKey(c)
                if k not in seen: seen[k] = c; nxt.append(c)
            if d < DEPTH:
                for kind, v, tag_ in jcands(x):
                    if tag_ != 'kept': continue
                    try: c = child(x, kind, v)
                    except Exception: continue
                    steps.append(('Jne0', qx, qinv(c)))
        front = nxt; print(root, 'depth', d + 1, 'nodes', len(seen), 'steps so far', len(steps), '%.0fs' % (time.time() - t0), flush=True)
    return seen, steps
res = {}
for root in ROOTS:
    seen, steps = ballrun(root); res[root] = (seen, steps)
    cnt = collections.Counter((s[0], s[1] == s[2]) for s in steps)
    print('ROOT', root, 'nodes', len(seen), 'root invariant hist', qinv(seen[fingerprint.canonicalKey(nk.LinearNakayamaAlgebra.fromClassName(root))])[2],
          'steps (kind, invariant equal):', dict(cnt), flush=True)
    print('  distinct node invariants in ball:', len({qinv(a) for a in seen.values()}), flush=True)
# planted control: cross-ball pairs
inv_of = {r: {qinv(a) for a in res[r][0].values()} for r in res}
for a in ROOTS:
    for b in ROOTS:
        if a < b: print('planted cross pair', a, b, 'invariant sets disjoint:', inv_of[a].isdisjoint(inv_of[b]), '(NO = differ on all pairs)')
