"""Round 137 (toolsmith), T10 decision data: `refuseNonTilting` on vs off, library F-only guarded walks (E-332 roots, n = 11).

Pre-registration: workshop/rounds/137/toolsmith.md (written before any run).

Walk = search.mutationSearchDepthFirst(start, D, printOutput=False, visited=fingerprint.Visited(), refuseNonTilting=FLAG),
defaults otherwise (key guard on, base key = start's).  start = LNA ROOT (`root`) or pathAlgebra.dualPathAlgebra of it (`dual`).
Step number of a gate call = depth of the child (1..D), tracked by wrapping search.mutationSearchDepthFirst (the recursion looks
the name up in the module, so the wrapper sees every level).  Every non-tilting admitted step (off: computed by a spy on the gate;
on: the option's own `isTiltingAtVertex` calls that return False) is recorded with step, vertex, legality, key kept, child class.

modes (repository root; node files go to $T10OUT, default /tmp/toolsmith137):
  walk ROOT SIDE D off|on [--plant FILE]   one walk; JSON summary on stdout; nodes {canonicalKey: (inClass, path)} pickled
  compare ROOT SIDE D [TAG]                 L_in, L_x, G between the off and on (or on-TAG) node files; lost in-class paths
  pick ROOT SIDE D K FILE                   power: write K in-class nodes at step D-1 of the off ball to FILE (pickle of keys)
  shard ROOT SIDE D off|on K M              shard M of K (by step-1 child) of the same walk; then
  merge ROOT SIDE D off|on K                union of the K shard files = the walk's node file (checked against `walk` at D = 5)
  time ROOT SIDE D off|on                   the same walk with no spy or wrapper: node count and seconds only (cost of the option)
  ROOT = e316-I (I = 0..18), SIDE = step:   start = parent of E-316's I-th J != 0 step (dual if an R step), so the step is step 1
  regress                                   R2 (E-318), R3 (E-323), R4 (E-327 seeds 514, 3277)
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py walk 303500030 root 7 off
"""
import sys, os, json, time, copy, pickle, collections, random
sys.path.insert(0, '.')
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import search, mutation, procedure, reduction, pathAlgebra, fingerprint, invariants
from quivermutation import nakayama as nk, arrowPaths as ap

OUT = os.environ.get('T10OUT', '/tmp/toolsmith137')

def inv(A):
    C = invariants.integerCartanMatrix(A); n = len(C)
    S = [[int(C[i][j]) + int(C[j][i]) for j in range(n)] for i in range(n)]
    M = smith_normal_form(sp.Matrix(S), domain=sp.ZZ); d = tuple(sorted(abs(int(M[i, i])) for i in range(n)))
    ev = np.linalg.eigvalsh(np.array(S, dtype=float)); return d, (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))

E316 = [  # E-316's 19 key-kept perI J != 0 steps (workshop/rounds/121/skeptic_jcount.py replay): (record, root, path to parent, step)
    ('E-233', '3060000', 'F1 R9 R8 R8', 'R8'),
    ('E-234 1275', '3060000', 'F1 F4 R9 R8 R8', 'R8'), ('E-234 1433', '3060000', 'F1 F5 R9 R8 R8', 'R8'),
    ('E-234 1746', '3060000', 'F1 R5 R9 R8 R8', 'R8'), ('E-234 1818', '3060000', 'F1 R6 R9 R8 R8', 'R8'),
    ('E-234 1915', '3060000', 'F1 R8 R9 F8 R9', 'R9'), ('E-234 1932', '3060000', 'F1 R8 R9 R9 R8', 'R8'),
    ('E-234 1970', '3060000', 'F1 R9 F1 R8 R8', 'R8'), ('E-234 2003', '3060000', 'F1 R9 F3 R8 R8', 'R8'),
    ('E-234 4551', '3060000', 'R8 R9 R9 F1 R8', 'R8'), ('E-234 4722', '3060000', 'R9 R8 R8 F1 F8', 'R1'),
    ('E-239 2026', '3304000', 'F2 R8 R9 R7 R7', 'R7'),
    ('E-244 4606', '3304000', 'F2 F1 R8 R9 R7 R7', 'R7'), ('E-244 8063', '3304000', 'F2 R8 R7 R9 F7 R9', 'R9'),
    ('E-244 8170', '3304000', 'F2 R8 R9 F2 R7 R7', 'R7'), ('E-244 8196', '3304000', 'F2 R8 R9 F4 R7 R7', 'R7'),
    ('E-244 8427', '3304000', 'F2 R9 R8 R7 R9 R7', 'R7'), ('E-244 11123', '3304000', 'R8 R9 R7 R7 F2 F7', 'R2'),
    ('E-254', '2400230', 'F1 F1 F2 F3 R9 R7 R9 F3', 'F3')]

def e316Start(i):
    """Start of the review-137 item-2 walk: E-316 step i's parent P (rebuilt along its recorded path), or dual(P) when the step is an R
    move, so that the recorded J != 0 step is the library walk's F step at its vertex from the start (step 1 of the walk)."""
    rec, root, path, step = E316[i]
    P = stepsFrom(root, [(t[0], int(t[1:])) for t in path.split()])
    return pathAlgebra.dualPathAlgebra(P) if step[0] == 'R' else P

def startOf(root, side):
    if root.startswith('e316-'): return e316Start(int(root[5:]))
    r = nk.LinearNakayamaAlgebra.fromClassName(root)
    return pathAlgebra.dualPathAlgebra(r) if side == 'dual' else r

def childOf(a, v):
    raw = mutation.quiverMutationAtVertex(copy.deepcopy(a), v)
    if any(ap.isIllegalRelation(raw.quiver, r) for r in procedure.relationsFrom(raw)): return False, None
    return True, reduction.reducePathAlgebra(raw)

_gate, _tilt, _dfs = mutation.mutationIsPossibleAtVertex, mutation.isTiltingAtVertex, search.mutationSearchDepthFirst
LEVEL = [0]
def _levelled(*a, **k):
    LEVEL[0] += 1
    try: return _dfs(*a, **k)
    finally: LEVEL[0] -= 1
search.mutationSearchDepthFirst = _levelled

def fpath(root, side, D, tag): return os.path.join(OUT, '%s_%s_d%d_%s.pkl' % (root, side, D, tag))

def walk(root, side, D, flag, plantFile=None):
    os.makedirs(OUT, exist_ok=True)
    s = startOf(root, side); i0 = inv(s); base = search._coxeterKeyOrNone(s)
    plant = pickle.load(open(plantFile, 'rb')) if plantFile else set()
    rows = []; C = collections.Counter()
    def record(alg, v):
        legal, c = childOf(alg, v)
        kept = legal and search._coxeterKeyOrNone(c) == base
        rows.append(dict(step=LEVEL[0], v=v, legal=legal, kept=kept, childKey=fingerprint.canonicalKey(c) if kept else None,
                         childClass=('in' if inv(c) == i0 else 'cross') if kept else None))
    def gateSpy(alg, v):
        ok = _gate(alg, v)
        if ok:
            C['admitted'] += 1
            if flag == 'off' and not _tilt(alg, v): record(alg, v)
        return ok
    def tiltSpy(alg, v):
        t = _tilt(alg, v)
        if t and plant:
            legal, c = childOf(alg, v)
            if legal and fingerprint.canonicalKey(c) in plant: C['plantedRefusals'] += 1; return False
        if not t: record(alg, v)
        return t
    mutation.mutationIsPossibleAtVertex = gateSpy; mutation.isTiltingAtVertex = tiltSpy
    got = {}; visits = [0]
    def visit(A, path):
        visits[0] += 1; k = fingerprint.canonicalKey(A)
        if k not in got: got[k] = (inv(A) == i0, tuple(path))
    t0 = time.time()
    search.mutationSearchDepthFirst(copy.deepcopy(s), D, printOutput=False, visitor=visit, visited=fingerprint.Visited(),
                                    refuseNonTilting=(flag == 'on'))
    secs = round(time.time() - t0, 1)
    mutation.mutationIsPossibleAtVertex = _gate; mutation.isTiltingAtVertex = _tilt
    tag = flag if not plantFile else 'on-plant'
    pickle.dump(got, open(fpath(root, side, D, tag), 'wb'))
    keptRows = [r for r in rows if r['kept']]
    print(json.dumps(dict(root=root, side=side, D=D, flag=tag, snfTail=i0[0][-3:], nodes=len(got), visits=visits[0],
                          crossNodes=sum(1 for c, _ in got.values() if not c), admitted=C['admitted'],
                          plantedRefusals=C['plantedRefusals'], nonTiltingSteps=len(rows),
                          nonTiltingByStep=dict(sorted(collections.Counter(r['step'] for r in rows).items())),
                          keyKeptNonTilting=[(r['step'], r['v'], r['childClass']) for r in keptRows],
                          keyKeptNonTiltingBelowHorizon=sum(r['step'] < D for r in keptRows),
                          keyKeptDistinctChildrenByStep=dict(sorted(collections.Counter(st for st, _ in
                                                              {(r['step'], r['childKey']) for r in keptRows}).items())),
                          keyKeptDistinctChildren=len({r['childKey'] for r in keptRows}), secs=secs)), flush=True)

def timed(root, side, D, flag):
    """Cost of the option: the library walk with no spy or wrapper (review 137 item 4)."""
    mutation.mutationIsPossibleAtVertex, mutation.isTiltingAtVertex, search.mutationSearchDepthFirst = _gate, _tilt, _dfs
    s = startOf(root, side); keys = set(); t0 = time.time()
    _dfs(copy.deepcopy(s), D, printOutput=False, visitor=lambda A, p: keys.add(fingerprint.canonicalKey(A)),
         visited=fingerprint.Visited(), refuseNonTilting=(flag == 'on'))
    print(json.dumps(dict(root=root, side=side, D=D, flag=flag, unwrapped=True, nodes=len(keys), secs=round(time.time() - t0, 1))), flush=True)

def shard(root, side, D, flag, K, m):
    """Shard m of K of the depth-D walk: the root's step-1 children (gate, legal, key kept, and tilting if `on`) with index = m mod K,
    each walked to depth D - 1 (own Visited per shard, base key = the root's).  The union over m of the node sets is the depth-D ball
    (a ball is the union of its children's balls; Visited only dedups).  Writes the node file with tag FLAG-sK-m."""
    os.makedirs(OUT, exist_ok=True)
    s = startOf(root, side); i0 = inv(s); base = search._coxeterKeyOrNone(s)
    kids = []; rows = []
    for v in sorted(s.quiver.nodes, reverse=True):
        if not _gate(s, v): continue
        t = _tilt(s, v); legal, c = childOf(s, v)
        kept = legal and search._coxeterKeyOrNone(c) == base
        if not t: rows.append(dict(step=1, v=v, legal=legal, kept=kept, childClass=None))
        if kept and (t or flag == 'off'): kids.append((v, c))
    got = {fingerprint.canonicalKey(s): (True, ())} if m == 0 else {}
    C = collections.Counter(); LEVEL[0] = 1
    def gateSpy(alg, v):
        ok = _gate(alg, v)
        if ok and flag == 'off' and not _tilt(alg, v): record(alg, v)
        return ok
    def tiltSpy(alg, v):
        t = _tilt(alg, v)
        if not t: record(alg, v)
        return t
    def record(alg, v):
        legal, c = childOf(alg, v); kept = legal and search._coxeterKeyOrNone(c) == base
        rows.append(dict(step=LEVEL[0], v=v, legal=legal, kept=kept, childKey=fingerprint.canonicalKey(c) if kept else None, childClass=('in' if inv(c) == i0 else 'cross') if kept else None))
    mutation.mutationIsPossibleAtVertex = gateSpy; mutation.isTiltingAtVertex = tiltSpy
    vis = fingerprint.Visited(); t0 = time.time(); mine = [kv for j, kv in enumerate(kids) if j % K == m]
    for v, c in mine:
        def visit(A, path, v=v):
            k = fingerprint.canonicalKey(A)
            if k not in got: got[k] = (inv(A) == i0, (v,) + tuple(path))
        search.mutationSearchDepthFirst(copy.deepcopy(c), D - 1, printOutput=False, visitor=visit, visited=vis, baseKey=base,
                                        refuseNonTilting=(flag == 'on'))
    mutation.mutationIsPossibleAtVertex = _gate; mutation.isTiltingAtVertex = _tilt
    if m != 0: rows = [r for r in rows if r['step'] > 1]
    pickle.dump(got, open(fpath(root, side, D, '%s-s%d-%d' % (flag, K, m)), 'wb'))
    keptRows = [r for r in rows if r['kept']]
    print(json.dumps(dict(root=root, side=side, D=D, flag=flag, shard='%d/%d' % (m, K), children=[v for v, _ in mine], nodes=len(got),
                          nonTiltingByStep=dict(sorted(collections.Counter(r['step'] for r in rows).items())),
                          keyKeptNonTilting=[(r['step'], r['v'], r['childClass']) for r in keptRows],
                          keyKeptNonTiltingBelowHorizon=sum(r['step'] < D for r in keptRows), secs=round(time.time() - t0, 1))), flush=True)

def merge(root, side, D, flag, K):
    got = {}
    for m in range(K):
        for k, val in pickle.load(open(fpath(root, side, D, '%s-s%d-%d' % (flag, K, m)), 'rb')).items():
            if k not in got or len(val[1]) < len(got[k][1]): got[k] = val
    pickle.dump(got, open(fpath(root, side, D, flag), 'wb')); print('merged %s %s D=%d %s: %d nodes' % (root, side, D, flag, len(got)))

def compare(root, side, D, tag='on'):
    off = pickle.load(open(fpath(root, side, D, 'off'), 'rb')); on = pickle.load(open(fpath(root, side, D, tag), 'rb'))
    lost = [k for k in off if k not in on]; gained = [k for k in on if k not in off]
    Lin = [k for k in lost if off[k][0]]
    print(json.dumps(dict(root=root, side=side, D=D, tag=tag, off=len(off), on=len(on), L_in=len(Lin), L_x=len(lost) - len(Lin),
                          G=len(gained), G_in=sum(on[k][0] for k in gained), crossOn=sum(not c for c, _ in on.values()),
                          lostCrossByStep=dict(sorted(collections.Counter(len(off[k][1]) for k in lost if not off[k][0]).items())),
                          lostInPaths=[' '.join(map(str, off[k][1])) for k in Lin][:10])), flush=True)

def pick(root, side, D, K, out):
    off = pickle.load(open(fpath(root, side, D, 'off'), 'rb'))
    cand = sorted(k for k, (c, p) in off.items() if c and len(p) == D - 1)
    keys = set(random.Random(137).sample(cand, K)); pickle.dump(keys, open(out, 'wb'))
    print('picked %d of %d in-class step-%d nodes' % (K, len(cand), D - 1))

def stepsFrom(name, steps):
    """LNA `name`, then F/R steps as the E-318 test builds them (R = mutate the dual, read back)."""
    a = nk.LinearNakayamaAlgebra.fromClassName(name)
    for kind, v in steps:
        side = a if kind == 'F' else pathAlgebra.dualPathAlgebra(a)
        assert _gate(side, v) and _tilt(side, v), (kind, v)
        c = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(side, v))
        a = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
    return a

def regress():
    for tag, name, steps, v in (('R2 E-318', '003033030', (('F', 3), ('R', 8), ('F', 3)), 3),
                                ('R3 E-323', '0040033030', (('F', 3), ('R', 9), ('F', 3)), 3)):
        P = stepsFrom(name, steps); base = search._coxeterKeyOrNone(nk.LinearNakayamaAlgebra.fromClassName(name))
        legal, c = childOf(P, v); ck = fingerprint.canonicalKey(c)
        res = {}
        for flag in (False, True):
            keys = []
            search.mutationSearchDepthFirst(copy.deepcopy(P), 1, printOutput=False, baseKey=base, refuseNonTilting=flag,
                                            visitor=lambda A, p: keys.append(fingerprint.canonicalKey(A)))
            res[flag] = ck in keys
        print('%s: gate %s, key kept %s, tilting %s, SNF tail %s -> %s; depth-1 walk reaches child: off %s, on %s -> %s'
              % (tag, _gate(P, v), search._coxeterKeyOrNone(c) == base, _tilt(P, v), inv(P)[0][-3:], inv(c)[0][-3:], res[False], res[True],
                 'PASS' if (_gate(P, v) and not _tilt(P, v) and res[False] and not res[True]) else 'FAIL'), flush=True)
    sys.path.insert(0, 'workshop/rounds/118'); sys.path.insert(0, 'workshop/rounds/091')
    import skeptic_conv3 as c3, skeptic_convention as sc
    P = sc.P; sg = lambda x: x - P if x > P // 2 else x
    for seed in (514, 3277):
        A = c3.gen('N', 6, random.Random(seed)); q, rels = c3.repo_input(A)
        rels = [{k: sg(x) for k, x in r.items()} for r in rels]
        alg = procedure.toPathAlgebra(q, rels); dual = pathAlgebra.dualPathAlgebra(alg)
        for i in range(A.n):
            for side, B, hs in (('F', alg, 'R'), ('R', dual, 'L')):
                if not any((a if side == 'F' else b) == i for a, b, _ in A.arrows): continue
                if not _gate(B, i + 1): continue
                hom = sc.invariants(A, i, hs)[0]; lib = _tilt(B, i + 1)
                if not hom or not lib:
                    print('R4 seed %d v%d %s: gate True, Hom-tilting %s, library isTiltingAtVertex %s -> %s'
                          % (seed, i + 1, side, hom, lib, 'PASS' if lib == hom else 'FAIL'), flush=True)

if __name__ == '__main__':
    A = sys.argv[1:]
    if A[0] == 'walk': walk(A[1], A[2], int(A[3]), A[4], A[A.index('--plant') + 1] if '--plant' in A else None)
    elif A[0] == 'compare': compare(A[1], A[2], int(A[3]), *(A[4:5]))
    elif A[0] == 'pick': pick(A[1], A[2], int(A[3]), int(A[4]), A[5])
    elif A[0] == 'regress': regress()
    elif A[0] == 'time': timed(A[1], A[2], int(A[3]), A[4])
    elif A[0] == 'shard': shard(A[1], A[2], int(A[3]), A[4], int(A[5]), int(A[6]))
    elif A[0] == 'merge': merge(A[1], A[2], int(A[3]), A[4], int(A[5]))
