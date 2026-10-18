"""Round 073 (toolsmith): ground-truth NO controls at n = 7. Same key-guarded walk as rounds/070/toolsmith_cyccount.py (complete balls, sharded), but every
KEY-MOVING legal child (child Coxeter key != class key, hence certainly NOT derived equivalent to its parent; ground truth = Coxeter polynomial, no tilting
code) is also tested with the independent Hom test (rounds/050/skeptic_tilt.py) and every key-kept FAIL child likewise.
  toolsmith_keymoving.py n CLS DEPTH BUDGET_SEC [I/K]   (repo root).  Prints: key-moving children tested / rejected by Hom(T,T[+-1]) != 0 or Cartan mismatch,
  and for key-kept steps: tiltingPlus true/false vs Hom verdict (agreement)."""

import sys, time
import numpy as np
from collections import Counter
sys.path.insert(0, '.')
_argv = sys.argv; sys.argv = ['x']
exec(compile(open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n"), 'h015', 'exec'))
sys.argv = _argv
exec(open('workshop/rounds/050/skeptic_tilt.py').read())
from quivermutation import invariants
H = Counter()
n, cls, depth, budget = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
SH = tuple(int(x) for x in sys.argv[5].split('/')) if len(sys.argv) > 5 else (0, 1)
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
def cyclicVertices(Q):
    bad = set()
    for comp in nx.strongly_connected_components(Q):
        if len(comp) > 1: bad |= comp
    bad |= {u for u, w in Q.edges() if u == w}
    return bad
seen = set(); cur = []
for s in classes[base]:
    k = fingerprint.canonicalKey(s)
    if k is not None:
        if k in seen: continue
        seen.add(k)
    cur.append(s)
T = Counter(); cycNodes = set(); failNodes = []; t0 = time.time(); stop = False
for d in range(depth):
    nxt = []
    for alg in (cur[SH[0]::SH[1]] if d == depth - 1 else cur):
        if time.time() - t0 > budget: stop = True; break
        T['expanded'] += 1
        if cyclicVertices(alg.quiver): T['cyclic-parent-expanded'] += 1; continue
        rels = procedure.relationsFrom(alg)
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            t = tiltingPlus(alg.quiver, rels, v)
            raw = mutation.quiverMutationAtVertex(alg, v)
            rawCyc = bool(cyclicVertices(raw.quiver))
            if any(ap.isIllegalRelation(raw.quiver, r) for r in procedure.relationsFrom(raw)):
                T['illegal'] += 1; T['illegal-cyclic-raw'] += rawCyc; continue
            child = reduction.reducePathAlgebra(raw)
            cv = cyclicVertices(child.quiver); cyc = bool(cv)
            ck = search._coxeterKeyOrNone(child); guard = ck == base
            T['edges'] += 1; T['edges-raw-cyclic'] += rawCyc
            T['guard' if guard else 'noguard'] += 1
            T['x-' + ('guard' if guard else 'noguard') + '-' + ('tilt' if t else 'FAIL')] += 1   # tilt x guard cross-tally (claim 3)
            T['cyclic-' + ('guard' if guard else 'noguard') + ('-tilt' if t else '-FAIL')] += cyc
            if cyc and guard:
                meets = [w for w in child.vertices() if mutation.mutationIsPossibleAtVertex(child, w)
                         and (cv & (nx.descendants(child.quiver, w) | {w}))]
                T['cyclic-guard-gate-vertices'] += sum(1 for w in child.vertices() if mutation.mutationIsPossibleAtVertex(child, w))
                T['cyclic-guard-cycle-meets-D(w)'] += len(meets)
                key = fingerprint.canonicalKey(child)
                cycNodes.add(key if key is not None else id(child))
            if not t and guard: failNodes.append((d, v, cyc, rawCyc))
            if not cyc and not rawCyc and sorted(child.vertices()) == sorted(alg.vertices()):
                Vs, res = stepTest(alg, v, -1)
                zero = not np.array(res[-1]).any() and not np.array(res[1]).any()
                Cc = np.array(invariants.cartanMatrix(child, exact=True).tolist(), dtype=int)
                acc_ = zero and bool((np.array(res[0], dtype=int) == Cc.T).all())
                H[('keykept' if guard else 'keymoving', 'tiltingPlus' if t else 'tpFALSE', 'HomACC' if acc_ else 'HomREJ')] += 1
            elif guard is False: H[('keymoving', 'not-vertex-preserving-or-cyclic')] += 1
            if guard:
                k = fingerprint.canonicalKey(child)
                if k is not None:
                    if k in seen: continue
                    seen.add(k)
                nxt.append(child)
        if stop: break
    print('depth', d + 1, 'expanded', T['expanded'], 'next', len(nxt), '%.0fs' % (time.time() - t0), flush=True)
    cur = nxt
    if stop or not cur: break
print('n', n, 'class', cls, 'key', base, 'depth reached', d + 1, 'budget-stopped', stop, 'distinct nodes', len(seen))
for k, v in sorted(T.items()): print(k, v)
print('distinct cyclic key-kept nodes', len(cycNodes))
print('failing key-kept children', len(failNodes), 'of which cyclic after reduce', sum(1 for f in failNodes if f[2]), 'cyclic before reduce', sum(1 for f in failNodes if f[3]))

print('HOM TALLY'); [print(k, v) for k, v in sorted(H.items())]
