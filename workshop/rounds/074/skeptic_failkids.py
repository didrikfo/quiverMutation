"""Round 074 (skeptic), T10: collect the key-kept tiltingPlus-false edges of the complete n = 7 balls (E-204).
Same walk as workshop/rounds/070/toolsmith_cyccount.py (copied logic: gate, tiltingPlus, illegal drop, reduce, key guard, canonicalKey dedupe).
  skeptic_failkids.py n CLS DEPTH BUDGET_SEC I/K OUT.pkl   (repo root). DEPTH 9 = parents of depth <= 8. Last level expands frontier[I::K] only.
Saves list of (depth, parentKey, v, parentAlg, childAlg) for every key-kept failing edge, and the tally."""
import sys, time, pickle
from collections import Counter
sys.path.insert(0, '.')
_argv = sys.argv; sys.argv = ['x']
exec(compile(open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n"), 'h015', 'exec'))
sys.argv = _argv
n, cls, depth, budget = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
SH = tuple(int(x) for x in sys.argv[5].split('/')); OUT = sys.argv[6]
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
seen = set(); cur = []
for s in classes[base]:
    k = fingerprint.canonicalKey(s)
    if k is not None:
        if k in seen: continue
        seen.add(k)
    cur.append(s)
import os
T = Counter(); fails = []; t0 = time.time(); stop = False; d0 = 0
CK = OUT + '.frontier'   # checkpoint: frontier before the last level (+ failures of earlier levels); written once, read by every shard
if os.path.exists(CK):
    T, fails, cur, seen = pickle.load(open(CK, 'rb')); d0 = depth - 1
    if SH[0] != 0: T = Counter(); fails = []
for d in range(d0, depth):
    if d == depth - 1 and d0 == 0:
        pickle.dump((T, fails, cur, seen), open(CK + '.tmp', 'wb')); os.replace(CK + '.tmp', CK)
    nxt = []
    for alg in (cur[SH[0]::SH[1]] if d == depth - 1 else cur):
        if time.time() - t0 > budget: stop = True; break
        T['expanded'] += 1
        rels = procedure.relationsFrom(alg)
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            t = tiltingPlus(alg.quiver, rels, v)
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, r) for r in procedure.relationsFrom(raw)): continue
            child = reduction.reducePathAlgebra(raw)
            guard = search._coxeterKeyOrNone(child) == base
            T['edges'] += 1; T['guard'] += guard
            if guard and not t:
                T['fail'] += 1
                if d < depth - 1 or True: fails.append((d, fingerprint.canonicalKey(alg), v, alg, child))
            if guard:
                k = fingerprint.canonicalKey(child)
                if k is not None:
                    if k in seen: continue
                    seen.add(k)
                nxt.append(child)
        if stop: break
    print('depth', d + 1, 'expanded', T['expanded'], 'next', len(nxt), 'fails', T['fail'], '%.0fs' % (time.time() - t0), flush=True)
    cur = nxt
    if stop or not cur: break
print('n', n, 'class', cls, 'depth reached', d + 1, 'budget-stopped', stop, dict(T))
pickle.dump((dict(T), fails, stop), open(OUT, 'wb'))
