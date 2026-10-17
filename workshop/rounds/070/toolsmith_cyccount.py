"""Round 070 (toolsmith), T10/T5 item 1: count cyclic children of the key-guarded E-192 balls (hence the cyclic parents the walk could meet).
Same walk as workshop/rounds/047/experimentalist_deepreplay.py (gate, tiltingPlus, illegal-relation drop, reduce, key guard, canonicalKey dedupe),
except that every child is tested for a loop / oriented cycle (networkx is_directed_acyclic_graph), before and after reduction.
A cyclic node is NOT expanded (deepreplay skips it too: tab 'parent-cyclic'); for it we count the gate-admitted vertices v whose D(v) = {v} + descendants meets a cycle.
  toolsmith_cyccount.py n CLS DEPTH BUDGET_SEC [I/K]  (repo root)   I/K: at the last level expand only frontier[I::K] (shard; the earlier levels are repeated by every shard)
Prints a tally ('edges' = all legal children; 'guard' = key-kept = the counter quoted as 'key-kept edges'; 'noguard' = key-moving; 'x-*' = tiltingPlus x guard); 'failing' lines are the key-kept tiltingPlus-false children (J != 0 steps)."""
import sys, time
from collections import Counter
sys.path.insert(0, '.')
_argv = sys.argv; sys.argv = ['x']
exec(compile(open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n"), 'h015', 'exec'))
sys.argv = _argv
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
