"""Round 038 (toolsmith): close the n = 6 LNA derived classes (key-preserving BFS, as rounds/035/toolsmith_closure.py) and tabulate (d_i, dim J_i) over every gate-admitted (algebra, v, i).
Usage: toolsmith_n6close.py --plan SECONDS        per-class size sample, no verdict
       toolsmith_n6close.py [--budget-hours H] [--only IDX]   full closure; exit 2 if budget spent
Also reports whether the 16 key-coinciding fans of rounds/037/theorist_d3_hits_n6.json (E-134) have their canonicalKey in a closed class."""
import sys, time, json, argparse
from collections import Counter, defaultdict
ap_ = argparse.ArgumentParser(); ap_.add_argument('--plan', type=float, default=None, help='seconds per class, sizing only')
ap_.add_argument('--budget-hours', type=float, default=0.15); ap_.add_argument('--only', type=int, default=None)
args = ap_.parse_args()
sys.path.insert(0, '.'); sys.argv = ['x', 'none']
g = {'__name__': 'bd'}
exec(compile(open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0], 'bd', 'exec'), g)
nk, pathAlgebra, search, fingerprint, mutation, procedure, reduction, ap, nx, perI, build = (g[k] for k in
    'nk pathAlgebra search fingerprint mutation procedure reduction ap nx perI build'.split())
n = 6
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
print('n = 6: %d LNA key classes, sizes %s' % (len(order), [len(classes[k]) for k in order]), flush=True)
def dims(alg, v):
    Q = alg.quiver; rels = procedure.relationsFrom(alg); res = {}; J = perI(alg, v)
    for i in Q.nodes:
        if i == v: continue
        P = ap.allPathsBetween(Q, i, v)
        if not P: continue
        res[i] = (len(P) - len(ap.idealBasis(Q, rels, i, v)), J.get(i, 0))
    return res
hits = json.load(open('workshop/rounds/037/theorist_d3_hits_n6.json'))
hitkeys = {}
for h in hits:
    A = build([tuple(x) for x in h[0]], h[1], list(range(1, 7)))
    hitkeys[fingerprint.canonicalKey(A)] = (search._coxeterKeyOrNone(A), h[3])
budget = args.plan if args.plan is not None else args.budget_hours * 3600
tab = Counter(); tabkeyed = {}; summary = []; spent_any = False
for idx, key in enumerate(order):
    if args.only is not None and idx != args.only: continue
    seen = set(); frontier = []; nexp = 0; lvl = 0; t0 = time.time(); spent = False; ctab = Counter(); maxd = 0; edges = 0
    def add(alg):
        ck = fingerprint.canonicalKey(alg)
        if ck is None: ck = ('none', id(alg))
        if ck in seen: return
        seen.add(ck); frontier.append(alg)
    for a in classes[key]: add(a)
    while frontier and not spent:
        cur, frontier[:] = frontier[:], []; lvl += 1
        for i, alg in enumerate(cur):
            if time.time() - t0 > budget: frontier.extend(cur[i:]); spent = True; break
            nexp += 1
            if list(nx.simple_cycles(alg.quiver)): continue
            for v in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, v): continue
                edges += 1
                for dj in dims(alg, v).values(): ctab[dj] += 1
                raw = mutation.quiverMutationAtVertex(alg, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if search._coxeterKeyOrNone(ch) == key: add(ch)
        print('  class %d level %d: seen %d frontier %d (%.0fs)' % (idx, lvl, len(seen), len(frontier), time.time() - t0), flush=True)
    closed = not frontier
    h = [hk for hk in hitkeys if hk in seen]
    summary.append((idx, len(classes[key]), len(seen), nexp, edges, closed, len(h), hitkeys and sum(1 for v in hitkeys.values() if v[0] == key)))
    print('CLASS idx %d LNAs %d seen %d expanded %d gate-edges %d %s | 16 hits in class: %d (of %d with this key) | (d,J) %s' % (
        idx, len(classes[key]), len(seen), nexp, edges, 'CLOSED' if closed else 'CAPPED', len(h), summary[-1][-1], dict(sorted(ctab.items()))), flush=True)
    tab.update(ctab); spent_any |= not closed
print('TOTAL (d,J) over all gate-admitted (alg, v, i):', dict(sorted(tab.items())))
print('RESULT', 'ALL CLOSED' if not spent_any else 'SOME CAPPED')
sys.exit(2 if spent_any else 0)
