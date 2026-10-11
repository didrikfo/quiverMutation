"""Round 035 (toolsmith): size / run the n = 7 BFS closure of LNA key classes against the both-die targets (E-126).
Usage: toolsmith_closure.py plan CLASSIDX SECONDS      -- per-level growth + rate for SECONDS, no verdict
       toolsmith_closure.py run  CLASSIDX [--budget-hours H]  -- full BFS; exit 2 if budget spent, 0 if closed
Same BFS as rounds/033/experimentalist_bothdie.py reach (acyclic-quiver algebras expanded, children kept iff key == class key)."""
import sys, time, argparse
from collections import Counter
mode = sys.argv[1]; cidx = int(sys.argv[2]); n = 7
p = argparse.ArgumentParser(); p.add_argument('--budget-hours', type=float, default=0.15)
args, rest = p.parse_known_args(sys.argv[3:])
secs = float(rest[0]) if (mode == 'plan' and rest) else args.budget_hours * 3600
sys.argv = ['x', 'none']
g = {'__name__': 'bd'}
exec(compile(open('workshop/rounds/033/experimentalist_bothdie.py').read(), 'bd', 'exec'), g)
nk, pathAlgebra, search, fingerprint, mutation, procedure, reduction, ap, nx = (g[k] for k in
    'nk pathAlgebra search fingerprint mutation procedure reduction ap nx'.split())
classes = {}
for lna in nk.LinearNakayamaAlgebra.allOfLength(n):
    for alg in (lna, pathAlgebra.dualPathAlgebra(lna)):
        classes.setdefault(search._coxeterKeyOrNone(alg), []).append(alg)
order = sorted(classes, key=lambda k: (len(classes[k]), str(k)))
key = order[cidx]
targets = set()
for cname, A, v, gate, kd, k in g['gen'](n, classes):
    if A is None or not gate or not kd or k != key: continue
    targets.add(fingerprint.canonicalKey(A))
print('class', cidx, 'size', len(classes[key]), 'targets', len(targets), flush=True)
seen = set(); frontier = []
def add(alg):
    ck = fingerprint.canonicalKey(alg)
    if ck is not None:
        if ck in seen: return
        seen.add(ck)
    frontier.append(alg)
for a in classes[key]: add(a)
t0 = time.time(); nexp = 0; lvl = 0; spent = False
while frontier and not spent:
    cur, frontier[:] = frontier[:], []; lvl += 1; tl = time.time(); ne = 0
    for i, alg in enumerate(cur):
        if time.time() - t0 > secs: frontier.extend(cur[i:]); spent = True; break
        nexp += 1; ne += 1
        if list(nx.simple_cycles(alg.quiver)): continue
        for v in sorted(alg.vertices()):
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if search._coxeterKeyOrNone(ch) == key: add(ch)
    print('level %d: expanded %d of %d (%.0fs, %.1f exp/s) seen %d frontier %d hits %d' % (
        lvl, ne, len(cur), time.time() - tl, ne / max(time.time() - tl, 1e-9), len(seen), len(frontier), len(targets & seen)), flush=True)
print('RESULT', 'CLOSED' if not frontier else 'CAPPED', 'expanded', nexp, 'seen', len(seen), 'hits', len(targets & seen), 'of', len(targets), '%.0fs' % (time.time() - t0))
sys.exit(0 if not frontier else 2)
