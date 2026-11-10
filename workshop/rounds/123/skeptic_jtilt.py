# Round 123 (skeptic): verbatim copy of workshop/rounds/121/maverick_jtilt.py (E-317 classifier), as the round-123 call asks; skeptic_mergereplay.py uses the same step classification.
# Copied verbatim from the skeptic referee's review of round 121 (workshop/rounds/121/maverick.review.md); author: skeptic. Usage: .venv/bin/python -W ignore workshop/rounds/121/maverick_jtilt.py ROOT DEPTH
# skeptic r121 scratch: classify every gate-admitted step in a guarded BFS ball from ROOT (and dual) by (J, key, tiltingPlus)
import sys, time, copy, os
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
root, D = ARGV[1], int(ARGV[2])
r = nk.LinearNakayamaAlgebra.fromClassName(root)
t0 = time.time(); seen = set(); frontier = []; tab = Counter()
def add(alg):
    k = fingerprint.canonicalKey(alg)
    if k is not None:
        if k in seen: return
        seen.add(k)
    frontier.append(alg)
add(copy.deepcopy(r)); add(pathAlgebra.dualPathAlgebra(r))
pk0 = search._coxeterKeyOrNone(r)
for level in range(D):
    cur, frontier[:] = frontier[:], []
    for alg in cur:
        if list(nx.simple_cycles(alg.quiver)): continue
        V = sorted(alg.quiver.nodes)
        for v in V:
            if not mutation.mutationIsPossibleAtVertex(alg, v): continue
            raw = mutation.quiverMutationAtVertex(alg, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): tab['illegal'] += 1; continue
            ch = reduction.reducePathAlgebra(raw)
            ck = search._coxeterKeyOrNone(ch)
            if ck is None: tab['cyclic-child'] += 1; continue
            J = bool(perI(alg, v)); moved = ck != pk0
            tp = bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v))
            tab[('J!=0' if J else 'J=0', 'moved' if moved else 'kept', 'tilting' if tp else 'NOT-tilting')] += 1
            if not moved: add(ch)
    print('level', level + 1, 'distinct', len(seen), dict(tab), '%.0fs' % (time.time() - t0), flush=True)
