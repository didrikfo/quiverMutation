"""Round 047 (toolsmith): witness paths of the 10 F-041 n = 8 merges and per-edge tiltingPlus.
For each pair (lna, lna with one 2-arrow relation deleted; as experimentalist_t10b.pairs) take
search.meetingPoints(A, B, depth, alsoDual=True) (shortest total first), and replay each recorded
half-path from its own start (the algebra, or its opposite when the steps are negative) checking
gate, tiltingPlus (J = 0) and the Coxeter key at every edge; the replay must end on the meeting key.
Edges are tallied in the forward direction from each side's start (the B half is NOT inverted).
Usage: .venv/bin/python workshop/rounds/047/toolsmith_witness.py [--depth 3] [--pairs 0,1] [--all]"""
import sys, time, argparse, copy
sys.path.insert(0, '.')
from quivermutation import nakayama as nk, mutation, procedure, reduction, search, pathAlgebra
from quivermutation import freeMoves as fm, lnaMoves as lm
exec(compile(open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0], 'h015', 'exec'))

argp = argparse.ArgumentParser(); argp.add_argument('--depth', type=int, default=3)
argp.add_argument('--pairs', default='all'); argp.add_argument('--all', action='store_true', help='replay every shortest-total meeting, not only the first')
A = argp.parse_args()

def pairs():
    lnas, orbits = fm.derivedOrbits(8, rules=lm.ALL_MOVES, free=False, edges=True, doubles=True)
    where = {m: k for k, ms in orbits.items() for m in ms}
    out = []
    for lna in lnas:
        for s, a in enumerate(lna):
            if a != 2: continue
            w = list(lna); w[s] = 0; w = tuple(w)
            if where[lna] != where[w]: out.append((lna, w))
    return out

def replay(alg0, path, target):
    """-> list of (v, gate, J0, keyKept) per edge, and whether the end quiver key == target."""
    dual = bool(path) and path[0] < 0
    alg = pathAlgebra.dualPathAlgebra(alg0) if dual else copy.deepcopy(alg0)
    base = search._coxeterKeyOrNone(alg); edges = []
    for s in path:
        v = abs(s); assert (s < 0) == dual
        gate = mutation.mutationIsPossibleAtVertex(alg, v)
        J0 = bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v))
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(alg), v))
        edges.append((v, gate, J0, search._coxeterKeyOrNone(alg) in (None, base)))
    end = search.quiverKey(pathAlgebra.dualPathAlgebra(alg) if dual else alg)
    return edges, end == target

def inverseEdges(alg0, path):
    """Inverse of each forward step alg_i -> alg_{i+1} at v, taken as the forward step at v of the OPPOSITE of alg_{i+1}
    (quiver mutation is an involution and commutes with op). -> list of (v, gate, J0) in the frame the path ran in."""
    dual = bool(path) and path[0] < 0
    alg = pathAlgebra.dualPathAlgebra(alg0) if dual else copy.deepcopy(alg0); out = []
    for s in path:
        v = abs(s)
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(alg), v))
        op = pathAlgebra.dualPathAlgebra(alg)
        out.append((v, mutation.mutationIsPossibleAtVertex(op, v), bool(tiltingPlus(op.quiver, procedure.relationsFrom(op), v))))
    return out

P = pairs(); sel = range(len(P)) if A.pairs == 'all' else [int(x) for x in A.pairs.split(',')]
print('pairs', len(P), flush=True)
tot = dict(edges=0, J0=0, gate=0, key=0); inv = dict(edges=0, J0=0, gate=0); depths = []
for i in sel:
    a, b = P[i]; t0 = time.time()
    X, Y = (nk.LinearNakayamaAlgebra(8, list(t)) for t in (a, b))
    ms = search.meetingPoints(X, Y, A.depth, alsoDual=True)
    if not ms: print(i, 'NO MEETING'); continue
    best = len(ms[0][1]) + len(ms[0][2]); tied = [m for m in ms if len(m[1]) + len(m[2]) == best]
    for j, (key, pl, pr) in enumerate(tied if A.all else tied[:1]):
        el, okl = replay(X, pl, key); er, okr = replay(Y, pr, key)
        es = el + er
        row = (i, ''.join(map(str, a)), ''.join(map(str, b)), len(pl), len(pr), best, len(tied),
               sum(e[1] for e in es), sum(e[2] for e in es), sum(e[3] for e in es), len(es), okl and okr)
        print('pair %d %s->%s halves %d+%d total %d (tied meetings %d): gate %d J0 %d keykept %d of %d edges; replay ends on meeting key: %s' % row,
              'paths', pl, pr, '%.0fs' % (time.time() - t0), flush=True)
        ie = inverseEdges(Y, pr); inv['edges'] += len(ie); inv['gate'] += sum(e[1] for e in ie); inv['J0'] += sum(e[2] for e in ie)
        print('   inverse B-half edges: gate %d J0 %d of %d' % (sum(e[1] for e in ie), sum(e[2] for e in ie), len(ie)), flush=True)
        for k, v in zip(('gate', 'J0', 'key'), row[7:10]): tot[k] += v
        tot['edges'] += len(es)
    depths.append(best)
print('INVERSE B-half', inv); print('TOTAL', tot, 'shortest totals', sorted(depths))
