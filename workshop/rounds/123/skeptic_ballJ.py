"""Round 123 (skeptic), response to referee item 2 and 5.  For each F-041 pair (as skeptic_f041.py) at n = N, walk the WHOLE depth-3
key-guarded gate ball that search.quiversReachedFrom(., 3, alsoDual=True) walks, on both ends (LNA and its opposite each), with a visitor,
and classify EVERY walk edge (parent node -> child node, vertex = last entry of the child's path): J != 0 by `perI`, `tiltingPlus`.
An edge is a walk edge: the search took it (gate admitted, key kept, no discard), so it is exactly the pool F-041 witnesses are drawn from.
Edges are counted per route (the undeduped walk), and also as distinct (parent quiverKey, vertex, side) triples.
Also (item 5): number of meeting keys per pair (shared quiverKeys of the two found dicts, as meetingPoints), and that the shortest
path stored for the first meeting key replays to a node whose quiverKey is that key on both sides.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_ballJ.py N [first last]"""
import sys, time, copy, collections
ARG0 = sys.argv[:]; N = int(ARG0[1])
sys.path.insert(0, '.'); sys.argv = ['x', 'none']
_f = open('workshop/rounds/123/skeptic_f041.py').read()
exec(compile(_f.split("P = pairs();")[0].replace("N = int(MYARGV[1])", "N = %d" % N), 'f041', 'exec'))   # machinery, INV, pairs(), half()

def ball(alg0):
    """found (key -> shortest signed path, as quiversReachedFrom) and the edge counts of both walks."""
    found = {}; cnt = collections.Counter(); distinct = {}
    for dual in (False, True):
        start = pathAlgebra.dualPathAlgebra(alg0) if dual else copy.deepcopy(alg0); nodes = {}
        def visit(q, mv, dual=dual, nodes=nodes):
            mv = tuple(mv); nodes[mv] = copy.deepcopy(q)
            k = search.quiverKey(pathAlgebra.dualPathAlgebra(q) if dual else q)
            if k is not None:
                p = [(-s if dual else s) for s in mv]
                if k not in found or len(p) < len(found[k]): found[k] = p
            if mv:
                par = nodes[mv[:-1]]; v = mv[-1]
                J = bool(perI(par, v)); tp = bool(tiltingPlus(par.quiver, procedure.relationsFrom(par), v))
                cnt['edges'] += 1; cnt['J!=0'] += J; cnt['tp false'] += (not tp)
                dk = (search.quiverKey(par), v, dual)
                if dk[0] is not None: distinct[dk] = (J, tp)
        search.mutationSearchDepthFirst(start, 3, [], 'meet', printOutput=False, visitor=visit)
    cnt['distinct edges'] = len(distinct); cnt['distinct J!=0'] = sum(j for j, _ in distinct.values())
    return found, cnt

P = pairs(); lo, hi = (int(ARG0[2]), int(ARG0[3])) if len(ARG0) > 3 else (0, len(P))
print('N', N, 'pairs', len(P), 'running', lo, hi, flush=True)
tot = collections.Counter(); t0 = time.time(); flagged = []
for i in range(lo, min(hi, len(P))):
    a, b = P[i]; X, Y = (nk.LinearNakayamaAlgebra(N, list(t)) for t in (a, b))
    fX, cX = ball(X); fY, cY = ball(Y); c = cX + cY
    shared = sorted(set(fX) & set(fY), key=lambda k: len(fX[k]) + len(fY[k]))
    tot['pairs'] += 1; tot['meeting'] += bool(shared); tot['meeting keys'] += len(shared); tot['pairs with >1 meeting'] += len(shared) > 1
    for k in c: tot[k] += c[k]
    keyok = None
    if shared:
        k = shared[0]; _, el = half(X, fX[k]); _, er = half(Y, fY[k])
        keyok = search.quiverKey(el) == k == search.quiverKey(er); tot['first meet key check ok'] += keyok
    if c['J!=0'] or c['tp false'] or (shared and not keyok): flagged.append((i, a, b, dict(c), keyok))
    print('%d %s %s  meetings %d keyok %s  edges %d (distinct %d)  J!=0 %d  tp false %d  %.0fs' % (i, ''.join(map(str, a)), ''.join(map(str, b)),
          len(shared), keyok, c['edges'], c['distinct edges'], c['J!=0'], c['tp false'], time.time() - t0), flush=True)
for k in ('edges', 'J!=0', 'tp false', 'distinct edges', 'distinct J!=0'): tot[k] += 0
print('TOTAL', {k: tot[k] for k in tot}); print('FLAGGED', flagged); print('time %.0fs' % (time.time() - t0))
