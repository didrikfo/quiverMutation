"""Round 045 (experimentalist): T10 guard audit (b).  Do the n = 8 F-041 merges (the 10 single two-arrow deletions the
known moves leave open, tests/test_other_families.py) depend on a J != 0 step?
Each pair (lna, lna with one 2-arrow relation deleted) is met in the middle (both algebras + relation duals, depth D per side,
nodes keyed by search.quiverKey as in search.meetingPoints) in three modes:
  guard : gate + Coxeter-key guard (the recorded search).  Also tallies per expanded edge whether tiltingPlus holds.
  tilt  : gate + tiltingPlus (J = 0) only, key guard OFF; also tallies whether the child key equals the base key.
  both  : gate + tiltingPlus + key guard (control: must be inside guard).
Usage: .venv/bin/python workshop/rounds/045/experimentalist_t10b.py [--plan] [--depth 3] [--pairs 0,1,..]"""
import sys, time, argparse, copy
from collections import Counter
sys.path.insert(0, '.')
import networkx as nx
from quivermutation import nakayama as nk, mutation, procedure, reduction, search, pathAlgebra
from quivermutation import arrowPaths as ap
from quivermutation import freeMoves as fm, lnaMoves as lm
exec(compile(open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0], 'h015', 'exec'))

ap_ = argparse.ArgumentParser(); ap_.add_argument('--plan', action='store_true'); ap_.add_argument('--depth', type=int, default=3)
ap_.add_argument('--pairs', default='all'); ap_.add_argument('--modes', default='guard,tilt,both')
A = ap_.parse_args()

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

def children(alg, mode, base, st):
    """yield (child, v) ; st counts."""
    if list(nx.simple_cycles(alg.quiver)): return
    for v in sorted(alg.vertices()):
        if not mutation.mutationIsPossibleAtVertex(alg, v): continue
        st['gate'] += 1
        J0 = tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v)
        if mode in ('tilt', 'both') and not J0: st['J!=0 dropped'] += 1; continue
        raw = mutation.quiverMutationAtVertex(copy.deepcopy(alg), v)
        if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
        ch = reduction.reducePathAlgebra(raw)
        k = search._coxeterKeyOrNone(ch)
        keep = (k is None or k == base)
        if mode in ('guard', 'both') and not keep: st['key moved dropped'] += 1; continue
        st['J0' if J0 else 'J!=0'] += 1
        if mode == 'tilt': st['tilt step key kept' if keep else 'tilt step KEY MOVED'] += 1
        yield ch, v

def reach(alg0, depth, mode, st, dual=False):
    start = pathAlgebra.dualPathAlgebra(alg0) if dual else alg0
    base = search._coxeterKeyOrNone(start)
    def tk(x): return search.quiverKey(pathAlgebra.dualPathAlgebra(x) if dual else x)
    found = {}; seen = {}; fr = [start]
    k0 = tk(start)
    if k0 is not None: found[k0] = 0
    for d in range(1, depth + 1):
        nxt = []
        for alg in fr:
            for ch, v in children(alg, mode, base, st):
                sk = search.quiverKey(ch)
                if sk is not None:
                    if sk in seen: continue
                    seen[sk] = 1
                nxt.append(ch)
                k = tk(ch)
                if k is not None and k not in found: found[k] = d
        fr = nxt
    return found, len(seen)

P = pairs()
print('pairs: %d' % len(P), flush=True)
sel = range(len(P)) if A.pairs == 'all' else [int(x) for x in A.pairs.split(',')]
if A.plan:
    for i in sel: print(i, P[i]); 
    sys.exit()
for mode in A.modes.split(','):
    print('=== mode', mode, 'depth', A.depth, flush=True)
    tot = Counter()
    for i in sel:
        a, b = P[i]; t0 = time.time(); st = Counter(); res = []
        for alg in (nk.LinearNakayamaAlgebra(8, list(a)), nk.LinearNakayamaAlgebra(8, list(b))):
            f = {}
            for dual in (False, True):
                g, ns = reach(alg, A.depth, mode, st, dual)
                for k, d in g.items():
                    if k not in f or d < f[k]: f[k] = d
            res.append(f)
        sh = set(res[0]) & set(res[1])
        best = min((res[0][k] + res[1][k] for k in sh), default=None)
        print('pair %d %s -> %s: meetings %d, shortest total %s, %.0fs, edges %s' % (i, ''.join(map(str, a)), ''.join(map(str, b)), len(sh), best, time.time() - t0, dict(st)), flush=True)
        tot.update(st)
    print('TOTAL edges', dict(tot), flush=True)
