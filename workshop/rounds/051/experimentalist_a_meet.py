"""Round 051: n = 10 polynomial group A (orbits 03345000 [69 members], 33460000 [42]; E-032: they merge at one-sided
depth 6).  Meet in the middle: search.quiversReachedFrom(member, 3, alsoDual=True) for every member, hash-join across the
two orbits, shortest total, replay both halves edge by edge (gate, tiltingPlus J = 0, key kept).  Members and orbit
names are read from the merges.py checkpoint of this round.
Usage: timeout 10m .venv/bin/python workshop/rounds/051/experimentalist_a_meet.py   (prints; no output file)"""
import sys, json, copy, time, multiprocessing
sys.path.insert(0, '.')
from quivermutation import nakayama as nk, mutation, procedure, reduction, search, pathAlgebra
exec(compile(open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0], 'h015', 'exec'))
CK = 'workshop/rounds/051/experimentalist_merges10.jsonl'
members = {}
for l in open(CK):
    r = json.loads(l)
    if r['orbit'] in ('03345000', '33460000'): members[tuple(r['member'])] = r['orbit']

def reach(m):
    return m, search.quiversReachedFrom(nk.LinearNakayamaAlgebra(10, list(m)), 3, True)

def replay(alg0, path):
    dual = bool(path) and path[0] < 0
    alg = pathAlgebra.dualPathAlgebra(alg0) if dual else copy.deepcopy(alg0)
    base = search._coxeterKeyOrNone(alg); out = []
    for s in path:
        v = abs(s)
        gate = mutation.mutationIsPossibleAtVertex(alg, v)
        J0 = bool(tiltingPlus(alg.quiver, procedure.relationsFrom(alg), v))
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(alg), v))
        out.append((v, gate, J0, search._coxeterKeyOrNone(alg) in (None, base)))
    return out

if __name__ == '__main__':
    t0 = time.time()
    with multiprocessing.Pool(4) as pool:
        R = dict(pool.map(reach, sorted(members), chunksize=1))
    print('members', len(R), 'reach sets done in %.0fs' % (time.time() - t0), flush=True)
    byKey = {}
    for m, d in R.items():
        for k, p in d.items(): byKey.setdefault(k, []).append((m, p))
    best = []
    for k, lst in byKey.items():
        A = [x for x in lst if members[x[0]] == '03345000']; B = [x for x in lst if members[x[0]] == '33460000']
        for a in A:
            for b in B: best.append((len(a[1]) + len(b[1]), a[0], b[0], a[1], b[1]))
    best.sort(key=lambda x: x[0])
    print('meeting triples (key, A-member, B-member):', len(best), 'shortest total', best[0][0] if best else None)
    tot = [0, 0, 0, 0]; shown = 0
    for L, a, b, pa, pb in best[:3] if best else []:
        ea = replay(nk.LinearNakayamaAlgebra(10, list(a)), pa); eb = replay(nk.LinearNakayamaAlgebra(10, list(b)), pb)
        es = ea + eb
        print(''.join(map(str, a)), pa, '|', ''.join(map(str, b)), pb, 'edges', len(es),
              'gate', sum(e[1] for e in es), 'J0', sum(e[2] for e in es), 'key', sum(e[3] for e in es), flush=True)
