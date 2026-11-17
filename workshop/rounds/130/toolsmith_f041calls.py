"""Round 130 (toolsmith), referee answer: the F-041 walks of toolsmith_f041refuse.py with two more counters on
`mutation.isTiltingAtVertex`: calls (every gate-admitted step reaches it, option on) and, per call, the step's
d = max_i dim e_i A e_v (paths i ~> v modulo I, i != v; v the mutated vertex) and the out-degree of v.
Referee's thinness reading: if d <= 1 at a gate-admitted step, J = 0 follows from the gate; so J != 0 needs d >= 2.

  run N --idx LIST [--direct] : one JSON line per pair (same fields as toolsmith_f041refuse.py run, plus calls / d2 / out2 per side).
  control                     : d at E-318's refused step (vertex 3 of `_e318Parent`).
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041calls.py run 11 --idx 5,131 --direct"""
import sys, json, time, importlib.util
sys.path.insert(0, '.')
spec = importlib.util.spec_from_file_location('tfr', 'workshop/rounds/130/toolsmith_f041refuse.py')
tfr = importlib.util.module_from_spec(spec); spec.loader.exec_module(tfr)     # installs the refusal counter
from quivermutation import search, mutation, procedure, nakayama as nk, arrowPaths as ap

C = dict(calls=0, d2=0, out2=0, dmax=0)
def maxDim(alg, v):
    q = alg.quiver; rel = procedure.relationsFrom(alg); d = 0
    for i in q.nodes:
        if i == v: continue
        c = ap.allPathsBetween(q, i, v)
        if c: d = max(d, len(c) - len(ap.idealBasis(q, rel, i, v)))
    return d
_inner = mutation.isTiltingAtVertex
def _wrap(alg, v):
    C['calls'] += 1; d = maxDim(alg, v); C['dmax'] = max(C['dmax'], d)
    if d >= 2: C['d2'] += 1
    if len(ap.arrowsOutOf(alg.quiver, v)) >= 2: C['out2'] += 1
    return _inner(alg, v)
mutation.isTiltingAtVertex = _wrap

def ball(alg, flag):
    for k in C: C[k] = 0
    tfr.REFUSED[0] = 0
    keys = set(search.quiversReachedFrom(alg, 3, alsoDual=True, refuseNonTilting=flag))
    return keys, tfr.REFUSED[0], dict(C)

def onePair(N, a, b, direct):
    X, Y = (nk.LinearNakayamaAlgebra(N, list(t)) for t in (a, b)); t0 = time.time()
    (Lon, rL, cL), (Ron, rR, cR) = ball(X, True), ball(Y, True)
    rec = dict(a=''.join(map(str, a)), b=''.join(map(str, b)), keysL=len(Lon), keysR=len(Ron), refusedL=rL, refusedR=rR,
               meet_on=len(Lon & Ron), secs_on=round(time.time() - t0, 1),
               calls=cL['calls'] + cR['calls'], d2=cL['d2'] + cR['d2'], out2=cL['out2'] + cR['out2'], dmax=max(cL['dmax'], cR['dmax']))
    if direct or rL or rR:
        t1 = time.time(); Loff, Roff = ball(X, False)[0], ball(Y, False)[0]
        rec.update(off='direct', meet_off=len(Loff & Roff), removed=len(Loff - Lon) + len(Roff - Ron),
                   gained=len(Lon - Loff) + len(Ron - Roff), secs_off=round(time.time() - t1, 1))
    else:
        rec.update(off='inferred', meet_off=rec['meet_on'], removed=0, gained=0)
    return rec

def control():
    sys.path.insert(0, 'tests'); import test_search_soundness as T
    P = T._e318Parent()
    for v in P.quiver.nodes:
        if mutation.mutationIsPossibleAtVertex(P, v):
            print('control E-318 parent, vertex %s: d = %d, out-degree %d, tilting %s'
                  % (v, maxDim(P, v), len(ap.arrowsOutOf(P.quiver, v)), _inner(P, v)))

if __name__ == '__main__':
    A = sys.argv[1:]
    if A[0] == 'control': control()
    else:
        N = int(A[1]); _, _, P = tfr.orbitsAndPairs(N)
        for i in [int(x) for x in A[A.index('--idx') + 1].split(',')]:
            print(json.dumps(dict(onePair(N, *P[i], '--direct' in A), i=i)), flush=True)
