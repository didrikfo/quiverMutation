"""Round 130 (toolsmith), T10: F-041 (meeting in the middle, `search.meetingPoints(A, B, 3, alsoDual=True)`) at n = 10, 11 with
`refuseNonTilting=True` against the default.  Pairs = every LNA and its single two-arrow deletion that `lm.ALL_MOVES` leaves in different
orbits (as 123/skeptic_f041.py, 125/toolsmith_refuse.py): n = 10 304, n = 11 1513.

Shortcut (exact): with no `visited` set and no `deeperWhen`, the DFS tree with the option on is the default tree minus the subtrees under
refused steps.  So each side is walked ONCE with the option on and a counter on `mutation.isTiltingAtVertex`; if it returned False
nowhere, the default ball is the same set of keys.  Where it did (or with --direct), the default ball is walked too and compared.

  plan N [K]          : pair count, time K sample pairs (default 2), shards of <= 8 min for 4 cores.
  run N LO HI [--direct] : pairs LO..HI-1, one JSON line per pair to stdout (append to workshop/rounds/130/toolsmith_f041refuse_N.jsonl).
  control             : planted: E-318's parent P (n = 11, `tests/test_search_soundness._e318Parent`) and its non-tilting child X;
                        meetingPoints(P, X, 1) off vs on, counter, and the same through `run`'s code path.
  missing N SHARDS PER FILES... : indices not yet done, as SHARDS comma lists (feed to `run N 0 0 --idx LIST`; resume after timeouts).
  summary N           : totals from workshop/rounds/130/toolsmith_f041refuse_N.tsv (the run lines, as a table: jsonl -> tsv by hand); joined pairs off/on; orbits of all LNAs under ALL_MOVES + meeting joins, off/on.
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041refuse.py plan 10"""
import sys, os, time, json, copy
sys.path.insert(0, '.')
from quivermutation import search, mutation, pathAlgebra, nakayama as nk, freeMoves as fm, lnaMoves as lm

REFUSED = [0]
_orig = mutation.isTiltingAtVertex
def _counting(alg, v):
    r = _orig(alg, v)
    if not r: REFUSED[0] += 1
    return r
mutation.isTiltingAtVertex = _counting      # search.py calls mutation.isTiltingAtVertex by attribute lookup

def orbitsAndPairs(N):
    lnas, orbits = fm.derivedOrbits(N, rules=lm.ALL_MOVES, free=False, edges=True, doubles=True)
    where = {m: k for k, ms in orbits.items() for m in ms}; out = []
    for lna in lnas:
        for s, a in enumerate(lna):
            if a != 2: continue
            w = list(lna); w[s] = 0; w = tuple(w)
            if where[lna] != where[w]: out.append((lna, w))
    return lnas, where, out

def ball(alg, flag):
    REFUSED[0] = 0
    keys = set(search.quiversReachedFrom(alg, 3, alsoDual=True, refuseNonTilting=flag))
    return keys, REFUSED[0]

def onePair(N, a, b, direct):
    X, Y = (nk.LinearNakayamaAlgebra(N, list(t)) for t in (a, b)); t0 = time.time()
    (Lon, rL), (Ron, rR) = ball(X, True), ball(Y, True)
    rec = dict(a=''.join(map(str, a)), b=''.join(map(str, b)), keysL=len(Lon), keysR=len(Ron), refusedL=rL, refusedR=rR,
               meet_on=len(Lon & Ron), secs_on=round(time.time() - t0, 1))
    if direct or rL or rR:
        t1 = time.time(); Loff, Roff = ball(X, False)[0], ball(Y, False)[0]
        rec.update(off='direct', meet_off=len(Loff & Roff), removed=len(Loff - Lon) + len(Roff - Ron),
                   gained=len(Lon - Loff) + len(Ron - Roff), secs_off=round(time.time() - t1, 1))
    else:
        rec.update(off='inferred', meet_off=rec['meet_on'], removed=0, gained=0)
    return rec

def plan(N, K=2):
    t = time.time(); _, _, P = orbitsAndPairs(N); tp = time.time() - t
    step = max(1, len(P) // K); secs = []
    for i in range(0, len(P), step)[:K]:
        r = onePair(N, *P[i], False); secs.append(r['secs_on']); print('  sample', i, r, flush=True)
    per = sum(secs) / len(secs); shard = max(1, int(480 / (per * 1.3)))   # 1.3: parallel slowdown on 4 cores (notebook)
    print('PLAN n=%d pairs %d (orbits+pairs %.0fs); %.1f s/pair (option-on only); total %.1f CPU-h; shard %d pairs (<= 8 min); '
          '%d shards, %.1f h wall on 4 cores' % (N, len(P), tp, per, per * len(P) / 3600, shard, -(-len(P) // shard),
                                                 -(-len(P) // shard) / 4 * 8 / 60))

def run(N, lo, hi, direct, idx=None):
    _, _, P = orbitsAndPairs(N)
    for i in (idx if idx is not None else range(lo, min(hi, len(P)))):
        print(json.dumps(dict(onePair(N, *P[i], direct), i=i)), flush=True)

def missing(N, files, shards, per):
    """indices of pairs not yet in any of `files` (jsonl), cut into `shards` comma lists of <= `per` (resume after a timeout)."""
    _, _, P = orbitsAndPairs(N); done = set()
    for f in files:
        for line in open(f):
            r = json.loads(line); done.add((r['a'], r['b']))
    todo = [i for i, (a, b) in enumerate(P) if (''.join(map(str, a)), ''.join(map(str, b))) not in done]
    print(len(todo), file=sys.stderr)
    for k in range(shards):
        part = todo[k::shards][:per]
        if part: print(','.join(map(str, part)))

def control():
    sys.path.insert(0, 'tests'); import test_search_soundness as T
    from quivermutation import reduction
    P = T._e318Parent(); X = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(lm._copy(P), 3)); kx = search.quiverKey(X)
    for flag in (False, True):
        REFUSED[0] = 0; m = search.meetingPoints(P, X, 1, refuseNonTilting=flag)
        print('control meetingPoints(P, X, 1, refuseNonTilting=%s): meetings %d, child key met %s, refusals counted %d'
              % (flag, len(m), kx in {k for k, _, _ in m}, REFUSED[0]))
    REFUSED[0] = 0; Lon = set(search.quiversReachedFrom(P, 1, refuseNonTilting=True)); r = REFUSED[0]
    Loff = set(search.quiversReachedFrom(P, 1)); print('control ball(P, 1): on %d keys, refusals %d -> shortcut must NOT infer; '
          'off %d keys, removed %d (child among them: %s)' % (len(Lon), r, len(Loff), len(Loff - Lon), kx in Loff - Lon))

def summary(N):
    lnas, where, P = orbitsAndPairs(N); recs = {}
    rows = open('workshop/rounds/130/toolsmith_f041refuse_%d.tsv' % N).read().split('\n'); head = rows[0].split('\t')
    for line in filter(None, rows[1:]):
        r = dict(zip(head, line.split('\t'))); r = {k: (v if k in ('a', 'b', 'off') else float(v)) for k, v in r.items()}
        recs[(r['a'], r['b'])] = r
    tot = dict(pairs=len(P), done=len(recs))
    for k in ('keysL', 'keysR', 'refusedL', 'refusedR', 'removed', 'gained'): tot[k] = sum(r[k] for r in recs.values())
    tot['direct'] = sum(r['off'] == 'direct' for r in recs.values())
    tot['pairs with a refusal'] = sum((r['refusedL'] + r['refusedR']) > 0 for r in recs.values())
    tot['joined off'] = sum(r['meet_off'] > 0 for r in recs.values()); tot['joined on'] = sum(r['meet_on'] > 0 for r in recs.values())
    tot['joined off not on'] = [k for k, r in recs.items() if r['meet_off'] and not r['meet_on']]
    tot['meetings off'] = sum(r['meet_off'] for r in recs.values()); tot['meetings on'] = sum(r['meet_on'] for r in recs.values())
    for flag in ('off', 'on'):
        par = {o: o for o in set(where.values())}
        def f(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for (a, b), r in recs.items():
            if r['meet_' + flag]: par[f(where[tuple(map(int, a))])] = f(where[tuple(map(int, b))])
        tot['classes ' + flag] = len({f(o) for o in par})
    tot['orbits (moves only)'] = len(set(where.values())); print('SUMMARY n=%d' % N, tot)

if __name__ == '__main__':
    A = sys.argv[1:]
    if A[0] == 'plan': plan(int(A[1]), *(int(x) for x in A[2:3]))
    elif A[0] == 'run': run(int(A[1]), int(A[2]), int(A[3]), '--direct' in A,
                            [int(x) for x in A[A.index('--idx') + 1].split(',')] if '--idx' in A else None)
    elif A[0] == 'missing': missing(int(A[1]), A[4:], int(A[2]), int(A[3]))
    elif A[0] == 'control': control()
    elif A[0] == 'summary': summary(int(A[1]))
