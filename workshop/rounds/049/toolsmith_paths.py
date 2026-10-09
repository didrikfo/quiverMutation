"""Round 049 (toolsmith), response to referee: path recovery and parent check on top of toolsmith_tiltpath.py (same moves F / R, J = 0 etc.).
Usage (repo root): toolsmith_paths.py in.pkl cls dC mode     mode = buildball (target ball, depth 5, WITH parent pointers; once per class, cached in /tmp/tsm) | paths | parents
 paths  : for each E-154 child (one per distinct canonical key) search to depth dC against the target ball and PRINT the moves (F v / R v) on both
          sides of the meeting key: child -> meet (moves from the child) and LNA/dual -> meet (moves from the LNA, traversed backwards = inverse of a
          tilting step). Both sides are replayed with fresh move generation and the keys checked (replay ok).
 parents: the same for each distinct parent (by canonicalKey) of the children (is the parent itself tilting-path connected to an LNA?).
SLICE=lo:hi env restricts to children lo..hi-1.   A move F v / R v is the forward step at v of the algebra / of its opposite (carried back)."""
import sys, os, time, pickle, hashlib
A = sys.argv; PK, CLS, DC, MODE = A[1], A[2], int(A[3]), A[4]
sys.argv = ['x', PK, CLS, '5', str(DC), '400000', 'paths']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
cache = '/tmp/tsm/tball_c%d_d5_par.pkl' % cls     # key -> (depth, parentKey, kind, v); the move is made from the parent algebra

def buildBall():
    seen = {}; frontier = []
    for s_ in classes[base]:
        k = fingerprint.canonicalKey(s_)
        if k in seen: continue
        seen[k] = (0, None, None, None); frontier.append((k, s_))
    for d in range(1, 6):
        nxt = []
        for pk_, alg in frontier:
            for kind, v, c in moves(alg):
                k = fingerprint.canonicalKey(c)
                if k is None or k in seen: continue
                seen[k] = (d, pk_, kind, v); nxt.append((k, c))
        print('ball level', d, len(nxt), flush=True); frontier = nxt
    return seen

if MODE == 'buildball':
    t0 = time.time(); pickle.dump(buildBall(), open(cache, 'wb')); print('built cls', cls, '%.0fs' % (time.time() - t0)); sys.exit()
tpar = pickle.load(open(cache, 'rb')); tball = {k: v[0] for k, v in tpar.items()}
seedkeys = {}
for i, a_ in enumerate(classes[base]): seedkeys.setdefault(fingerprint.canonicalKey(a_), i)
def tag(alg): return hashlib.md5(repr(fingerprint.canonicalKey(alg)).encode()).hexdigest()[:6]

def reachPath(start, depth):
    k0 = fingerprint.canonicalKey(start); par = {k0: None}; objs = {k0: start}; frontier = [(k0, start)]
    if k0 in tball: return k0, par, objs
    for d in range(1, depth + 1):
        nxt = []; best = None
        for pk_, alg in frontier:
            for kind, v, c in moves(alg):
                k = fingerprint.canonicalKey(c)
                if k is None or k in par: continue
                par[k] = (pk_, kind, v); objs[k] = c; nxt.append((k, c))
                if k in tball and (best is None or d + tball[k] < best[0]): best = (d + tball[k], k)
        frontier = nxt
        if best: return best[1], par, objs
    return None, par, objs

def follow(start, seq):
    """apply the moves seq = [(kind, v, expected key)] from start with fresh move generation; return the end key or None."""
    cur = start
    for kind, v, kx in seq:
        cand = [c for kd, vv, c in moves(cur) if kd == kind and vv == v and fingerprint.canonicalKey(c) == kx]
        if not cand: return None
        cur = cand[0]
    return fingerprint.canonicalKey(cur)

def pathOf(start, depth):
    meet, par, objs = reachPath(start, depth)
    if meet is None: return None
    cs = []; k = meet
    while par[k] is not None: pk_, kind, v = par[k]; cs.append((kind, v, k)); k = pk_
    cs.reverse()
    ts = []; k = meet
    while tpar[k][1] is not None: d_, pk_, kind, v = tpar[k]; ts.append((kind, v, k)); k = pk_
    ts.reverse(); seed = k
    ok = follow(start, cs) == meet and follow(classes[base][seedkeys[seed]], ts) == meet
    return ('ok' if ok else 'REPLAY-FAILED', cs, ts, seedkeys[seed])

recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
LO, HI = (int(x) for x in os.environ.get('SLICE', '0:%d' % len(recs)).split(':'))
done = {}
f = lambda sq: ' '.join('%s%d' % (kd, v) for kd, v, _ in sq) or '-'
for i, r in enumerate(recs):
    if not LO <= i < HI: continue
    obj = r['childObj'] if MODE == 'paths' else r['parentObj']
    t = tag(obj)
    if t in done: print(MODE, i, 'key', t, 'same key as', done[t]); continue
    done[t] = i; t1 = time.time(); res = pathOf(obj, DC)
    if res is None: print(MODE, i, 'key', t, 'parentdepth', r['depth'], 'MISS at dC', DC, '%.0fs' % (time.time() - t1), flush=True); continue
    print(MODE, i, 'key', t, 'parentdepth', r['depth'], 'replay', res[0], 'total', len(res[1]) + len(res[2]), '| child->meet:', f(res[1]), '| LNA/dual #%d ->meet:' % res[3], f(res[2]), '%.0fs' % (time.time() - t1), flush=True)
print('DONE', MODE, 'cls', cls, 'stats', stats)
