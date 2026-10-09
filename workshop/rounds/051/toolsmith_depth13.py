"""Round 051 (toolsmith, from 050's toolsmith_depth7.py: env TBD target-ball file depth, TD truncation, CAP key cap), T10 (i): depth-7 child ball for the E-155 misses, sharded.  Same moves as rounds/049/toolsmith_tiltpath.py (F and R, J = 0,
tiltingPlus, gate, legal, vertex-preserving, class key kept), same canonical-key meet against the class target ball (LNAs + duals, depth TD, cache
/tmp/tsm/tball_c<cls>_d5.pkl built by 049's script).  A meeting at child depth d and target depth t = a tilting-only path of length d + t.
Usage (repo root; every command is bounded, see --plan):
  toolsmith_depth7.py plan  in.pkl cls childIdx D          print the sizes known so far, nothing expanded
  toolsmith_depth7.py front in.pkl cls childIdx D tag      BFS the child ball to depth D-1 (hits checked on the way), pickle frontier -> /tmp/tsm/fr_<tag>.pkl
  toolsmith_depth7.py shard tag lo:hi                      expand frontier[lo:hi] (level D), check each new key against the target ball -> /tmp/tsm/sh_<tag>_<lo>.pkl
  toolsmith_depth7.py merge tag                            combine the shards: nodes per level, distinct level-D keys, hits (min total)
env TD=<t> truncates the target ball to depth <= t (positive control: a child known to meet at 6 + 5 must be found at 7 + 4)."""
import sys, os, time, pickle, hashlib, glob
A = sys.argv; MODE = A[1]
if MODE in ('plan', 'front'):
    PK, CLS, IDX, D = A[2], int(A[3]), int(A[4]), int(A[5]); TAG = A[6] if len(A) > 6 else None
elif MODE == 'shard': TAG = A[2]; LO, HI = (int(x) for x in A[3].split(':'))
else: TAG = A[2]
if MODE in ('shard', 'merge'):
    fr = pickle.load(open('/tmp/tsm/fr_%s.pkl' % TAG, 'rb')) if MODE == 'shard' else pickle.load(open('/tmp/tsm/fr_%s.meta.pkl' % TAG, 'rb'))
    PK, CLS = fr['pk'], fr['cls']
sys.argv = ['x', PK, str(CLS), '5', '6', '400000', 'none']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
TD = int(os.environ.get('TD', '6')); TBD = int(os.environ.get('TBD', '6'))   # TD truncates; TBD = depth of the cached target ball file
CAP = int(os.environ.get('CAP', '720'))   # canonicalKey relabeling cap (default 720 = fingerprint.DEFAULT_CAP); round 051: raise it
if CAP != 720:
    _ck = fingerprint.canonicalKey
    fingerprint.canonicalKey = lambda a, cap=CAP, gauge=True: _ck(a, cap=cap, gauge=gauge)
def loadball():
    tb, tl, tc = pickle.load(open('/tmp/tsm/tball_c%d_d%d.pkl' % (CLS, TBD), 'rb'))
    return {k: d for k, d in tb.items() if d <= TD}
def h(k): return hashlib.md5(repr(k).encode()).hexdigest()[:12]

if MODE == 'plan':
    tb = loadball(); recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
    print('target ball depth<=%d: %d keys; child %d key %s; D %d; 049 level sizes for the misses: [1,10,56,229,982,4237,18034] (c1 f7abe9), c2 child 6 5703 nodes at d6' % (TD, len(tb), IDX, h(fingerprint.canonicalKey(recs[IDX]['childObj'])), D))
if MODE == 'front':
    tb = loadball(); recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
    start = recs[IDX]['childObj']; k0 = fingerprint.canonicalKey(start); seen = {k0: 0}; frontier = [start]; levels = [1]; best = None; t1 = time.time()
    if k0 in tb: best = (tb[k0], 0, tb[k0])
    for d in range(1, D):
        nxt = []
        for alg in frontier:
            for kind, v, c in moves(alg):
                k = fingerprint.canonicalKey(c)
                if k is None: nxt.append(c); continue
                if k in seen: continue
                seen[k] = d; nxt.append(c)
                if k in tb and (best is None or d + tb[k] < best[0]): best = (d + tb[k], d, tb[k])
        levels.append(len(nxt)); frontier = nxt
        print('level', d, 'size', len(nxt), '%.0fs' % (time.time() - t1), 'best', best, flush=True)
    meta = dict(pk=PK, cls=CLS, idx=IDX, D=D, levels=levels, best=best, seenN=len(seen), TD=TD, key=h(k0))
    pickle.dump(dict(meta, frontier=frontier, seen=set(seen)), open('/tmp/tsm/fr_%s.pkl' % TAG, 'wb'))
    pickle.dump(meta, open('/tmp/tsm/fr_%s.meta.pkl' % TAG, 'wb'))
    print('FRONT', TAG, meta, 'frontier', len(frontier), '%.0fs' % (time.time() - t1), flush=True)
if MODE == 'shard':
    tb = loadball(); seen = fr['seen']; sub = fr['frontier'][LO:HI]; keys = set(); hits = {}; nodes = 0; t1 = time.time(); D = fr['D']
    import signal
    NODELIM = int(os.environ.get('NODELIM', '0')); skipped = []; nonekey = 0
    class TO(Exception): pass
    def _h(s_, f_): raise TO()
    signal.signal(signal.SIGALRM, _h)
    for ai, alg in enumerate(sub):
        try:
            if NODELIM: signal.alarm(NODELIM)
            mv = moves(alg); signal.alarm(0)
        except TO: skipped.append(LO + ai); continue
        for kind, v, c in mv:
            k = fingerprint.canonicalKey(c); nodes += 1
            if k is None: nonekey += 1
            if k is None or k in seen: continue
            keys.add(h(k))
            if k in tb: hits[h(k)] = min(hits.get(h(k), 99), D + tb[k])
    pickle.dump(dict(keys=keys, hits=hits, nodes=nodes, n=len(sub), skipped=skipped, nonekey=nonekey), open('/tmp/tsm/sh_%s_%d.pkl' % (TAG, LO), 'wb'))
    print('SHARD', TAG, LO, HI, 'expanded', len(sub), 'children generated', nodes, 'new keys', len(keys), 'no-key children', nonekey, 'nodes skipped (per-node time limit)', skipped, 'hits', hits, 'stats', stats, '%.0fs' % (time.time() - t1), flush=True)
if MODE == 'merge':
    m = fr; keys = set(); hits = {}; ex = 0; skipped = []; nok = 0; sh = sorted(glob.glob('/tmp/tsm/sh_%s_*.pkl' % TAG))
    for f in sh:
        s = pickle.load(open(f, 'rb')); keys |= s['keys']; ex += s['n']; skipped += s.get('skipped', []); nok += s.get('nonekey', 0)
        for k, t in s['hits'].items(): hits[k] = min(hits.get(k, 99), t)
    print('MERGE', TAG, 'child key', m['key'], 'D', m['D'], 'TD', m['TD'], 'levels<D', m['levels'], 'ball nodes through D-1', m['seenN'], 'shards', len(sh),
          'level-D frontier expanded', ex, 'new level-D keys', len(keys), 'ball nodes through D', m['seenN'] + len(keys), 'earlier best', m['best'],
          'skipped nodes', skipped, 'children without canonical key', nok, 'HITS at level D', hits if hits else 'none')
