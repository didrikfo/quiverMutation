"""Round 051 (toolsmith): printed + replayed total-13 path for c1 child 14 (child depth 7 + target depth 6).  Three stages, each < 10 min:
  tpar   : target ball depth 6 with parent pointers -> /tmp/tsm/tpar6_c<cls>.pkl   (resumable, env DEADLINE)
  cball  : child ball to depth 6 with parent pointers (nodes without key kept, not deduplicated) -> /tmp/tsm/cb_<idx>.pkl
  find   : expand the level-6 nodes in order (per-node alarm NODELIM, env LO:HI slice SL), take the first child whose key lies in the target ball, print child->meet
           and LNA/dual->meet move lists, replay both with fresh move generation.
Usage (repo root): toolsmith_path13.py in.pkl cls idx stage.   env CAP = canonicalKey cap (default 720)."""
import sys, os, time, pickle, hashlib, signal
A = sys.argv; PK, CLS, IDX, STAGE = A[1], int(A[2]), int(A[3]), A[4]
sys.argv = ['x', PK, str(CLS), '5', '6', '400000', 'none']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
CAP = int(os.environ.get('CAP', '720'))
if CAP != 720:
    _ck = fingerprint.canonicalKey; fingerprint.canonicalKey = lambda a, cap=CAP, gauge=True: _ck(a, cap=cap, gauge=gauge)
TP = '/tmp/tsm/tpar6_c%d.pkl' % CLS; DL = float(os.environ.get('DEADLINE', '500')); t0 = time.time()
if STAGE == 'tpar':
    ST = TP + '.state'
    if os.path.exists(ST): st = pickle.load(open(ST, 'rb'))
    else:
        seen = {}; fr = []
        for s in classes[base]:
            k = fingerprint.canonicalKey(s)
            if k in seen: continue
            seen[k] = (0, None, None, None); fr.append((k, s))
        st = dict(seen=seen, d=1, done=0, cur=[], frontier=fr)
    seen, d, frontier, cur, done = st['seen'], st['d'], st['frontier'], st['cur'], st['done']
    while d <= 6:
        while done < len(frontier):
            pk_, alg = frontier[done]
            for kind, v, c in moves(alg):
                k = fingerprint.canonicalKey(c)
                if k is None or k in seen: continue
                seen[k] = (d, pk_, kind, v); cur.append((k, c))
            done += 1
            if time.time() - t0 > DL:
                pickle.dump(dict(seen=seen, d=d, done=done, cur=cur, frontier=frontier), open(ST, 'wb')); print('CHECKPOINT', d, done, len(frontier), flush=True); sys.exit(3)
        print('level', d, len(cur), flush=True); frontier, cur, done = cur, [], 0; d += 1
    pickle.dump(seen, open(TP, 'wb')); print('TPAR', len(seen)); sys.exit()
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']; start = recs[IDX]['childObj']
CB = '/tmp/tsm/cb_%d.pkl' % IDX
if STAGE == 'cball':
    k0 = fingerprint.canonicalKey(start); seen = {k0}; nodes = [(start, None, None, None)]; level = [0]
    for d in range(1, 7):
        nl = []
        for i in level:
            for kind, v, c in moves(nodes[i][0]):
                k = fingerprint.canonicalKey(c)
                if k is not None:
                    if k in seen: continue
                    seen.add(k)
                nodes.append((c, i, kind, v)); nl.append(len(nodes) - 1)
        level = nl; print('level', d, len(level), '%.0fs' % (time.time() - t0), flush=True)
    pickle.dump((nodes, level), open(CB, 'wb')); sys.exit()
# find
tpar = pickle.load(open(TP, 'rb')); tball = {k: v[0] for k, v in tpar.items()}
nodes, level = pickle.load(open(CB, 'rb'))
seedkeys = {}
for i, a_ in enumerate(classes[base]): seedkeys.setdefault(fingerprint.canonicalKey(a_), i)
def tag(alg): return hashlib.md5(repr(fingerprint.canonicalKey(alg)).encode()).hexdigest()[:6]
def follow(st, seq):
    cur = st
    for kind, v, kx in seq:
        cand = [c for kd, vv, c in moves(cur) if kd == kind and vv == v and fingerprint.canonicalKey(c) == kx]
        if not cand: return None
        cur = cand[0]
    return fingerprint.canonicalKey(cur)
class TO(Exception): pass
def _h(s_, f_): raise TO()
signal.signal(signal.SIGALRM, _h); NODELIM = int(os.environ.get('NODELIM', '100'))
LO, HI = (int(x) for x in os.environ.get('SL', '0:%d' % len(level)).split(':'))
for n_, i in enumerate(level[LO:HI]):
    signal.alarm(NODELIM)
    try: mv = moves(nodes[i][0]); signal.alarm(0)
    except TO: print('node', LO + n_, 'skipped'); continue
    for kind, v, c in mv:
        k = fingerprint.canonicalKey(c)
        if k is not None and k in tball:
            cs = [(kind, v, k)]; j = i
            while nodes[j][1] is not None: cs.append((nodes[j][2], nodes[j][3], fingerprint.canonicalKey(nodes[j][0]))); j = nodes[j][1]
            cs.reverse(); ts = []; kk = k
            while tpar[kk][1] is not None: d_, pk_, kd, vv = tpar[kk]; ts.append((kd, vv, kk)); kk = pk_
            ts.reverse(); seed = kk
            ok = follow(start, cs) == k and follow(classes[base][seedkeys[seed]], ts) == k
            f = lambda sq: ' '.join('%s%d' % (kd, v) for kd, v, _ in sq) or '-'
            print('child', IDX, 'key', tag(start), 'replay', 'ok' if ok else 'REPLAY-FAILED', 'total', len(cs) + len(ts), '| child->meet:', f(cs), '| LNA/dual #%d ->meet:' % seedkeys[seed], f(ts), '%.0fs' % (time.time() - t0)); sys.exit()
print('no hit in slice')
