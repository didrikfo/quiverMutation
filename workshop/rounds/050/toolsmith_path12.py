"""Round 050 (toolsmith): printed + replayed path for c1 child 12 (its ball contains nodes without a canonical key, so 049's reachPath BFS does not finish).
Here: BFS to depth 3 with parent pointers (nodes without key kept, not deduplicated, as in the ball), then expand each depth-3 node (per-node alarm NODELIM s) and
take the first child whose (env DEP = depth of the node level expanded, default 3) canonical key lies in the target ball (parent-pointer cache from toolsmith_paths.py buildball).  Prints the child -> meet and LNA/dual -> meet
move lists and replays both with fresh move generation.   Usage (repo root): toolsmith_path12.py in.pkl 1 idx"""
import sys, os, time, pickle, hashlib, signal
A = sys.argv; PK, CLS, IDX = A[1], A[2], int(A[3])
sys.argv = ['x', PK, CLS, '7', 'paths']
src = open('workshop/rounds/049/toolsmith_paths.py').read().split("recs = [r for r in pickle.load")[0]
exec(compile(src.replace("sys.exit()", "pass"), 'pp', 'exec'))
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']; start = recs[IDX]['childObj']
class TO(Exception): pass
def _h(s_, f_): raise TO()
signal.signal(signal.SIGALRM, _h); NODELIM = int(os.environ.get('NODELIM', '150'))
k0 = fingerprint.canonicalKey(start); seen = {k0}; nodes = [(start, None, None, None)]; level = [0]; found = []
DEP = int(os.environ.get('DEP', '3'))
for d in range(1, DEP + 1):
    nl = []
    for i in level:
        for kind, v, c in moves(nodes[i][0]):
            k = fingerprint.canonicalKey(c)
            if k is not None:
                if k in seen: continue
                seen.add(k)
            nodes.append((c, i, kind, v)); nl.append(len(nodes) - 1)
    level = nl
print('depth-%d nodes' % DEP, len(level), flush=True)
t0 = time.time()
for i in level:
    signal.alarm(NODELIM)
    try: mv = moves(nodes[i][0]); signal.alarm(0)
    except TO: print('node', i, 'skipped (time limit)'); continue
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
print('no hit')
