"""Round 050 (toolsmith): why does c1 child 12 (parallel arrows) not finish level 4?  Times moves() per level-3 node, per-node alarm LIM seconds (env, default 20),
prints the distribution and the slowest nodes (arrow count, kind F or R).  Usage (repo root): toolsmith_slow12.py in.pkl idx"""
import sys, os, time, pickle, signal
A = sys.argv; PK, IDX = A[1], int(A[2])
sys.argv = ['x', PK, '1', '5', '6', '400000', 'none']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src, 'tp', 'exec'))
LIM = int(os.environ.get('LIM', '20'))
class TO(Exception): pass
def h_(s, f): raise TO()
signal.signal(signal.SIGALRM, h_)
recs = [r for r in pickle.load(open(PK, 'rb')) if r['kind'] == 'fail']
start = recs[IDX]['childObj']; seen = {fingerprint.canonicalKey(start)}; frontier = [start]
for d in range(1, 4):
    nxt = []
    for alg in frontier:
        for kind, v, c in moves(alg):
            k = fingerprint.canonicalKey(c)
            if k in seen: continue
            seen.add(k); nxt.append(c)
    frontier = nxt
print('level 3 frontier', len(frontier), flush=True)
times = []
for i, alg in enumerate(frontier):
    t1 = time.time(); signal.alarm(LIM)
    try: n = len(moves(alg)); st = 'ok %d moves' % n
    except TO: st = 'TIMEOUT'
    signal.alarm(0); times.append((time.time() - t1, i, st, alg.quiver.number_of_edges()))
times.sort(reverse=True)
print('total %.0fs' % sum(t for t, *_ in times), 'slowest:', [(round(t, 1), i, s, e) for t, i, s, e in times[:8]], 'median %.2f' % times[len(times) // 2][0], flush=True)
# second part: children with canonicalKey None (not deduplicated by the 049 ball, cannot meet anything)
fr2 = []; seen = {fingerprint.canonicalKey(start)}; frontier = [start]; nn = 0
for d in range(1, 4):
    nxt = []
    for alg in frontier:
        for kind, v, c in moves(alg):
            k = fingerprint.canonicalKey(c)
            if k is None: nn += 1; nxt.append(c); continue
            if k in seen: continue
            seen.add(k); nxt.append(c)
    print('depth', d, 'size', len(nxt), 'cumulative None-key children', nn, flush=True); frontier = nxt
nonek = [a for a in frontier if fingerprint.canonicalKey(a) is None]
print('level-3 nodes with key None:', len(nonek), [a.quiver.number_of_edges() for a in nonek])
for a in nonek[:4]:
    t1 = time.time(); signal.alarm(LIM)
    try: ms = moves(a); st = 'ok %d moves (None-key children %d)' % (len(ms), sum(fingerprint.canonicalKey(c) is None for _, _, c in ms))
    except TO: st = 'TIMEOUT'
    signal.alarm(0); print('None-key node moves:', st, '%.1fs' % (time.time() - t1), flush=True)
