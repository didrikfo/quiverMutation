"""Round 051 (toolsmith): resumable target ball (LNAs + duals of key class cls, moves F and R as in rounds/049/toolsmith_tiltpath.py) to depth D.
Round 049's one-shot builder (no checkpoint) takes > 10 min at depth 6, so this saves state every DEADLINE seconds (exit 3) and resumes on rerun.
Usage (repo root): DEADLINE=500 timeout 10m .venv/bin/python -u workshop/rounds/051/toolsmith_tball.py cls D    -> /tmp/tsm/tball_c<cls>_d<D>.pkl = (seen{key:depth}, levels, False)"""
import sys, os, time, pickle
ARGV = sys.argv; CLS, D = int(ARGV[1]), int(ARGV[2])
sys.argv = ['x', '/tmp/tsm/none.pkl', str(CLS), '5', '6', '400000', 'none']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("def ball(")[0]
exec(compile(src, 'tp', 'exec'))
OUT = '/tmp/tsm/tball_c%d_d%d.pkl' % (CLS, D); ST = OUT + '.state'; DL = float(os.environ.get('DEADLINE', '500')); t0 = time.time()
if os.path.exists(ST): st = pickle.load(open(ST, 'rb'))
else:
    seen = {}; fr = []
    for s in classes[base]:
        k = fingerprint.canonicalKey(s)
        if k in seen: continue
        seen[k] = 0; fr.append(s)
    st = dict(seen=seen, d=1, done=0, cur=[], frontier=fr, levels=[len(fr)])
seen, d, frontier, cur, levels = st['seen'], st['d'], st['frontier'], st['cur'], st['levels']; done = st['done']
while d <= D:
    while done < len(frontier):
        alg = frontier[done]
        for kind, v, c in moves(alg):
            k = fingerprint.canonicalKey(c)
            if k is None or k in seen: continue
            seen[k] = d; cur.append(c)
        done += 1
        if time.time() - t0 > DL:
            pickle.dump(dict(seen=seen, d=d, done=done, cur=cur, frontier=frontier, levels=levels), open(ST, 'wb'))
            print('CHECKPOINT depth', d, 'expanded', done, 'of', len(frontier), 'keys', len(seen), flush=True); sys.exit(3)
    levels.append(len(cur)); print('level', d, len(cur), '%.0fs' % (time.time() - t0), flush=True)
    frontier, cur, done = cur, [], 0; d += 1
pickle.dump((seen, levels, False), open(OUT, 'wb')); print('TBALL', CLS, D, 'levels', levels, 'keys', len(seen), flush=True)
