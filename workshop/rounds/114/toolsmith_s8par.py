"""Round 114 (toolsmith), T10 step 8: parallel (4 cores), resumable, balanced bidirectional chain search for E-201 step 8 (n=7 class 1).
Reuses 113/toolsmith_c1_chain.py (moves = SEARCH filter, verify = Hom/Cartan/K0/loop/meet-iso verifier).
  .venv/bin/python -u workshop/rounds/114/toolsmith_s8par.py [STEP=8] [DEADLINE_S=480]   (rerun while it prints PAUSED; ckpt /tmp/w/s8_STEP.pkl)
Expands always the side with the smaller whole ball (level by level, each level finished before the next). Needs /tmp/skeptic_c1_13_ck1.pkl (069/skeptic_c1_13.py 1 540 23000, twice)."""
import sys, os, time, pickle
STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 8; DL = float(sys.argv[2]) if len(sys.argv) > 2 else 480
sys.argv = ['x', 'list']
src = open('workshop/rounds/113/toolsmith_c1_chain.py').read().split("t0 = time.time()\nnexp, Z")[0]
exec(compile(src, 'c1', 'exec'))
import multiprocessing as mp
def expand(alg):
    try:
        return [(kind, v, c, key(c)) for kind, v, c in moves(alg)], key(alg)
    except Exception as e:
        return [], key(alg)
if __name__ == '__main__':
    T0 = time.time(); nexp, Z = load_steps(); A, B, v0 = Z[STEP]['alg'], Z[STEP]['ch'], Z[STEP]['v']
    ck = '/tmp/w/s8_%d.pkl' % STEP; os.makedirs('/tmp/w', exist_ok=True)
    if os.path.exists(ck): st = pickle.load(open(ck, 'rb')); print('resumed', [len(s) for s in st['S']], 'levels', st['lv'], flush=True)
    else: st = dict(S=[{key(A): (None, None, None, A)}, {key(B): (None, None, None, B)}], fr=[[A], [B]], pos=[0, 0], nxt=[[], []], lv=[0, 0], s=None)
    S = st['S']; pool = mp.Pool(4); found = None
    while found is None:
        if st['s'] is None or st['pos'][st['s']] == 0 and not st['nxt'][st['s']]:
            pass
        if st['pos'][0] == 0 and st['pos'][1] == 0 or st['s'] is None:   # choose a side at a level boundary
            st['s'] = 0 if len(S[0]) <= len(S[1]) else 1
        s = st['s']; fr = st['fr'][s]
        while st['pos'][s] < len(fr) and found is None:
            if time.time() - T0 > DL:
                pickle.dump(st, open(ck, 'wb')); print('PAUSED', [len(x) for x in S], 'levels', st['lv'], '%.0fs' % (time.time() - T0)); sys.exit()
            batch = fr[st['pos'][s]: st['pos'][s] + 16]; st['pos'][s] += len(batch)
            for res, ka in pool.map(expand, batch, chunksize=1):
                for kind, v, c, kc in res:
                    if kc is None or kc in S[s]: continue
                    S[s][kc] = (ka, kind, v, c); st['nxt'][s].append(c)
                    if kc in S[1 - s] and found is None: found = (s, kc)
        if found: break
        st['fr'][s] = st['nxt'][s]; st['nxt'][s] = []; st['pos'][s] = 0; st['lv'][s] += 1; st['s'] = None
        print('level done side', s, 'levels', st['lv'], 'balls', [len(x) for x in S], '%.0fs' % (time.time() - T0), flush=True)
        if not st['fr'][s]: print('EXHAUSTED side', s); sys.exit()
    mvA, mvB = trace(S, found[0], found[1], A, B)
    ok, la, lb, meet, par = verify(A, B, mvA, mvB)
    print('A-side', ' '.join(la), '| B-side', ' '.join(lb), '| edges', len(mvA) + len(mvB), 'meet-isomorphic', meet, 'parallel-at-meet', par)
    print('STEP', STEP, 'CHAIN VERDICT', 'YES' if ok else 'NO (verifier rejected)', 'fails', sorted(VERIFY_FAILS), '%.0fs' % (time.time() - T0))
