"""Round 114 (toolsmith): profile one level of search at step 8. Run: .venv/bin/python workshop/rounds/114/toolsmith_s8prof.py"""
import sys, cProfile, pstats
sys.argv = ['x', 'list']
src = open('workshop/rounds/113/toolsmith_c1_chain.py').read().split("t0 = time.time()\nnexp, Z")[0]
exec(compile(src, 'c1', 'exec'))
nexp, Z = load_steps(); A, B = Z[8]['alg'], Z[8]['ch']
def lvl(x):
    t = time.time(); m = moves(x); t1 = time.time(); ks = [key(c) for _, _, c in m]; return m, t1 - t, time.time() - t1
for nm, x in (('A', A), ('B', B)):
    m, a, b = lvl(x); print(nm, 'moves', [(k, v) for k, v, _ in m], 'moves time %.1f key time %.1f' % (a, b), flush=True)
    for k, v, c in m:
        m2, a2, b2 = lvl(c); print('   child', k, v, 'parallel', parallel(c), 'n moves', len(m2), 'times %.1f %.1f' % (a2, b2), flush=True)
cProfile.run('moves(A)', '/tmp/w/prof'); pstats.Stats('/tmp/w/prof').sort_stats('cumtime').print_stats(14)
