"""Round 045 response: all 10 pairs, library search.meetingPoints(A,B,3,alsoDual=True) vs author's reach(...,'guard').
Usage: .venv/bin/python workshop/rounds/045/experimentalist_crosscheck.py"""
import sys, time
sys.argv = ['x']
src = open('workshop/rounds/045/experimentalist_t10b.py').read().split("\nP = pairs()")[0]
exec(compile(src, 't10b', 'exec'))
P = pairs()
ok = 0
for i, (a, b) in enumerate(P):
    t0 = time.time()
    A_ = nk.LinearNakayamaAlgebra(8, list(a)); B_ = nk.LinearNakayamaAlgebra(8, list(b))
    lib = search.meetingPoints(A_, B_, 3, alsoDual=True)
    libk = {m[0]: (len(m[1]), len(m[2])) for m in lib}
    st = Counter(); res = []
    for alg in (A_, B_):
        f = {}
        for dual in (False, True):
            g, _ = reach(alg, 3, 'guard', st, dual)
            for k, d in g.items():
                if k not in f or d < f[k]: f[k] = d
        res.append(f)
    sh = {k: (res[0][k], res[1][k]) for k in set(res[0]) & set(res[1])}
    same = set(sh) == set(libk)
    ok += same
    print('pair %d %s->%s: library %d meeting(s) split %s | author %d split %s | same keys %s | %.0fs' % (
        i, ''.join(map(str, a)), ''.join(map(str, b)), len(lib), sorted(libk.values()), len(sh), sorted(sh.values()), same, time.time() - t0), flush=True)
print('same meeting-key set on %d of %d pairs' % (ok, len(P)))
