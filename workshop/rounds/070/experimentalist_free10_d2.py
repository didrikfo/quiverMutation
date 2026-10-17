"""The two n = 10 classes where D >= 2 fails: split of image class by side and by the length of the nearest relation.
usage: .venv/bin/python workshop/rounds/070/experimentalist_free10_d2.py  (needs experimentalist_free10.py's definitions; ~15 s)"""
import sys, collections, io, contextlib
sys.argv = ['x']
with contextlib.redirect_stdout(io.StringIO()):
    import runpy; g = runpy.run_path('workshop/rounds/070/experimentalist_free10.py')
members, img, dist, kind, N = (g[k] for k in 'members img dist kind N'.split())
for c, mem in members.items():
    S = collections.defaultdict(collections.Counter)
    for m in mem:
        for i in range(1, N + 1):
            if dist(m, i) == 2:
                near = [(r, v) for v, r in enumerate(m, 1) if r and min(abs(i - v), abs(i - v - r)) == 2]
                S[(kind(m, i), tuple(sorted({r for r, v in near})))][img(m, i)] += 1
    allim = {t for h in S.values() for t in h}
    d2 = {t for m in mem for i in range(1, N + 1) if dist(m, i) >= 2 for t in [img(m, i)]}
    if len(d2) > 1:
        print(str(c)[-30:], len(mem))
        for k, h in sorted(S.items()): print('   ', k, dict(h))
