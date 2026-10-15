"""Follow-up: classes whose d>=2 free-stretch cells split; list (row, i, image) to see what separates them. ~5 s"""
import runpy, sys, io, contextlib, collections
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path('workshop/rounds/067/experimentalist_deltable.py')
rows, img, dist, R, C, L, N = (g[k] for k in 'rows img dist R C L N'.split())
for c in ['P^(4)_(1,3)']:
    for d in (2,):
        out = collections.defaultdict(list)
        for rs in rows:
            if C[rs] != c: continue
            for i in range(1, N+1):
                if dist(R[rs], i) >= d:
                    nz = [(v, r) for v, r in enumerate(R[rs], 1) if r]
                    left = i < nz[0][0]; right = i > nz[-1][0] + nz[-1][1]
                    out[(L[img[rs, i][0]], 'left' if left else 'right' if right else 'gap', dist(R[rs], i))].append(rs)
        for k, v in sorted(out.items()): print(c, k, len(v), v[:3])
# free END only (left of all or right of all), K = number of free vertices at that end, vs image class purity
print('free-end strip: K = free vertices at the end containing i; classes, purity of image class')
res = collections.defaultdict(lambda: collections.defaultdict(set))
for rs in rows:
    nz = [(v, r) for v, r in enumerate(R[rs], 1) if r]
    if not nz: continue
    K0 = nz[0][0] - 1; K1 = N - (nz[-1][0] + nz[-1][1])
    for i in range(1, N+1):
        if i <= K0: res[C[rs]][('L', K0)].add(img[rs, i][0])
        if i > N - K1: res[C[rs]][('R', K1)].add(img[rs, i][0])
for K in (1, 2, 3, 4):
    cl = {c for c, d in res.items() if any(k[1] >= K for k in d)}
    pure = [c for c in cl if len({x for k, s in res[c].items() if k[1] >= K for x in s}) == 1]
    print('K>=%d: classes with such an end %d, one image class %d' % (K, len(cl), len(pure)))
