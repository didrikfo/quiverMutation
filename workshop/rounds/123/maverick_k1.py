import sys, collections
sys.path.insert(0, 'workshop/rounds/123'); sys.path.insert(0, 'workshop/rounds/117')
import experimentalist_classlevel as E
n, h, moves, kmin = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
cls, _ = E.keyclass(n, h); K0 = min(h, n-4-h)
sub = [m for m in cls if E.Kof(m, n) >= kmin]; S = set(sub)
par = {m: m for m in sub}
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for m in sub:
    for r in E.neighbours(n, m, moves):
        if r in S: par[f(m)] = f(r)
lab = {m: E.key(n-1, m[1:]) for m in sub}
comp = collections.defaultdict(list)
for m in sub: comp[f(m)].append(m)
for K in range(kmin, K0+1):
    hit = [c for c in comp.values() if any(E.Kof(m, n) == K for m in c)]
    mixed = [c for c in hit if len({lab[m] for m in c if E.Kof(m, n) == K}) > 1]
    print('n=%d h=%d %s, rows K>=%d: %d rows, %d orbits; K=%d: hit %d mixed %d' % (n, h, moves, kmin, len(sub), len(comp), K, len(hit), len(mixed)))
print('  whole restricted orbit (all K >= %d) mixed: %d of %d; image keys over K>=%d rows: %d' % (kmin, sum(1 for c in comp.values() if len({lab[m] for m in c}) > 1), len(comp), kmin, len(set(lab.values()))))
