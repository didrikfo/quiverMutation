"""Round 126 (skeptic), response to referee: the same test as skeptic_inj.py on two variants of the N generator.
K: every binomial p - q of skeptic_conv3.gen('N') becomes p + c q with c drawn from +-{2..9} (rng seeded seed + 10**6).
T: as N, then ONE binomial that has a third parallel path of the same length gets that path with coefficient c in {2..9};
   parents with no such binomial are skipped (so 'parents' counts only 3-term parents).
Pre-registered: HIT = injectivity (LIB and OWN) == Hom on every (vertex, side) case.
usage: skeptic_inj2.py KIND N SEED0 COUNT   (KIND in K, T, N; from the repo root). N = plain generator (for n = 8)."""
import sys, random
from collections import Counter
sys.path.insert(0, 'workshop/rounds/126')
import skeptic_inj as si
c3, P = si.c3, si.P

def variant(kind, n, seed):
    A = c3.gen('N', n, random.Random(seed))
    if A is None or kind == 'N': return A
    rng = random.Random(seed + 10 ** 6)
    arrows = [(a, b) for a, b, _ in A.arrows]
    rels = [dict(r) for r in A.rels]
    if kind == 'K':
        for r in rels:
            if len(r) == 2:
                p, q = sorted(r)
                r[p] = 1; r[q] = (rng.choice([-1, 1]) * rng.randint(2, 9)) % P
    else:
        def ends(p): return (A.arrows[p[0]][0], A.arrows[p[-1]][1], len(p))
        par = {}
        for u in range(n):
            for w in range(n):
                if u != w:
                    for p in A.allpaths(u, w): par.setdefault(ends(p), []).append(p)
        opts = [(k, x) for k, r in enumerate(rels) if len(r) == 2
                for x in par.get(ends(next(iter(r))), []) if x not in r]
        if not opts: return None
        k, x = rng.choice(opts); rels[k][x] = rng.randint(2, 9)
    return c3.AlgG(n, arrows, rels, True)

def sweep(kind, n, s0, cnt):
    tab = Counter(); dis = []
    for seed in range(s0, s0 + cnt):
        A = variant(kind, n, seed)
        if A is None: continue
        tab['parents'] += 1
        for i, side, gate, lib, own, hom in si.cases(A):
            tab['gate %s Hom %s LIB %s OWN %s' % (gate, hom, lib, not own)] += 1
            if lib != hom or (not own) != hom: dis.append((seed, i + 1, side, gate, hom, lib, not own))
            if gate and not hom:
                sh = si.e080(A, i, side, own); tab['gate over-admits; shape %s' % sh.split(':')[0]] += 1
                if sh.startswith('other'): print('  other: seed', seed, 'v%d' % (i + 1), side, sh)
    for k in sorted(tab): print(tab[k], k)
    print('injectivity != Hom:', len(dis), dis[:10])

if __name__ == '__main__':
    sweep(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
