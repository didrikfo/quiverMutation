"""Round 123 (experimentalist): H-026 at class level.  Repo root:
  timeout 10m .venv/bin/python workshop/rounds/123/experimentalist_classlevel.py plan N H [N H ...]
  timeout 10m .venv/bin/python workshop/rounds/123/experimentalist_classlevel.py run  N H [moves]
For the Coxeter-key class of the lone-3 row at (n, h) -- ALL rows at length n carrying that key, every K --
partition the class into orbits under the symmetric (undirected, union-find) closure of the move families
  moves = rules            lnaMoves.rewritesOf only (what E-307's forward components used, but symmetric and on the whole class)
          mut   (default)  rules + edgeMoves + doubleMutation    (mutation-class orbits, as freeMoves.derivedOrbits(free=False, edges, doubles))
          der              mut + free move (stripLengthTwo)      (derived orbits)
and ask: is the image key (Coxeter key of the vertex-1 deletion, the n-1 row m[1:]) constant on an orbit,
(a) at each fixed K (K = free head stretch), (b) over the whole orbit.  Also q_S (sig, |det S|, q_S histogram of S = C + C^T,
toolsmith_qsscan.inv) per orbit, to say what the q_S classes of the key class are.  Moves preserve the derived class, so a
move out of the key class is a bug and is asserted against.  Controls: shuffled image labels (same counts, same orbits), and a
planted invariant (the K0-split image key replaced by a function of the orbit) which must come out constant.
"""
import sys, time, collections, random, os
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '117'))
from quivermutation import nakayama, lnaMoves, edgeMoves, doubleMutation, freeMoves


def cartan(n, r):
    rel = [(i+1, i+1+x) for i, x in enumerate(r) if x]
    C = np.zeros((n, n))
    for s in range(1, n+1):
        stop = min((e for a, e in rel if a >= s), default=n+1)
        for t in range(s, min(n, stop-1)+1): C[t-1, s-1] = 1
    return C


def key(n, r):   # identical to rounds/117/experimentalist_kpower.key
    C = cartan(n, r); Phi = -np.linalg.solve(C.T, C); p = np.poly(Phi); q = np.rint(p)
    assert np.max(np.abs(p-q)) < 1e-6
    return tuple(int(v) for v in q)


def lone(n, h): return tuple([0]*h + [3] + [0]*(n-3-h))
def Kof(m, n): return next((i for i, x in enumerate(m) if x), n)

_N = None; _TARGET = None
def _init(n, target): global _N, _TARGET; _N, _TARGET = n, target
def _chunk(rows): return [r for r in rows if key(_N, r) == _TARGET]


def keyclass(n, h):
    target = key(n, lone(n, h)); allr = list(nakayama.allRelationLengths(n))
    chunks = [allr[i::64] for i in range(64)]
    with Pool(4, _init, (n, target)) as p: out = p.map(_chunk, chunks)
    return sorted(r for c in out for r in c), len(allr)


def neighbours(n, r, moves):
    out = list(lnaMoves.rewritesOf(n, list(r)))
    if moves in ('mut', 'der'):
        out += [x for x, _ in edgeMoves.rewritesOf(n, r)]
        out += [x for x, _ in doubleMutation.rewritesOf(n, r)]
    if moves == 'der':
        out.append(freeMoves.stripLengthTwo(r))
    return [tuple(x) for x in out]


def orbits(n, cls, moves):
    S = set(cls); parent = {m: m for m in cls}; edges = 0
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for m in cls:
        for r in neighbours(n, m, moves):
            assert r in S, ('move left the key class', m, r)
            edges += 1; a, b = find(m), find(r)
            if a != b: parent[a] = b
    comp = collections.defaultdict(list)
    for m in cls: comp[find(m)].append(m)
    return list(comp.values()), edges


def analyse(n, h, moves):
    t = time.time(); law = min(h, n-4-h)
    cls, total = keyclass(n, h); print('n=%d h=%d (class of lone h=%d, law K0=%d): %d of %d rows carry the key  %.0fs' % (n, h, h, law, len(cls), total, time.time()-t), flush=True)
    imgs = {}; lab = {}
    for m in cls:
        k = key(n-1, m[1:]); imgs.setdefault(k, len(imgs)); lab[m] = imgs[k]
    tabK = collections.defaultdict(collections.Counter)
    for m in cls: tabK[Kof(m, n)][lab[m]] += 1
    print(' image keys on the whole class by K (label: rows):'); [print('   K', K, dict(sorted(c.items()))) for K, c in sorted(tabK.items())]
    comps, edges = orbits(n, cls, moves)
    print(' moves=%s: %d undirected move edges, %d orbits, sizes %s  %.0fs' % (moves, edges, len(comps), sorted(map(len, comps), reverse=True)[:12], time.time()-t), flush=True)
    # (a) per K, (b) whole orbit
    for K in sorted(tabK):
        mixed = [c for c in comps if len({lab[m] for m in c if Kof(m, n) == K}) > 1]
        hit = sum(1 for c in comps if any(Kof(m, n) == K for m in c))
        line = '   K=%d: orbits containing K-rows %d, mixed (>= 2 image keys at this K) %d' % (K, hit, len(mixed))
        if mixed:
            c = max(mixed, key=len); byl = collections.defaultdict(list)
            for m in c:
                if Kof(m, n) == K: byl[lab[m]].append(m)
            line += '  e.g. orbit of %d rows: %s' % (len(c), {l: ''.join(map(str, v[0])) for l, v in byl.items()})
        print(line, flush=True)
    whole = sum(1 for c in comps if len({lab[m] for m in c}) > 1)
    print('   whole orbit (any K): %d of %d orbits carry >= 2 image keys' % (whole, len(comps)))
    # the split at K = law: first counterexample, as a chain of moves if it exists
    K0 = law; mixedK0 = [c for c in comps if len({lab[m] for m in c if Kof(m, n) == K0}) > 1]
    if mixedK0:
        c = min(mixedK0, key=len); a = next(m for m in sorted(c) if Kof(m, n) == K0)
        b = next(m for m in sorted(c) if Kof(m, n) == K0 and lab[m] != lab[a])
        path = shortest(n, set(c), a, b, moves)
        print('   FIRST COUNTEREXAMPLE at K=%d (smallest mixed orbit, %d rows): %s (img %d) ~ %s (img %d), move path length %d:' % (K0, len(c), ''.join(map(str, a)), lab[a], ''.join(map(str, b)), lab[b], len(path)-1))
        print('     ' + ' - '.join('%s[K%d,i%d]' % (''.join(map(str, m)), Kof(m, n), lab[m]) for m in path))
    # same question on reduced rows only (no two-arrow relation, freeMoves.isReduced): K is then not changed by a free move
    for K in sorted(tabK):
        mixed = [c for c in comps if len({lab[m] for m in c if Kof(m, n) == K and freeMoves.isReduced(m)}) > 1]
        nr = sum(1 for c in comps for m in c if Kof(m, n) == K and freeMoves.isReduced(m))
        if K >= K0 - 1 and nr:
            ex = ''
            if mixed:
                c = max(mixed, key=len); byl = {}
                for m in sorted(c):
                    if Kof(m, n) == K and freeMoves.isReduced(m): byl.setdefault(lab[m], ''.join(map(str, m)))
                ex = '  e.g. %s' % byl
                if K == K0:
                    a, b = [tuple(map(int, s)) for s in list(byl.values())[:2]]
                    path = shortest(n, set(c), a, b, moves)
                    ex += '\n     path (%d moves): ' % (len(path)-1) + ' - '.join('%s[K%d,i%d]' % (''.join(map(str, m)), Kof(m, n), lab[m]) for m in path)
            print('   reduced rows only, K=%d: %d rows, image keys %s, mixed orbits %d%s' % (K, nr, sorted({lab[m] for c in comps for m in c if Kof(m, n) == K and freeMoves.isReduced(m)}), len(mixed), ex))
    # q_S per orbit (one representative; plus 3 more members as a consistency check)
    from toolsmith_qsscan import inv as qinv
    qs = collections.defaultdict(list); bad = 0; rnd = random.Random(1)
    for c in comps:
        v = qinv(n, c[0]); qs[v].append(c)
        for m in rnd.sample(c, min(3, len(c)-1)) if len(c) > 1 else []:
            bad += qinv(n, m) != v
    print(' q_S classes of the key class: %d (q_S inconsistent inside an orbit: %d checks failed)' % (len(qs), bad))
    for v, cs in qs.items():
        rowsq = [m for c in cs for m in c]; kk = collections.Counter(lab[m] for m in rowsq if Kof(m, n) == K0)
        print('   sig %s |det S| %s hist %s: %d orbits, %d rows; image keys at K=%d: %s' % (v[0], v[1], v[2], len(cs), len(rowsq), K0, dict(kk)))
    # controls
    rnd = random.Random(7)
    for K in (K0, K0-1):
        rowsK = [(m, i) for i, c in enumerate(comps) for m in c if Kof(m, n) == K]
        if not rowsK: continue
        def nm(labels):
            d = collections.defaultdict(set)
            for (m, i), l in zip(rowsK, labels): d[i].add(l)
            return sum(1 for s in d.values() if len(s) > 1)
        L = [lab[m] for m, _ in rowsK]; nulls = []
        for _ in range(200): L2 = L[:]; rnd.shuffle(L2); nulls.append(nm(L2))
        planted = nm([i % 2 for _, i in rowsK])   # a function of the orbit: must give 0 mixed
        print(' control K=%d: real mixed %d; shuffled labels mixed mean %.1f (min %d, max %d of 200); planted orbit-function mixed %d' % (K, nm(L), sum(nulls)/200, min(nulls), max(nulls), planted))
    print('%.0fs' % (time.time()-t))


def shortest(n, S, a, b, moves):
    adj = collections.defaultdict(set)
    for m in S:
        for r in neighbours(n, m, moves): adj[m].add(r); adj[r].add(m)
    prev = {a: None}; fr = [a]
    while fr and b not in prev:
        nf = []
        for x in fr:
            for y in adj[x]:
                if y not in prev: prev[y] = x; nf.append(y)
        fr = nf
    p = [b]
    while prev[p[-1]] is not None: p.append(prev[p[-1]])
    return p[::-1]


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'plan':
        a = sys.argv[2:]
        for n, h in zip(a[::2], a[1::2]):
            n, h = int(n), int(h); t = time.time(); allr = list(nakayama.allRelationLengths(n)); te = time.time()-t
            smp = random.Random(0).sample(allr, 2000); t = time.time(); [key(n, r) for r in smp]; tk = (time.time()-t)/2000
            t = time.time(); [neighbours(n, r, 'mut') for r in smp[:300]]; tn = (time.time()-t)/300
            print('n=%d h=%d: %d rows (enumerate %.0fs); key %.0f us/row -> class scan about %.0fs on 4 cores; move neighbours %.1f ms/row'
                  ' (class of c rows: c * %.1f ms single core)' % (n, h, len(allr), te, tk*1e6, len(allr)*tk/4*1.3, tn*1e3, tn*1e3))
    else:
        analyse(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4] if len(sys.argv) > 4 else 'mut')
