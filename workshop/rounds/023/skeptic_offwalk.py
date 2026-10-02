"""Round 023 (skeptic, T5): hunt for a counterexample to "reject iff long square" OFF the guarded walks.
Random hand-built parents (acyclic quiver on n vertices, no parallel arrows, 0..4 random homogeneous relations, mode B
seeds a long-sided square on purpose). At every vertex v with an out-arrow: gate (mutationIsPossibleAtVertex), tiltingPlus,
hasLongSquare(alg.rels) [copied from rounds/022] and on procedure.relationsFrom(alg), and minimality of the long relation
(is the long relation's truncation a_x..v already in ideal?). Classify gate-admitted steps:
  A  = long square and tiltingPlus True   (counterexample to 'long => reject')
  B  = no long square and tiltingPlus False (counterexample to 'reject => long')
Usage: skeptic_offwalk.py n samples seed mode(A|B|C)"""
import sys, random, time, itertools
from collections import Counter
sys.path.insert(0, '.'); _argv = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/001/scholar_h015.py').read().replace("\nmain()\n", "\n")
exec(compile(src, 'h015', 'exec'))
sys.argv = _argv
n, N, seed, mode = int(_argv[1]), int(_argv[2]), int(_argv[3]), _argv[4]
rnd = random.Random(seed)
def vpaths(adj, s, t, minlen=1):
    out = []
    def go(c, p):
        if c == t and len(p) > 1: out.append(tuple(p)); return
        for d in adj.get(c, ()): go(d, p + [d])
    go(s, [s]); return [p for p in out if len(p) - 1 >= minlen]
def hasLong(rels, quiver, v):
    outs = ap.arrowsOutOf(quiver, v)
    if len(outs) != 1: return False
    e = outs[0][1]
    for rel in rels:
        rel = [list(q) for q in rel]
        if len(rel) >= 2 and all(len(q) >= 3 and q[-1] == e and q[-2] == v for q in rel) and len({q[0] for q in rel}) == 1 and len({q[-3] for q in rel}) == len(rel):
            return True
    return False
def gen():
    while True:
        edges = set()
        for i in range(1, n):  # connected-ish: each vertex i>1 gets a predecessor
            pass
        p = rnd.uniform(0.25, 0.6)
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                if rnd.random() < p and j - i <= 3: edges.add((i, j))
        adj = {}
        for (i, j) in edges: adj.setdefault(i, []).append(j)
        rels = []
        pairs = [(s, t) for s in range(1, n + 1) for t in range(s + 1, n + 1)]
        # seeded long square (mode A,B): same start a, two distinct paths of >=3 arrows into v then a single arrow v->e
        if mode in 'AB':
            cand = []
            for s, t in pairs:
                ps = vpaths(adj, s, t)
                for L in {len(x) for x in ps}:
                    same = [x for x in ps if len(x) == L and L - 1 >= 4 - 1 + 0]
                cand.append((s, t, ps))
            rnd.shuffle(cand)
            for s, t, ps in cand:
                byL = {}
                for x in ps: byL.setdefault(len(x), []).append(x)
                for L, xs in byL.items():
                    if L >= 4 and len(xs) >= 2:
                        k = rnd.choice([2, 2, 3]) if len(xs) >= 3 else 2
                        w = rnd.choice(xs)[-2]; xs2 = [x for x in xs if x[-2] == w]
                        if len(xs2) < 2: continue
                        sel = rnd.sample(xs2, min(k, len(xs2)))
                        if len({x[-3] for x in sel}) == len(sel): rels.append([list(x) for x in sel]); break
                if rels: break
            if not rels: continue
        for _ in range(rnd.choice([0, 0, 1, 1, 2, 3])):
            s, t = rnd.choice(pairs); ps = vpaths(adj, s, t)
            byL = {}
            for x in ps: byL.setdefault(len(x), []).append(x)
            opts = [(L, xs) for L, xs in byL.items() if L >= 3]
            if not opts: continue
            L, xs = rnd.choice(opts)
            if len(xs) >= 2 and rnd.random() < 0.7: rels.append([list(x) for x in rnd.sample(xs, min(len(xs), rnd.choice([2, 2, 3])))])
            else: rels.append([list(rnd.choice(xs))])
        return edges, rels
def longDict(quiver, rels, v):
    """(shaped, genuine) on relationsFrom-style dict relations: shaped = E-097 shape; genuine = truncated combination nonzero mod I."""
    outs = ap.arrowsOutOf(quiver, v)
    if len(outs) != 1: return (False, False)
    e = outs[0]; shaped = False; genuine = False
    for r in rels:
        ps = list(r)
        if len(ps) >= 2 and all(len(q) >= 3 and q[-1][0] == v and q[-1][1] == e[1] and q[-2][1] == v for q in ps) and len({q[0][0] for q in ps}) == 1 and len({q[-2][0] for q in ps}) == len(ps):
            shaped = True
            a = ps[0][0][0]
            trunc = ap.combination([(q[:-1], r[q]) for q in ps])
            res = ap.reduceAgainstPivots(trunc, ap.idealBasis(quiver, rels, a, v))
            if res: genuine = True
    return (shaped, genuine)
def perArrow(quiver, rels, v):
    """for each out-arrow of v: does a shaped relation (E-097 shape, ending in that arrow) exist?"""
    res = []
    for e in ap.arrowsOutOf(quiver, v):
        ok = False
        for r in rels:
            ps = list(r)
            if len(ps) >= 2 and all(len(q) >= 3 and q[-1][0] == v and q[-1][1] == e[1] and q[-2][1] == v for q in ps) and len({q[0][0] for q in ps}) == 1 and len({q[-2][0] for q in ps}) == len(ps): ok = True
        res.append(ok)
    return tuple(res)
cnt = Counter(); ex = {'A': [], 'B': []}; t0 = time.time()
for it in range(N):
    edges, rels0 = gen()
    A = pathAlgebra.PathAlgebra(); A.add_vertices_from(range(1, n + 1))
    for (i, j) in sorted(edges): A.add_arrow(i, j)
    try:
        for r in rels0: A.add_rel(r)
        if any(ap.isIllegalRelation(A.quiver, [tuple(map(lambda a: a, q)) for q in r]) for r in []): continue
        rels = procedure.relationsFrom(A)
    except Exception as ex_: cnt['err'] += 1; continue
    cnt['parents'] += 1
    for v in sorted(A.vertices()):
        if not ap.arrowsOutOf(A.quiver, v): continue
        try:
            g = bool(mutation.mutationIsPossibleAtVertex(A, v)); t = bool(tiltingPlus(A.quiver, rels, v))
        except Exception as e: cnt['err2'] += 1; continue
        L1 = hasLong(A.rels, A.quiver, v); L2s, L2 = longDict(A.quiver, rels, v)
        key = ('gate' if g else 'refused', 'tilt' if t else 'rej', 'Lrels' if L1 else '-', 'shaped' if L2s else '-', 'genuine' if L2 else '-')
        cnt[key] += 1
        if g and (L2 == t):
            kind = 'A' if t else 'B'
            if len(ex[kind]) < 400: ex[kind].append((sorted(edges), [list(map(list, r)) for r in A.rels], v, key, 'outs', len(ap.arrowsOutOf(A.quiver, v)), 'perarrow', perArrow(A.quiver, rels, v)))
print('n', n, 'samples', N, 'seed', seed, 'mode', mode, '%.0fs' % (time.time() - t0))
for k, c in sorted(cnt.items(), key=str): print(k, c)
for kind in 'AB':
    print('examples', kind, len(ex[kind]))
    for e in ex[kind]: print('  ', e)
print('B by (outs, perarrow):', Counter((e[5], e[7]) for e in ex['B']))
