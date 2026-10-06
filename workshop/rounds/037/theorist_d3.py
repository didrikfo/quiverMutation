"""Round 037 (theorist): hand examples with d_i >= 3 and J_i != 0 (gate-admitted); LNA key test and bounded BFS for an LNA in the mutation class.
Usage: theorist_d3.py hand | bfs NAME maxexp"""
import sys
from collections import Counter
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
_a = ARGV; exec(compile(src, 'bd', 'exec')); _a = ARGV
def dims(alg, v):
    Q = alg.quiver; rels = procedure.relationsFrom(alg); res = {}
    for i in Q.nodes:
        if i == v: continue
        P = ap.allPathsBetween(Q, i, v)
        if not P: continue
        res[i] = (len(P) - len(ap.idealBasis(Q, rels, i, v)), perI(alg, v).get(i, 0))
    return res
a3 = [(1,2),(1,3),(1,4),(2,5),(3,5),(4,5),(5,6)]
F = {  # name: (arrows, rels, nverts, v)
 'T1 one comm 1-2-5-6=1-3-5-6': (a3, [[[1,2,5,6],[1,3,5,6]]], 6, 5),
 'T2 two comm (m=3, one out)':  (a3, [[[1,2,5,6],[1,3,5,6]],[[1,3,5,6],[1,4,5,6]]], 6, 5),
 'K1 enum hit (1->2,1-3-2,1-4-2; outs 5,6; all commute)': ([(1,2),(1,3),(3,2),(1,4),(4,2),(2,5),(2,6)], [[[1,2,5],[1,4,2,5]],[[1,2,6],[1,4,2,6]],[[1,3,2,5],[1,4,2,5]],[[1,3,2,6],[1,4,2,6]]], 6, 2),
 'T3 chain a-b: 2 comm + 1 mono? (3-5-6 zero)': (a3, [[[1,2,5,6],[1,3,5,6]],[[3,5,6]]], 6, 5),
}
def mk(name):
    arr, rl, n, v = F[name]; return build(arr, rl, list(range(1, n+1))), v
if _a[1] == 'hand':
    ks = lnaKeys(6)
    for name in F:
        A, v = mk(name); key = search._coxeterKeyOrNone(A)
        print(name, '| gate', mutation.mutationIsPossibleAtVertex(A, v), '| (d,J)', dims(A, v), '| key in LNA keys', key in ks, key)
elif _a[1] == 'bfs':
    name, maxexp = _a[2], int(_a[3])
    if name.startswith('file'):
        import json; h = json.load(open(name.split(':')[1]))[int(name.split(':')[2])]; A = build([tuple(x) for x in h[0]], h[1], list(range(1, 1 + len({x for e in h[0] for x in e})))); v = h[3]
    else: A, v = mk(name)
    ks = lnaKeys(len(A.vertices()))
    seen = set(); fr = [A]; n = 0; hits = []; tab = Counter()
    def isLNA(alg):
        return any(fingerprint.canonicalKey(alg) == fingerprint.canonicalKey(l) or fingerprint.canonicalKey(alg) == fingerprint.canonicalKey(pathAlgebra.dualPathAlgebra(l)) for l in LN)
    LN = list(nk.LinearNakayamaAlgebra.allOfLength(len(A.vertices())))
    LK = set()
    for l in LN:
        for x in (l, pathAlgebra.dualPathAlgebra(l)): LK.add(fingerprint.canonicalKey(x))
    while fr and n < maxexp:
        cur, fr[:] = fr[:], []
        for alg in cur:
            k = fingerprint.canonicalKey(alg)
            if k in seen: continue
            seen.add(k); n += 1
            if k in LK: hits.append(n)
            for w in sorted(alg.vertices()):
                if not mutation.mutationIsPossibleAtVertex(alg, w): continue
                for i, dj in dims(alg, w).items(): tab[dj] += 1
                raw = mutation.quiverMutationAtVertex(alg, w)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                fr.append(reduction.reducePathAlgebra(raw))
    print(name, 'expanded', n, 'LNA/dual-LNA hits at', hits[:5], '(d,J) tab', dict(sorted(tab.items())))
if _a[1] == 'enum':
    import itertools
    n = int(_a[2]); ks = lnaKeys(n); res = Counter(); ex = {}
    # cores: 3 parallel-ish paths i->v of lengths (2,2,2),(2,2,3),(2,3,3); out arrows 1 or 2; relations: commutations of paths then b
    cores = {}
    def core(lens, nout):
        arrows = []; paths = []; nxt = 3; v = 2
        for L in lens:
            nodes = [1] + list(range(nxt, nxt + L - 1)) + [v]; nxt += L - 1
            arrows += list(zip(nodes, nodes[1:])); paths.append(nodes)
        outs = list(range(nxt, nxt + nout)); arrows += [(v, o) for o in outs]
        return arrows, paths, v, outs
    for lens in ((2,2,2),(2,2,3),(2,3,3),(3,3,3),(1,2,2),(1,2,3)):
        for nout in (1, 2):
            arr, paths, v, outs = core(lens, nout)
            k0 = len({x for e in arr for x in e})
            if k0 > n: continue
            # relation sets: choose for each out b a partition of paths into commuting groups (as chain of pair relations) ; plus optional zeros
            pairs = [(a, b) for a in range(3) for b in range(a+1, 3)]
            relopts = []
            for mask in itertools.product(range(2), repeat=len(pairs) * nout):
                rl = []
                for j, (a, b) in enumerate(pairs):
                    for oi, o in enumerate(outs):
                        if mask[j * nout + oi]: rl.append([paths[a] + [o], paths[b] + [o]])
                if rl: relopts.append(rl)
            extra = n - k0; newv = list(range(k0 + 1, n + 1)); nodes0 = list(range(1, k0 + 1))
            def attach(i, nodes, arrows):
                if i == len(newv): yield arrows; return
                w = newv[i]
                for u in nodes:
                    for dirn in (0, 1): yield from attach(i + 1, nodes + [w], arrows + [(u, w) if dirn == 0 else (w, u)])
            for att in attach(0, nodes0, []):
                for rl in relopts:
                    for zmask in range(2 ** len(att) if extra else 1):
                        rels = list(rl); arrows = arr + att
                        for j, (x, y) in enumerate(att):
                            if zmask >> j & 1:
                                for (c, d) in arrows:
                                    if d == x and (c, d) != (x, y): rels.append([[c, x, y]])
                                    if c == y and (c, d) != (x, y): rels.append([[x, y, d]])
                        try:
                            A = build(arrows, rels, list(range(1, n + 1)))
                            if list(nx.simple_cycles(A.quiver)): continue
                            if any(ap.isIllegalRelation(A.quiver, rr) for rr in procedure.relationsFrom(A)): continue
                            if not mutation.mutationIsPossibleAtVertex(A, v): res['gate refuses'] += 1; continue
                            dj = dims(A, v); bad = [x for x in dj.values() if x[0] >= 3 and x[1] >= 1]
                            if not bad: res['gate ok, no (d>=3,J!=0)'] += 1; continue
                            inl = search._coxeterKeyOrNone(A) in ks
                            res[('d>=3,J!=0', 'LNAkey' if inl else 'noLNAkey')] += 1
                            if inl: ex.setdefault(fingerprint.canonicalKey(A), (arrows, rels, dj, v))
                        except Exception as e: res['err ' + type(e).__name__] += 1
    print('n', n, dict(res), 'distinct LNAkey hits', len(ex)); import json; json.dump(list(ex.values()), open('workshop/rounds/037/theorist_d3_hits_n%d.json' % n, 'w'))
