"""Round 109 (toolsmith), H-027 / T10: audit of the E-287 tilting-exchange paths B_n -> A_n.  For every intermediate T on every SHORTEST path
(all first-depth hits), present End(T) as a quiver with relations (thin, from chain-map composites in homotopy classes), say whether it is an LNA /
gentle LNA, and run the repo's right and left single mutation at the exchanged vertex: E-187 gate (mutationIsPossibleAtVertex), the reduced result,
whether it equals End(T') (same labelling, Cartan), and the Coxeter-key guard (key kept?).  Also the control: all one-step vertex mutations of B_n.
Usage: toolsmith_pathaudit.py n depth   (repository root)"""
import sys, itertools, json, time, io, contextlib
sys.path.insert(0, '.'); sys.path.insert(0, 'tests')
src = open('workshop/rounds/107/toolsmith_tiltsearch.py').read()
exec(src.split("\nif __name__ == '__main__':")[0])
exec('def homotopyImage' + src.split('def homotopyImage', 1)[1].split("def verifyAlgebraA")[0])
import quivermutation as qm
from quivermutation import mutation, pathAlgebra, reduction, invariants, search
from helpers import path_algebra, quiet

n = int(sys.argv[1]); depth = int(sys.argv[2]); W = 7
import os
def build(n):      # as 107's, but objs from 109 when 107 has no file (n = 7)
    alg = Bn(n); A = Alg(alg)
    f = 'workshop/rounds/107/toolsmith_objs_%d.json' % n
    if not os.path.exists(f): f = 'workshop/rounds/109/toolsmith_objs_%d.json' % n
    objs = [deser(S) for S in json.load(open(f))]; N = len(objs); J = 14
    h = [[{j: homdim(A, objs[r], shiftC(objs[s], j), 0) for j in range(-J, J + 1)} for s in range(N)] for r in range(N)]
    return A, objs, h, J
A, objs, h, J = build(n); N = len(objs); print('n', n, 'objects', N, 'Hom tables done'); sys.stdout.flush()
cls = [np.array(klassC(X, n)) for X in objs]
def hom(x, y, m=0):
    (r, a), (s, b) = x, y; j = b - a + m
    return h[r][s][j] if -J <= j <= J else 0
cache = {}
def cp(x, y):
    if (x, y) not in cache: cache[(x, y)] = all(hom(x, y, m) == 0 and hom(y, x, m) == 0 for m in range(-2 * W - 2, 2 * W + 3) if m != 0)
    return cache[(x, y)]
selfok = lambda x: all(hom(x, x, m) == 0 for m in range(-8, 9) if m != 0) and hom(x, x, 0) == 1
klass = lambda x: ((-1) ** x[1]) * cls[x[0]]
P = {}
for r, X in enumerate(objs):
    if len(X[0]) == 1 and len(X[0][0]) == 1: P[X[0][0][0]] = (r, 0)
T0 = tuple(P[i] for i in range(1, n + 1))          # ordered slots 0..n-1 = vertices 1..n
allz = [(r, a) for r in range(N) for a in range(-W, W + 1)]
def isTilt(T): return all(cp(x, y) for x, y in itertools.combinations(T, 2)) and abs(round(np.linalg.det(np.array([klass(x) for x in T])))) == 1
CA = np.array(invariants.cartanMatrix(An(n), exact=True).tolist(), dtype=int)
def Hm(T): return np.array([[hom(x, y) for y in T] for x in T])
def isA(T):
    M = Hm(T)
    return sorted(M.flatten()) == sorted(CA.flatten()) and (perm_canon_eq(M, CA) or perm_canon_eq(M.T, CA))
def nbrs(T):                      # ordered: replace slot i
    out = []
    for i, U in enumerate(T):
        rest = [x for x in T if x != U]
        for y in allz:
            if y == U or y in rest or not selfok(y): continue
            if all(cp(y, x) for x in rest):
                T2 = T[:i] + (y,) + T[i + 1:]
                if abs(round(np.linalg.det(np.array([klass(x) for x in T2])))) == 1: out.append((i, T2))
    return out
# BFS on unordered sets, remembering all shortest-path parents (as ordered tuples)
key_ = lambda T: frozenset(T)
layers = [{key_(T0): [T0]}]; seen = {key_(T0)}; parents = {key_(T0): []}; ordered = {key_(T0): T0}
hitdepth = None
for d in range(1, depth + 1):
    new = {}
    for ks, T in [(k, ordered[k]) for k in layers[-1]]:
        for i, T2 in nbrs(T):
            k2 = key_(T2)
            if k2 in seen and k2 not in new: continue
            if k2 not in new: new[k2] = T2; ordered.setdefault(k2, T2)
            parents.setdefault(k2, []).append((ks, i))
    seen |= set(new); layers.append(new)
    hits = [k for k in new if isA(new[k])]
    print('depth', d, 'new', len(new), 'A_n hits', len(hits)); sys.stdout.flush()
    if hits: hitdepth = d; break
hits = [k for k in layers[-1] if isA(ordered[k])] if hitdepth else []
# enumerate all shortest paths (as lists of (setkey, slot)) to each first-depth hit
def paths(k):
    if not parents[k]: return [[k]]
    return [p + [k] for pk, i in parents[k] for p in paths(pk)]
# ---- presentation of End(T)
hv = lambda T, g=None: None
def presentation(T):
    T = list(T); Z = [shiftC(objs[r], a) for r, a in T]; M = Hm(T)
    if M.max() > 1: return None, 'not thin: Cartan max %d' % M.max()
    gen = {}
    for x in range(n):
        for y in range(n):
            if x != y and M[x][y]:
                if M[y][x]: return None, 'cyclic'
                gen[(x, y)] = generator(A, Z[x], Z[y])
    comp = {}
    def cmp_(x, z, y):
        if (x, z, y) not in comp:
            g = compose(A, Z[x], Z[z], Z[y], gen[(x, z)], gen[(z, y)]); comp[(x, z, y)] = bool(g) and not isNullHomotopic(A, Z[x], Z[y], g)
        return comp[(x, z, y)]
    arrows = [(x, y) for (x, y) in gen if not any((x, z) in gen and (z, y) in gen and cmp_(x, z, y) for z in range(n))]
    out = {}
    for a in arrows: out.setdefault(a[0], []).append(a[1])
    def allpaths(x):
        res = []
        def go(p):
            for y in out.get(p[-1], []): res.append(p + (y,)); go(p + (y,))
        go((x,)); return res
    pm = {}
    def pmap(p):         # chain map of the path composite (or None if zero), memo
        if len(p) == 2: return gen[(p[0], p[1])]
        if p not in pm:
            f = pmap(p[:-1]); 
            if f is None: pm[p] = None
            else:
                g = compose(A, Z[p[0]], Z[p[-2]], Z[p[-1]], f, gen[(p[-2], p[-1])])
                pm[p] = g if (g and not isNullHomotopic(A, Z[p[0]], Z[p[-1]], g)) else None
        return pm[p]
    rels = []; sign_issue = False
    for x in range(n):
        ps = allpaths(x); byend = {}
        for p in ps: byend.setdefault(p[-1], []).append(p)
        for y, pl in byend.items():
            nz = [p for p in pl if pmap(p) is not None]
            for p in pl:
                if pmap(p) is None and (len(p) == 2 or (pmap(p[:-1]) is not None and pmap(p[1:]) is not None)): rels.append([p])
            for p in nz[1:]:
                # ratio of composites to the generator is +-1? compare p with nz[0]
                f, g = pmap(nz[0]), pmap(p)
                if not isNullHomotopic(A, Z[x], Z[y], {k: {ab: {q: v - w for q in set(c) | set(g.get(k, {}).get(ab, {})) for v, w in [(c.get(q, 0), g.get(k, {}).get(ab, {}).get(q, 0))]} for ab, c in m.items()} for k, m in f.items()}):
                    if not isNullHomotopic(A, Z[x], Z[y], {k: {ab: {q: v + w for q in set(c) | set(g.get(k, {}).get(ab, {})) for v, w in [(c.get(q, 0), g.get(k, {}).get(ab, {}).get(q, 0))]} for ab, c in m.items()} for k, m in f.items()}): sign_issue = True
                    else: sign_issue = 'minus'
                rels.append([nz[0], p])
    return (arrows, rels, sign_issue), None
def repoAlg(pres, transp):
    arrows, rels, _ = pres
    f = (lambda p: tuple(v + 1 for v in p[::-1])) if transp else (lambda p: tuple(v + 1 for v in p))
    ar = [(b + 1, a + 1) if transp else (a + 1, b + 1) for a, b in arrows]
    return path_algebra(ar, [[list(f(p)) for p in r] for r in rels])
def describe(pres):
    arrows, rels, si = pres
    ind = [sum(1 for a in arrows if a[1] == v) for v in range(n)]; outd = [sum(1 for a in arrows if a[0] == v) for v in range(n)]
    line = len(arrows) == n - 1 and max(ind) <= 1 and max(outd) <= 1
    mono = all(len(r) == 1 for r in rels)
    lens = sorted(len(r[0]) - 1 for r in rels if len(r) == 1)
    s = 'arrows %s rels %s' % (sorted(arrows), [[ ''.join(map(str, [v + 1 for v in p])) for p in r] for r in rels])
    cls_ = ('LNA' + (' gentle' if line and mono and set(lens) <= {2} else ' non-gentle')) if line and mono else ('line, non-monomial' if line else ('not a line (%d arrows, max out %d in %d)' % (len(arrows), max(outd), max(ind))))
    return cls_, s, ('' if si is False else ' [sign issue %s]' % si)
def key_of(a):
    return search._coxeterKeyOrNone(a)
def cartanOf(a):
    return np.array(invariants.integerCartanMatrix(a))
def stepAudit(T, i, T2, pres, label):
    """repo single mutation at vertex i+1 (right, left) from End(T): gate, matches End(T2)?, key kept?"""
    rows = []
    H2 = Hm(T2)
    for transp in (False,):      # conv0 = End(B_n) reads as the repo's B_n; conv1 is the dual reading
        a = repoAlg(pres, transp)
        if not (pres[2] is False): continue
        k0 = key_of(a)
        for side in ('right', 'left'):
            b = pathAlgebra.dualPathAlgebra(a) if side == 'left' else a
            try: gate = mutation.mutationIsPossibleAtVertex(b, i + 1)
            except Exception as e: gate = 'RAISED ' + type(e).__name__
            res = None; match = None; kept = None
            if gate is True:
                try:
                    f = mutation.leftQuiverMutationAtVertex if side == 'left' else mutation.quiverMutationAtVertex
                    raw = quiet(f, a, i + 1); red = quiet(reduction.reducePathAlgebra, raw)
                    C = cartanOf(red); nodes = sorted(red.quiver.nodes)
                    match = bool(C.shape == H2.shape and ((C == H2).all() or (C == H2.T).all())) if nodes == list(range(1, n + 1)) else False
                    k1 = key_of(red); kept = (k1 == k0) if k1 is not None else None
                except Exception as e: res = 'RAISED ' + type(e).__name__ + ' ' + str(e)[:60]
            rows.append(('conv%d' % transp, side, gate, match, kept, res))
    return rows

if __name__ == '__main__':
    print('CONTROL: one-step vertex mutations of B_%d in the repo (gate / equals a tilting neighbour End / key kept)' % n)
    pres0, why = presentation(T0); print(' End(B_n):', describe(pres0)[:2])
    tiltN = {i: T2 for i, T2 in nbrs(T0)}
    for i in range(n):
        nb = [T2 for ii, T2 in nbrs(T0) if ii == i]
        for side in ('right', 'left'):
            ms = [stepAudit(T0, i, T2, pres0, 'B') for T2 in nb]
            rr = [r for rows in ms for r in rows if r[1] == side]
            gate = rr[0][2] if rr else None
            match = any(r[3] for r in rr); kept = [r[4] for r in rr if r[3]]
            print('  vertex %d %s: gate %s | tilting exchanges at this slot: %d | repo result = End of one of them: %s | key kept (matched) %s' % (i + 1, side, gate, len(nb), match, kept[:1]))
    for k in hits:
        for pth in paths(k):
            print('=== PATH to hit', sorted(ordered[k]))
            cur = T0
            for step in range(len(pth) - 1):
                nxtk = pth[step + 1]; slot = [i for pk, i in parents[nxtk] if pk == pth[step]][0]
                Tn = ordered[pth[step]]; T2 = Tn[:slot] + (None,) + Tn[slot + 1:]
                # find ordered tuple for nxt consistent with slot
                newobj = [x for x in ordered[nxtk] if x not in Tn][0]; T2 = Tn[:slot] + (newobj,) + Tn[slot + 1:]
                pres, why = presentation(Tn)
                if pres is None: print(' step', step, 'End(T): ', why); cur = T2; continue
                c, s, si = describe(pres)
                print(' T_%d End(T): %s | %s%s' % (step, c, s, si))
                print('   exchange slot %d (vertex %d): %s -> %s, tilting %s' % (slot, slot + 1, Tn[slot], newobj, isTilt(T2)))
                for r in stepAudit(Tn, slot, T2, pres, ''): print('    ', r)
            pres, why = presentation(T2 if False else ordered[pth[-1]]); 
            if pres: print(' final End(T):', describe(pres)[0], describe(pres)[1])
