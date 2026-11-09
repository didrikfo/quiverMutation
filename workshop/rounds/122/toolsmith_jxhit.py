"""Round 122 (toolsmith): check the key-kept J != 0 steps that `toolsmith_jxsearch.py` found at n = 11 (group 5: root 003033030, 13-LNA class)
with INV(parent) != INV(child).  Rebuilds each parent from its printed path, re-takes the step, and checks independently of the library's Cartan code:
  (1) the Cartan matrix by brute-force linear algebra (paths i -> j modulo the span of u * rel * w), compared with `integerCartanMatrix`;
  (2) det C, acyclicity (finite global dimension), Coxeter polynomial of both;
  (3) SNF(S), S = C + C^T, by sympy and by determinantal divisors (gcd of k x k minors, k = n-2, n-1), and the LNA class each INV belongs to;
  (4) gate / perI / key at the step;
  (5) a J = 0 key-kept BFS from the child to radius R (default 4) listing the LNAs it reaches (and their INV class): a reached LNA of the OTHER
      class would mean the guarded walk joins two derived-inequivalent LNAs.
usage (repo root): R=4 timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxhit.py"""
import sys, os, itertools, math, time, collections
sys.path.insert(0, '.'); sys.argv = ['x', 'none']
_s = open('workshop/rounds/122/toolsmith_jxsearch.py').read(); exec(compile(_s.split("if __name__ == '__main__':")[0], 'jx', 'exec'))
N = 11; ROOT = '003033030'; R = int(os.environ.get('R', '4'))
STEPS = [('root', 'F3 R8 F3', 'F3'), ('root', 'F3 R8 F3', 'R5'), ('dual', 'F8 R3 R3', 'F5'), ('dual', 'F8 R3 R3', 'R3')]

def step(alg, mv):
    kind, v = mv[0], int(mv[1:]); a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
    assert mutation.mutationIsPossibleAtVertex(a, v), 'gate refuses'
    raw = mutation.quiverMutationAtVertex(a, v)
    assert not any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw))
    ch = reduction.reducePathAlgebra(raw)
    return (pathAlgebra.dualPathAlgebra(ch) if kind == 'R' else ch), bool(perI(a, v))

def paths(Q, i, j):
    if i == j: return [()]
    out = []
    for (x, y, k) in Q.out_edges(i, keys=True):
        out += [((x, y, k),) + p for p in paths(Q, y, j)]
    return out

def brute_cartan(alg):
    """C[j][i] = dim e_j A e_i as integerCartanMatrix (row j, column i): paths i -> j modulo span{u rel w}."""
    Q = alg.quiver; rels = procedure.relationsFrom(alg); V = sorted(Q.nodes); C = [[0] * len(V) for _ in V]
    for ii, i in enumerate(V):
        for jj, j in enumerate(V):
            P = paths(Q, i, j); idx = {p: t for t, p in enumerate(P)}
            if not P: continue
            rows = []
            for rel in rels:
                (p0, _), = list(rel.items())[:1]; s, t = p0[0][0], p0[-1][1]
                for u in paths(Q, i, s):
                    for w in paths(Q, t, j):
                        vec = [0] * len(P)
                        for p, c in rel.items(): vec[idx[u + tuple(p) + w]] += c
                        if any(vec): rows.append(vec)
            rk = sp.Matrix(rows).rank() if rows else 0
            C[jj][ii] = len(P) - rk
    return C

def detdiv(S, k):
    n = len(S); M = sp.Matrix(S); g = 0
    for r in itertools.combinations(range(n), k):
        for c in itertools.combinations(range(n), k):
            g = math.gcd(g, int(M.extract(list(r), list(c)).det()))
            if g == 1: return 1
    return g

lab = labels(N); cls = groups(N)[lab[ROOT]['group']][1]
clsinv = {ci: INV(ct.lnaCartanMatrix(N, tuple(int(c) for c in m[0]))) for ci, m in cls.items()}
which = lambda A: [ci for ci, v in clsinv.items() if v == INV(invariants.integerCartanMatrix(A))]
r = nk.LinearNakayamaAlgebra.fromClassName(ROOT); pk0 = search._coxeterKeyOrNone(r)
print('root', ROOT, 'class sizes', {ci: len(m) for ci, m in cls.items()}, 'class INV SNF tails', {ci: v[0][-3:] for ci, v in clsinv.items()})
t0 = time.time(); kids = []
for side, path, st in STEPS:
    x = copy.deepcopy(r) if side == 'root' else pathAlgebra.dualPathAlgebra(r); js = []
    for mv in path.split():
        x, J = step(x, mv); js.append(J); assert search._coxeterKeyOrNone(x) == pk0
    c, J = step(x, st); kids.append(c)
    print('\n%s %s then %s: path steps J != 0? %s; step J != 0 %s; key kept %s' % (side, path, st, js, J, search._coxeterKeyOrNone(c) == pk0))
    for nm, A in (('parent', x), ('child', c)):
        C = invariants.integerCartanMatrix(A); B = brute_cartan(A); S = [[C[i][j] + C[j][i] for j in range(N)] for i in range(N)]
        Cm = sp.Matrix(C); x_ = sp.symbols('x'); phi = (-Cm.T * Cm.inv()).charpoly(x_).as_expr() if Cm.det() != 0 else None
        print('  %s: acyclic %s, brute-force Cartan == library %s, det C %s, Coxeter poly %s' % (nm, nx.is_directed_acyclic_graph(A.quiver), B == C, Cm.det(), sp.factor(phi)))
        print('     SNF(S) %s; det divisors d_%d = %s, d_%d = %s; INV class %s' % (INV(C)[0], N - 1, detdiv(S, N - 1), N - 2, detdiv(S, N - 2), which(A)))
print('\nchildren pairwise isomorphic:', [(a, b, iso(kids[a], kids[b])) for a in range(4) for b in range(a + 1, 4)])
# (5) J = 0 key-kept BFS from each child: LNAs reached and their class
for t, c in enumerate(kids[:2]):
    seen = {fingerprint.canonicalKey(c)}; fr = [c]; found = collections.Counter()
    for d in range(1, R + 1):
        nf = []
        for x in fr:
            if list(nx.simple_cycles(x.quiver)): continue
            for kind, v in [(kd, v) for kd in 'FR' for v in sorted(x.quiver.nodes)]:
                a = x if kind == 'F' else pathAlgebra.dualPathAlgebra(x)
                if not mutation.mutationIsPossibleAtVertex(a, v) or perI(a, v): continue
                raw = mutation.quiverMutationAtVertex(a, v)
                if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
                ch = reduction.reducePathAlgebra(raw)
                if kind == 'R': ch = pathAlgebra.dualPathAlgebra(ch)
                if search._coxeterKeyOrNone(ch) != pk0: continue
                k = fingerprint.canonicalKey(ch)
                if k is None or k in seen: continue
                seen.add(k); nf.append(ch); nm = lna_name(ch, N)
                if nm: found[(nm, tuple(which(ch)), d)] += 1
        fr = nf
        print('child %d J = 0 ball radius %d: %d nodes, LNAs (name, class, radius) %s  %.0fs' % (t, d, len(seen), sorted(found), time.time() - t0), flush=True)

if os.environ.get('WALK', '0') == '1':
    # the explicit 5-step walk: root -(J=0)^3-> parent -(J != 0, key kept)-> child 0 -(J=0)-> an LNA; every step gate-admitted, key-kept
    c = kids[0]
    for kind, v in [(kd, v) for kd in 'FR' for v in sorted(c.quiver.nodes)]:
        a = c if kind == 'F' else pathAlgebra.dualPathAlgebra(c)
        if not mutation.mutationIsPossibleAtVertex(a, v): continue
        try: ch, J = step(c, kind + str(v))
        except AssertionError: continue
        nm = lna_name(ch, N)
        if nm and search._coxeterKeyOrNone(ch) == pk0:
            print('WALK root %s: F3 R8 F3 (J=0) F3 (J != 0) %s%d (J != 0? %s) -> LNA %s, INV class %s, SNF tail %s' % (ROOT, kind, v, J, nm, which(ch), INV(invariants.integerCartanMatrix(ch))[0][-3:]))
