"""Round 129 (maverick), revision of round 127 (T10 / H-032, E-331).  Three instruments, same step semantics as round 127
`maverick_pureF.py` (gate, illegal drop, reducePathAlgebra, Coxeter-key guard, labelled canonicalKey dedup), which this file reuses by exec.

  ctrl7 N CLS D : POWER control for E-331's forward detector.  Multi-seed F-only BFS from all LNAs + duals of key class CLS at length N
                  (classes ordered as rounds/045 toolsmith_guardaudit.py and rounds/070 toolsmith_cyccount.py: by (size, str(key))),
                  depth D, every key-kept F step followed (the library walk).  Prints every key-kept J != 0 step with its seed and path.
                  E-204 (independent code, tiltingPlus) counts 4 key-kept failing steps in this ball at N = 7, CLS = 1, D = 8.
  inv N D ROOT  : exhaustive edge-inversion count.  F-only ball (root + dual) to depth D; for EVERY gate-admitted key-kept J = 0 F edge
                  X -> Y out of a node at level d - 1, look among the gate-admitted key-kept R children of Y for one isomorphic to X
                  (canonicalKey equality, else explicit sign-aware iso).  Tally by d; for a failure, diagnose: R children of Y without the
                  J = 0 filter / key guard, and raw R mutations of Y with the gate and the illegal-relation drop switched off.
  cls N D ROOTS : forward F-only balls (root + dual, every key-kept F step followed) from the given LNAs; every node is iso-tested
                  against the 13 known F-type crossing parents of round 125 at length N; key-kept J != 0 steps are counted.

usage (repo root):
  timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py ctrl7 7 1 8
  timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py inv 11 3 003033030
  timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py cls 11 5 003033030,203033030
  timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py diag 11 003033030 "root F3 F3 + F4"   # one failing edge, all R children of Y
  timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py idealctrl                           # planted controls for sameIdeal
  MK=FR timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py cls 11 3 003033030     # iso-detector control"""
import sys, os, time, collections, copy
sys.path.insert(0, '.'); _MV = sys.argv[:]
sys.argv = ['x']
_p = open('workshop/rounds/127/maverick_pureF.py').read(); exec(compile(_p.split("\nif __name__ == '__main__':")[0], 'pureF', 'exec'))
sys.argv = _MV

def ctrl7(N, cls, D):
    classes = {}
    for lna in nk.LinearNakayamaAlgebra.allOfLength(N):
        for tag, alg in (('L', lna), ('D', pathAlgebra.dualPathAlgebra(lna))):
            classes.setdefault(search._coxeterKeyOrNone(alg), []).append((alg, tag + ':' + str(lna_name(lna, N))))
    order = sorted(classes, key=lambda k: (len(classes[k]), str(k))); base = order[cls]
    starts = classes[base]; print('N=%d class %d key %s: %d seeds, F-only depth %d' % (N, cls, base, len(starts), D), flush=True)
    hits = []
    def on_jx(p, c, path, d): hits.append((d, path)); print('  KEY-KEPT J!=0 at parent level %d: seed %s path %s' % (d - 1, path[0], ' '.join(path[1:])), flush=True)
    t0 = time.time(); sizes, T = bfs(starts, 'F', D, base, on_jx=on_jx)
    print('CTRL7 N=%d cls %d depth %d: nodes %d, steps %s, key-kept J!=0 %d  %.0fs' % (N, cls, D, sum(sizes), dict(T), len(hits), time.time() - t0))

def same(X, Y):
    kx, ky = fingerprint.canonicalKey(X), fingerprint.canonicalKey(Y)
    if kx is not None and kx == ky: return True
    if X.quiver.number_of_edges() != Y.quiver.number_of_edges() or len(X.quiver.nodes) != len(Y.quiver.nodes): return False
    return H(X) == H(Y) and iso(X, Y, limit=100000)

def qgraph(a):
    g = nx.DiGraph()
    for v in a.quiver.nodes(): g.add_node(('v', v), lab='v')
    for (x, y, k) in a.quiver.edges(keys=True):
        g.add_node(('a', x, y, k), lab='a'); g.add_edge(('a', x, y, k), ('v', x), lab='s'); g.add_edge(('a', x, y, k), ('v', y), lab='t')
    return g

def sameIdeal(X, Z, limit=5000):
    """X and Z are the same algebra in a weak, checkable sense: some quiver isomorphism (no arrow rescaling) carries the ideal of Z
    onto the ideal of X (every relation of each lies in the ideal of the other; arrowPaths.isInIdeal).  Generating sets may differ."""
    if X.quiver.number_of_edges() != Z.quiver.number_of_edges() or len(X.quiver.nodes) != len(Z.quiver.nodes): return False
    rX, rZ = procedure.relationsFrom(X), procedure.relationsFrom(Z)
    gm = nx.isomorphism.DiGraphMatcher(qgraph(Z), qgraph(X), node_match=lambda a, b: a['lab'] == b['lab'], edge_match=lambda a, b: a['lab'] == b['lab'])
    for t, m in enumerate(gm.isomorphisms_iter()):
        if t >= limit: return False
        am = {k[1:]: v[1:] for k, v in m.items() if k[0] == 'a'}; vm = {k[1]: v[1] for k, v in m.items() if k[0] == 'v'}
        rZX = [{tuple(am[a] for a in path): c for path, c in rel.items()} for rel in rZ]
        inv_am = {v: k for k, v in am.items()}
        rXZ = [{tuple(inv_am[a] for a in path): c for path, c in rel.items()} for rel in rX]
        if all(arrowPaths.isInIdeal(X.quiver, rX, r) for r in rZX if r) and all(arrowPaths.isInIdeal(Z.quiver, rZ, r) for r in rXZ if r): return True
    return False

def allR(Y, pk0, gate=True):
    """R children of Y: (move, child, J, keykept), no J filter, no key filter; gate=False also skips the gate and the illegal drop."""
    a = pathAlgebra.dualPathAlgebra(Y); out = []
    for v in sorted(a.quiver.nodes):
        try:
            if gate and not mutation.mutationIsPossibleAtVertex(a, v): continue
            raw = mutation.quiverMutationAtVertex(a, v)
            if gate and any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = pathAlgebra.dualPathAlgebra(reduction.reducePathAlgebra(raw))
            ck = search._coxeterKeyOrNone(ch)
            out.append(('R' + str(v), ch, bool(perI(a, v)) if gate else None, ck == pk0))
        except Exception as e: out.append(('R%s:exc %s' % (v, type(e).__name__), None, None, False))
    return out

def invcount(N, D, rt):
    r = nk.LinearNakayamaAlgebra.fromClassName(rt); pk0 = search._coxeterKeyOrNone(r); t0 = time.time()
    seen = set(); fr = []
    for a, tag in ((copy.deepcopy(r), 'root'), (pathAlgebra.dualPathAlgebra(r), 'dual')):
        k = fingerprint.canonicalKey(a)
        if k not in seen: seen.add(k); fr.append((a, (tag,)))
    tot = collections.Counter()
    for d in range(1, D + 1):
        T = collections.Counter(); nf = []
        for X, path in fr:
            for mv, Y, J in children(X, 'F', pk0):
                if J: T['F J!=0 edges (not tested)'] += 1; continue
                T['F J=0 edges'] += 1
                Rs = [(m, c) for m, c, Jr in children(Y, 'R', pk0) if not Jr]
                hit = [m for m, c in Rs if same(X, c)]
                if hit: T['inverse J=0 R found'] += 1
                elif any(sameIdeal(X, c) for m, c in Rs):
                    T['inverse only up to ideal (other generators)'] += 1
                    print('  INVERSE ONLY AS IDEAL d=%d: %s %s + %s at %s' % (d, path[0], ' '.join(path[1:]), mv, [m for m, c in Rs if sameIdeal(X, c)]), flush=True)
                else:
                    T['NO inverse'] += 1
                    full = allR(Y, pk0); raw = allR(Y, pk0, gate=False)
                    diag = dict(arrowsX=X.quiver.number_of_edges(), arrowsY=Y.quiver.number_of_edges(),
                                Rgated=['%s%s%s' % (m, '' if kk else ' keymoved', ' J!=0' if Jr else '') for m, c, Jr, kk in full],
                                isoAmongGatedAnyJ=[m for m, c, Jr, kk in full if c is not None and same(X, c)],
                                isoAmongUngated=[m for m, c, Jr, kk in raw if c is not None and same(X, c)])
                    print('  NO INVERSE d=%d: %s %s + %s  %s' % (d, path[0], ' '.join(path[1:]), mv, diag), flush=True)
                k = fingerprint.canonicalKey(Y)
                if k is None or k in seen: continue
                seen.add(k); nf.append((Y, path + (mv,)))
        fr = nf; tot.update(T)
        print('ROOT %s depth %d %s  new nodes %d  %.0fs' % (rt, d, dict(T), len(nf), time.time() - t0), flush=True)
    print('INV TOTAL %s N=%d depth %d: %s' % (rt, N, D, dict(tot)))

def sig(a):
    q = a.quiver
    return (q.number_of_edges(), tuple(sorted((q.in_degree(v), q.out_degree(v)) for v in q.nodes)))

MK = os.environ.get('MK', 'F')   # cls: moves followed; MK=FR is the iso-detector control (F/R balls reach the parents at depth 3)
def cls(N, D, roots):
    P = targets(N); print('%d crossing parents (n = %d) up to iso' % (len(P), N), flush=True)
    Psig = collections.defaultdict(list)
    for i, (p, tag) in enumerate(P): Psig[sig(p)].append((i, p, H(p)))
    tot = collections.Counter()
    for rt in roots:
        r = nk.LinearNakayamaAlgebra.fromClassName(rt); pk0 = search._coxeterKeyOrNone(r); t0 = time.time(); hits = []; jx = []; ntest = [0]
        def on_node(a, path, d):
            c = Psig.get(sig(a))
            if not c: return
            h = H(a); ntest[0] += 1   # struct iso, else the same ideal under a quiver iso (generating sets may differ)
            for i, p, hp in c:
                if (hp == h and iso(a, p, limit=100000)) or sameIdeal(a, p): hits.append((i, d, ' '.join(path))); print('   PARENT #%d reached at depth %d path %s' % (i, d, ' '.join(path)), flush=True)
        def on_jx(p, c, path, d): jx.append((d, ' '.join(path))); print('   key-kept J!=0 at depth %d path %s' % (d, ' '.join(path)), flush=True)
        sizes, T = bfs([(copy.deepcopy(r), 'root'), (pathAlgebra.dualPathAlgebra(r), 'dual')], MK, D, pk0, on_jx=on_jx, on_node=on_node, tlimit=560)
        cap = [k for k in T if k.startswith('TIME')]
        print(('  ROOT %s (n=%d) moves '+MK+' depth %d%s: nodes %d, steps %s, signature-matched nodes iso-tested %d, parents reached %s, key-kept J!=0 %d  %.0fs')
              % (rt, N, D, ' ' + cap[0] if cap else '', sum(sizes), {k: v for k, v in T.items() if not k.startswith('TIME')}, ntest[0], hits or 'NONE', len(jx), time.time() - t0), flush=True)
        tot['roots'] += 1; tot['nodes'] += sum(sizes); tot['hits'] += len(hits); tot['J!=0'] += len(jx)
    print('CLS TOTAL', dict(tot))

def replay(r, path):
    A = copy.deepcopy(r) if path[0] == 'root' else pathAlgebra.dualPathAlgebra(r)
    for mv in path[1:]:
        A = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(A, int(mv[1:])))
    return A

def diag(N, rt, spec):
    """spec 'root F3 F3 + F4': replay X = root.F3.F3, Y = X.F4; describe Y's R children against X."""
    r = nk.LinearNakayamaAlgebra.fromClassName(rt); pk0 = search._coxeterKeyOrNone(r)
    a, b = spec.split('+'); X = replay(r, a.split()); Y = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(X, int(b.strip()[1:])))
    CX = invariants.integerCartanMatrix(X); dimX = sum(int(x) for row in CX for x in row)
    print('EDGE %s %s: X arrows %d rels %d dim %d | Y arrows %d rels %d; J(F) = %s' % (rt, spec, X.quiver.number_of_edges(), len(procedure.relationsFrom(X)), dimX,
          Y.quiver.number_of_edges(), len(procedure.relationsFrom(Y)), bool(perI(X, int(b.strip()[1:])))))
    for gate in (True, False):
        for m, Z, Jr, kk in allR(Y, pk0, gate=gate):
            if Z is None: print('   gate=%s %s' % (gate, m)); continue
            CZ = invariants.integerCartanMatrix(Z); dimZ = sum(int(x) for row in CZ for x in row)
            qi = nx.is_isomorphic(nx.MultiDiGraph(X.quiver), nx.MultiDiGraph(Z.quiver))
            back = [w for w in sorted(Z.quiver.nodes) if mutation.mutationIsPossibleAtVertex(Z, w) and same(Y, reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(Z, w)))]
            print('   gate=%s %s J=%s keykept=%s: arrows %d rels %d dim %d quiver-iso-to-X %s SNF,sig=X %s iso-to-X %s same-ideal %s; F back to Y at %s' % (gate, m, Jr, kk,
                  Z.quiver.number_of_edges(), len(procedure.relationsFrom(Z)), dimZ, qi, inv(Z) == inv(X), same(X, Z), sameIdeal(X, Z), back))

def idealctrl():
    """Planted controls for sameIdeal on X = 003033030.root.F3.F3 (the first non-invertible edge's parent)."""
    r = nk.LinearNakayamaAlgebra.fromClassName('003033030'); X = replay(r, ['root', 'F3', 'F3']); R = procedure.relationsFrom(X)
    def withrels(A, rels): return procedure.toPathAlgebra(A.quiver.copy(), [dict(q) for q in rels])
    a = lambda s, t: (s, t, 0)
    old3 = {(a(5, 6), a(6, 7), a(7, 8)): 1}; alt3 = {(a(5, 6), a(6, 3), a(3, 8)): 1}
    m2 = {(a(3, 8), a(8, 9)): 1}; m2b = {(a(8, 9), a(9, 10)): 1}
    assert old3 in R and m2 in R
    perm = {v: (v % 11) + 1 for v in X.quiver.nodes}   # cyclic relabelling of the vertices
    def relabel(A, rels):
        return procedure.toPathAlgebra(nx.relabel_nodes(A.quiver, perm, copy=True),
                                       [{tuple((perm[x], perm[y], k) for (x, y, k) in path): c for path, c in rel.items()} for rel in rels])
    cases = [('same generators (truth YES)', withrels(X, list(R)), True),
             ('length-3 generator replaced by the one equal mod the commutativity (truth YES)', withrels(X, [alt3 if q == old3 else q for q in R]), True),
             ('same, vertices cyclically relabelled (truth YES)', relabel(X, [alt3 if q == old3 else q for q in R]), True),
             ('one relation dropped (truth NO)', withrels(X, [q for q in R if q != old3]), False),
             ('monomial (3,8)(8,9) moved to (8,9)(9,10) (truth NO)', withrels(X, [m2b if q == m2 else q for q in R]), False),
             ('commutativity sign flipped (NO for this strict test; isomorphic by rescaling)', withrels(X, [{k: 1 for k in q} if len(q) == 2 else q for q in R]), False)]
    ok = 0
    for nm, B, want in cases:
        assert procedure.relationsFrom(B) == [dict(q) for q in B.arrowRels]   # the planted generators are the ones read back
        got = sameIdeal(X, B); ok += got == want; print('  %-85s sameIdeal %s  struct-iso %s' % (nm, got, same(X, B)))
    print('IDEALCTRL %d of %d as expected' % (ok, len(cases)))

def jx(N, rt, spec):
    """spec 'root F1 F1 F4 F3 F1 + F8': replay the parent, re-check the step (gate, J, tiltingPlus, key) and compare SNF(S), signature."""
    r = nk.LinearNakayamaAlgebra.fromClassName(rt); pk0 = search._coxeterKeyOrNone(r)
    a, b = spec.split('+'); P = replay(r, a.split()); v = int(b.strip()[1:])
    ok = [mv for mv, c, J in children(P, 'F', pk0) if mv == 'F%d' % v]
    X = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(P, v))
    tp = None   # tiltingPlus (scholar_h015) = mutation.isTiltingAtVertex (E-327); the library one is printed
    print('JX %s %s: gate+legal+key-kept %s, J!=0 %s, tiltingPlus %s, isTiltingAtVertex %s | parent arrows %d SNF tail %s sig %s | child arrows %d SNF tail %s sig %s | CROSSING %s | parent iso to a known crossing parent %s'
          % (rt, spec, bool(ok), bool(perI(P, v)), tp, mutation.isTiltingAtVertex(P, v), P.quiver.number_of_edges(), inv(P)[0][-3:], inv(P)[1],
             X.quiver.number_of_edges(), inv(X)[0][-3:], inv(X)[1], inv(P) != inv(X), [i for i, (Q, t) in enumerate(targets(N)) if same(P, Q) or sameIdeal(P, Q)]))

def after(N, rt, spec, D, kinds):
    """From the child X of a crossing step (spec as for jx), a J = 0 ball (moves `kinds`) of depth D: which LNAs does it reach?"""
    r = nk.LinearNakayamaAlgebra.fromClassName(rt)
    a, b = spec.split('+'); P = replay(r, a.split()); X = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(P, int(b.strip()[1:])))
    pk0 = search._coxeterKeyOrNone(X); found = {}
    def on_node(A, path, d):
        nm = lna_name(A, N)
        if nm and nm not in found: found[nm] = (d, ' '.join(path[1:]), inv(A)[0][-3:]); print('   LNA %s at depth %d via %s SNF tail %s' % ((nm,) + found[nm]), flush=True)
    sizes, T = bfs([(X, 'X')], kinds, D, pk0, on_node=on_node, j0only=True)
    print('AFTER %s %s: child SNF tail %s; %s J=0 ball depth %d: nodes %d, LNAs %s' % (rt, spec, inv(X)[0][-3:], kinds, D, sum(sizes), found or 'NONE'))

def walk(N, rt, spec):
    """spec 'dual F6 F11 ...': replay step by step with the walk's own step semantics (children(): gate, illegal drop, reduce, key guard);
    print J per step, SNF tail per node, and the LNA name of the end; compare with the SNF of the named LNA built from scratch."""
    r = nk.LinearNakayamaAlgebra.fromClassName(rt); pk0 = search._coxeterKeyOrNone(r); p = spec.split()
    A = copy.deepcopy(r) if p[0] == 'root' else pathAlgebra.dualPathAlgebra(r); out = []
    for mv in p[1:]:
        ch = [(c, J) for m, c, J in children(A, 'F', pk0) if m == mv]
        assert ch, 'step %s not admitted (gate / illegal / key)' % mv
        A, J = ch[0]; out.append('%s%s%s' % (mv, '*' if J else '', inv(A)[0][-3:]))
    nm = lna_name(A, N); L = nk.LinearNakayamaAlgebra.fromClassName(nm) if nm else None
    print('WALK %s %s: %s | end LNA %s, SNF tail (built from name) %s, key = root key %s; start SNF tail %s (* = J != 0)' % (rt, spec, ' '.join(out), nm,
          inv(L)[0][-3:] if L else None, search._coxeterKeyOrNone(L) == pk0 if L else None, inv(r)[0][-3:]))

if __name__ == '__main__':
    m = _MV[1]
    if m == 'ctrl7': ctrl7(int(_MV[2]), int(_MV[3]), int(_MV[4]))
    elif m == 'inv': invcount(int(_MV[2]), int(_MV[3]), _MV[4])
    elif m == 'after': after(int(_MV[2]), _MV[3], _MV[4], int(_MV[5]), _MV[6])
    elif m == 'walk': walk(int(_MV[2]), _MV[3], _MV[4])
    elif m == 'jx': jx(int(_MV[2]), _MV[3], _MV[4])
    elif m == 'idealctrl': idealctrl()
    elif m == 'diag': diag(int(_MV[2]), _MV[3], _MV[4])
    elif m == 'cls': cls(int(_MV[2]), int(_MV[3]), _MV[4].split(','))
