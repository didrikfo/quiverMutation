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
    return (q.number_of_edges(), tuple(sorted((q.in_degree(v), q.out_degree(v)) for v in q.nodes)), len(procedure.relationsFrom(a)))

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
            h = H(a); ntest[0] += 1
            for i, p, hp in c:
                if hp == h and iso(a, p, limit=100000): hits.append((i, d, ' '.join(path))); print('   PARENT #%d reached at depth %d path %s' % (i, d, ' '.join(path)), flush=True)
        def on_jx(p, c, path, d): jx.append((d, ' '.join(path))); print('   key-kept J!=0 at depth %d path %s' % (d, ' '.join(path)), flush=True)
        sizes, T = bfs([(copy.deepcopy(r), 'root'), (pathAlgebra.dualPathAlgebra(r), 'dual')], MK, D, pk0, on_jx=on_jx, on_node=on_node, tlimit=560)
        cap = [k for k in T if k.startswith('TIME')]
        print(('  ROOT %s (n=%d) moves '+MK+' depth %d%s: nodes %d, steps %s, signature-matched nodes iso-tested %d, parents reached %s, key-kept J!=0 %d  %.0fs')
              % (rt, N, D, ' ' + cap[0] if cap else '', sum(sizes), {k: v for k, v in T.items() if not k.startswith('TIME')}, ntest[0], hits or 'NONE', len(jx), time.time() - t0), flush=True)
        tot['roots'] += 1; tot['nodes'] += sum(sizes); tot['hits'] += len(hits); tot['J!=0'] += len(jx)
    print('CLS TOTAL', dict(tot))

if __name__ == '__main__':
    m = _MV[1]
    if m == 'ctrl7': ctrl7(int(_MV[2]), int(_MV[3]), int(_MV[4]))
    elif m == 'inv': invcount(int(_MV[2]), int(_MV[3]), _MV[4])
    elif m == 'cls': cls(int(_MV[2]), int(_MV[3]), _MV[4].split(','))
