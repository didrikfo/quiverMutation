"""Round 127 (maverick): T10 / H-032.  Is there a PURE-F route (the library walk: right mutations only, from an LNA or from its dual,
gate + illegal-drop + reducePathAlgebra + Coxeter-key guard, every key-kept child expanded) to a key-kept J != 0 crossing step?

Two instruments (same step semantics as round-122 `toolsmith_jxsearch.ball`, KINDS restricted; labelled canonicalKey dedup = ISO=0):
  fwd  : F-only balls from LNA roots (and their duals); tally gate-admitted key-kept steps, key-kept J != 0 steps, and crossings
         (SNF(S) tail or signature of S = C + C^T differ between parent and child; flint SNF).
  back : R-only balls from the KNOWN F-type crossing parents of E-318 / E-323 (round-125 `toolsmith_jxclass_out.jsonl`; an R-type
         crossing parent is replaced by its dual, which has the mirrored F step), deduplicated up to isomorphism.  If an LNA L reaches the
         parent P by k tilting F steps, then P reaches L (or dual L) by k R steps -- PROVIDED the inverse of each gate-admitted F step is a
         gate-admitted R step in the repo's implementation.  That proviso is what `ctrl` tests.  J0=1 (default) follows only J = 0 R steps:
         a J != 0 R step need not invert an F step (J0=0 finds LNAs of the OTHER SNF class through the known mirror crossing step).
  ctrl : power control for `back`: take nodes at F-depth exactly k from an LNA root, run the R-only ball of radius k from each, and
         check that the root (or its dual) is found.  Must say YES on every node (truth YES).

usage (repo root):
  timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py fwd N DEPTH root1,root2,...
  [J0=0|1] [TGT=i,j] timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py back DEPTH [N]
  timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py ctrl ROOT K [SAMPLE]"""
import sys, os, time, json, collections, random, copy
sys.path.insert(0, '.'); _MA = sys.argv[:]
os.environ.setdefault('ISO', '0'); os.environ['EXT'] = '0'
sys.argv = ['x', 'none']
_s = open('workshop/rounds/125/toolsmith_jxclass.py').read(); exec(compile(_s.split("\nif __name__ == '__main__':")[0], 'jxc', 'exec'))
sys.argv = _MA
import flint

def inv(A):
    C = invariants.integerCartanMatrix(A); n = len(C)
    S = [[int(C[i][j]) + int(C[j][i]) for j in range(n)] for i in range(n)]
    M = flint.fmpz_mat(S).snf(); d = tuple(sorted(abs(int(M[i, i])) for i in range(n)))
    ev = np.linalg.eigvalsh(np.array(S, dtype=float))
    return d, (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))

def children(alg, kinds, pk0):
    """yield (kind+v, child, J) for gate-admitted, legal, acyclic-key, key-kept steps."""
    if list(nx.simple_cycles(alg.quiver)): return
    for kind in kinds:
        a = alg if kind == 'F' else pathAlgebra.dualPathAlgebra(alg)
        for v in sorted(a.quiver.nodes):
            if not mutation.mutationIsPossibleAtVertex(a, v): continue
            raw = mutation.quiverMutationAtVertex(a, v)
            if any(ap.isIllegalRelation(raw.quiver, rr) for rr in procedure.relationsFrom(raw)): continue
            ch = reduction.reducePathAlgebra(raw)
            if kind == 'R': ch = pathAlgebra.dualPathAlgebra(ch)
            ck = search._coxeterKeyOrNone(ch)
            if ck is None or ck != pk0: continue
            yield kind + str(v), ch, bool(perI(a, v))

def bfs(starts, kinds, D, pk0, on_jx=None, on_node=None, tlimit=None, j0only=False):
    """starts: list of (alg, tag).  Returns (layers sizes, steps counter)."""
    seen = set(); fr = []
    for a, tag in starts:
        k = fingerprint.canonicalKey(a)
        if k in seen: continue
        seen.add(k); fr.append((a, (tag,))); on_node and on_node(a, (tag,), 0)
    T = collections.Counter(); sizes = [len(fr)]; t0 = time.time()
    for d in range(1, D + 1):
        nf = []
        for alg, path in fr:
            for mv, ch, J in children(alg, kinds, pk0):
                T['J!=0 kept' if J else 'J=0 kept'] += 1
                if J and on_jx: on_jx(alg, ch, path + (mv,), d)
                if J and j0only: continue
                k = fingerprint.canonicalKey(ch)
                if k is None or k in seen: continue
                seen.add(k); nf.append((ch, path + (mv,))); on_node and on_node(ch, path + (mv,), d)
            if tlimit and time.time() - t0 > tlimit: T['TIME CAP at depth %d' % d] = 1; return sizes + [len(nf)], T
        fr = nf; sizes.append(len(fr))
        print('    depth %d: new nodes %d, cumulative %d, steps %s  %.0fs' % (d, len(fr), len(seen), dict(T), time.time() - t0), flush=True)
    return sizes, T

def fwd(N, D, roots):
    tot = collections.Counter()
    for rt in roots:
        r = nk.LinearNakayamaAlgebra.fromClassName(rt); pk0 = search._coxeterKeyOrNone(r); i0 = inv(r); t0 = time.time()
        rows = []
        def on_jx(p, c, path, d):
            ip, ic = inv(p), inv(c); rows.append((d, ' '.join(path), ip[0][-3:], ic[0][-3:], ip != ic, ic != i0))
        print('ROOT %s (n=%d) SNF tail %s sig %s, F-only depth %d' % (rt, N, i0[0][-3:], i0[1], D), flush=True)
        sizes, T = bfs([(copy.deepcopy(r), 'root'), (pathAlgebra.dualPathAlgebra(r), 'dual')], 'F', D, pk0, on_jx=on_jx)
        for row in rows: print('   J!=0 key-kept: depth %d path %s SNF %s -> %s cross %s child-vs-root differs %s' % row, flush=True)
        print('  ROOT %s: nodes %d, steps %s, J!=0 kept %d, crossings %d  %.0fs' % (rt, sum(sizes), dict(T), len(rows), sum(1 for x in rows if x[4]), time.time() - t0), flush=True)
        tot['roots'] += 1; tot['nodes'] += sum(sizes); tot['steps'] += T['J=0 kept'] + T['J!=0 kept']; tot['J!=0'] += len(rows); tot['cross'] += sum(1 for x in rows if x[4])
    print('TOTAL', dict(tot))

def targets(Nsel=None):
    out = []
    for line in open('workshop/rounds/125/toolsmith_jxclass_out.jsonl'):
        o = json.loads(line)
        if 'distinct' in o or o.get('planted') or not o.get('cross'): continue
        if Nsel and len(o['parent_struct']['arrows']) + 1 != Nsel and Nsel != o.get('n', Nsel): pass
        P = _deser(o['parent_struct'])
        if o['v'][0] == 'R': P = pathAlgebra.dualPathAlgebra(P)
        if Nsel and len(P.quiver.nodes) != Nsel: continue
        if not any(iso(P, Q) for Q, _ in out): out.append((P, o['root'] + ' ' + o['side'] + '|' + o['path'] + ' ' + o['v']))
    return out

J0 = os.environ.get('J0', '1') == '1'   # back: follow only J = 0 R steps (inverses of tilting F steps); J0=0 follows all key-kept R steps

def back(D, Nsel=None):
    T = targets(Nsel); print('%d F-type crossing parents up to iso (n = %s)' % (len(T), Nsel or 'all'), flush=True)
    sel = os.environ.get('TGT')   # comma list of target indices (default all)
    for ti, (P, tag) in enumerate(T):
        if sel and str(ti) not in sel.split(','): continue
        N = len(P.quiver.nodes); pk0 = search._coxeterKeyOrNone(P)
        # sanity: P has a key-kept J != 0 F step whose child changes SNF
        cr = [(mv, inv(P)[0][-3:], inv(c)[0][-3:]) for mv, c, J in children(P, 'F', pk0) if J and inv(P) != inv(c)]
        found = collections.defaultdict(list)
        def on_node(a, path, d):
            nm = lna_name(a, N)
            if nm:
                found[nm].append((d, ' '.join(path[1:]), inv(a)[0][-3:])); print('    LNA hit', nm, found[nm][-1], flush=True)
        print('TARGET #%d n=%d %s: F crossing steps at P %s; R-only ball depth %d, J=0 steps only %s' % (ti, N, tag, cr, D, J0), flush=True)
        sizes, St = bfs([(P, 'P')], 'R', D, pk0, on_node=on_node, j0only=J0)
        print('  LNAs reached by R-only from P (= pure-F routes L -> P): %s' % ({k: v[:2] for k, v in found.items()} or 'NONE'), flush=True)

def ctrl(rt, K, sample):
    r = nk.LinearNakayamaAlgebra.fromClassName(rt); N = len(r.quiver.nodes); pk0 = search._coxeterKeyOrNone(r)
    names = {lna_name(r, N), lna_name(pathAlgebra.dualPathAlgebra(r), N)}
    layer = []
    def on_node(a, path, d):
        if d == K: layer.append((a, path))
    bfs([(copy.deepcopy(r), 'root'), (pathAlgebra.dualPathAlgebra(r), 'dual')], 'F', K, pk0, on_node=on_node)
    random.seed(127); pick = random.sample(layer, min(sample, len(layer)))
    ok = 0; nonlna = 0
    for a, path in pick:
        hit = []
        def on_node2(b, p2, d):
            nm = lna_name(b, N)
            if nm in names: hit.append(d)
        if lna_name(a, N) is None: nonlna += 1
        bfs([(a, 'X')], 'R', K, pk0, on_node=on_node2, j0only=J0)
        ok += bool(hit); print('  node %s: root/dual found by R-only at depths %s' % (' '.join(path), sorted(set(hit))[:3]), flush=True)
    print('CTRL root %s K=%d: %d of %d depth-%d F-nodes (%d non-LNA) lead back to the root by R-only within %d' % (rt, K, ok, len(pick), K, nonlna, K))

if __name__ == '__main__':
    m = _MA[1]
    if m == 'fwd': fwd(int(_MA[2]), int(_MA[3]), _MA[4].split(','))
    elif m == 'back': back(int(_MA[2]), int(_MA[3]) if len(_MA) > 3 else None)
    elif m == 'ctrl': ctrl(_MA[2], int(_MA[3]), int(_MA[4]) if len(_MA) > 4 else 10)
