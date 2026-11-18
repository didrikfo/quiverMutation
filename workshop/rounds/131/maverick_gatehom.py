"""Round 131 (maverick), T5 / R-candidate "the gate decides tilting" (round-128 blocker): library gate
(`mutation.mutationIsPossibleAtVertex`) and `mutation.isTiltingAtVertex` against an independent Hom-level tilting test on
CYCLIC MONOMIAL parents reached by walks from LNAs, split by where the cycle sits relative to the mutation vertex i.

Hom engine = rounds/091/skeptic_convention.py + rounds/118/skeptic_conv3.py (the engine of E-327 / E-334), monomial mode
(paths = nonzero arrow-index tuples, so cycles are fine while the algebra is finite-dimensional).  For an F step at i
on algebra a: T = (+)_{u != i} P_u (+) [P_i -> (+)_{i->j} P_j] (sc.mutate_cx 'R', the E-327 / E-334 convention);
Hom-tilting <=> Hom(T, T[m]) = 0 for m = -2, -1, 1, 2.  An R step is the F step on the dual (as lnaMoves and the library do).
Third reference: the EXACT gate (every nonzero path into i from a vertex != i, not only vertex-simple ones).

The walk: from every LNA on n vertices, F and R steps (R = F on the dual), library gate, `quiverMutationAtVertices` (with its
reduction), key guard on acyclic children (Coxeter key of the root; the guard is inert on a cyclic child, as in the
library); states deduplicated by lnaMoves._stateKey.  Two walk semantics are recorded:
  LIB  = `search.mutationSearchDepthFirst`: never expands a cyclic node (`if depth > 0 and noCycles`), so the gate is
         never evaluated on a cyclic parent; its cyclic nodes are leaves (first cyclic node on a path);
  LNAM = `lnaMoves.movesFrom` / `localMutationSequences` (no key guard there; here with it): expand cyclic nodes too.
Every cyclic monomial node reached (either semantics) is a "walk parent"; at each vertex i with an out-arrow and on each
side (A, A^op) the three tests are compared.

usage (repo root):
  timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py controls
  timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py walk N DEPTH [MAXNODES]
  timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py libcount N DEPTH   (instrumented library walk)
"""
import sys, time, itertools, warnings, collections
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/118'); sys.path.insert(0, 'workshop/rounds/091')
warnings.simplefilter('ignore')
import networkx as nx
import skeptic_conv3 as c3          # patches skeptic_convention's Dmat / comp_classes
import skeptic_convention as sc
from quivermutation import mutation, procedure, pathAlgebra, search, lnaMoves, nakayama as nk

# ---------------------------------------------------------------- algebra <-> arrow model
def model(alg):
    """(V, arrows [(s,t,edge)], rels [tuple of arrow ids]) or None if some relation is not monomial"""
    V = sorted(alg.quiver.nodes); E = list(alg.quiver.edges(keys=True)); ex = {e: k for k, e in enumerate(E)}
    rels = []
    for r in procedure.relationsFrom(alg):
        if len(r) != 1: return None
        p = next(iter(r)); rels.append(tuple(ex[e] for e in p))
    return V, E, rels

def nz(p, rels):
    return not any(p[o:o + len(z)] == z for z in rels for o in range(len(p) - len(z) + 1))

def nonzero_paths(E, rels, maxlen=60):
    out = []; fr = [(k,) for k in range(len(E)) if nz((k,), rels)]
    while fr:
        out += fr; nxt = []
        for p in fr:
            t = E[p[-1]][1]
            for b, e in enumerate(E):
                if e[0] == t and nz(p + (b,), rels): nxt.append(p + (b,))
        if nxt and len(nxt[0]) > maxlen: return None
        fr = nxt
    return out

def exact_gate(E, rels, paths, i):
    outs = [b for b, e in enumerate(E) if e[0] == i]
    if not outs: return False
    for p in paths:
        if E[p[-1]][1] == i and E[p[0]][0] != i and not any(nz(p + (b,), rels) for b in outs): return False
    return True

def hom_tilting(V, E, rels, i):
    ix = {v: t for t, v in enumerate(V)}
    A = c3.AlgG(len(V), [(ix[e[0]], ix[e[1]]) for e in E], [{r: 1} for r in rels], False)
    I = ix[i]; summ = [sc.stalk(u) if u != I else sc.mutate_cx(A, I, 'R') for u in range(A.n)]
    deg = {}
    for m in (-2, -1, 1, 2):
        d = sum(sc.cx_hom(A, summ[x], summ[y], m).dim for x in range(A.n) for y in range(A.n) if x == I or y == I)
        if d: deg[m] = d
    return not deg, deg

def position(alg, i):
    """'on': i on an oriented cycle; 'down': i reaches a cycle (D(i) cyclic); 'up': D(i) acyclic and a cycle vertex reaches i;
    'off': no oriented cycle reaches i or is reached from i.  Only 'up' and 'off' steps are carried out by `mutateAtVertex`
    (it raises CyclicDescendantsError when D(i) has a cycle, E-222)."""
    q = nx.DiGraph(alg.quiver)
    cyc = set().union(*[c for c in nx.strongly_connected_components(q) if len(c) > 1 or q.has_edge(next(iter(c)), next(iter(c)))] or [set()])
    if i in cyc: return 'on'
    if any(c in nx.descendants(q, i) for c in cyc): return 'down'      # D(i) cyclic: mutateAtVertex refuses
    return 'up' if any(c in nx.ancestors(q, i) for c in cyc) else 'off'

def compare(alg, tally, examples, tag, plantedRef=None):
    """all (side, i) on one algebra; returns number of rows"""
    rows = 0
    for side in ('A', 'Aop'):
        a = alg if side == 'A' else pathAlgebra.dualPathAlgebra(alg)
        m = model(a)
        if m is None: tally['nonmonomial side'] += 1; continue
        V, E, rels = m
        paths = nonzero_paths(E, rels)
        if paths is None: tally['infinite-dim side'] += 1; continue
        for i in V:
            if not any(e[0] == i for e in E) or a.quiver.has_edge(i, i): continue
            pos = position(a, i)
            g = mutation.mutationIsPossibleAtVertex(a, i); j = mutation.isTiltingAtVertex(a, i)
            x = exact_gate(E, rels, paths, i); h, deg = hom_tilting(V, E, rels, i)
            rows += 1
            k = (tag, pos)
            tally[k + ('rows',)] += 1
            tally[k + ('gate admits',)] += g
            tally[k + ('Hom tilting',)] += h
            tally[k + ('gate != Hom',)] += g != h
            tally[k + ('isTilting != Hom',)] += j != h
            tally[k + ('exact != Hom',)] += x != h
            tally[k + ('gate&isTilting both != Hom',)] += (g != h) and (j != h)
            tally[k + ('gate over',)] += g and not h; tally[k + ('isTilting over',)] += j and not h
            if (g != h or x != h or (j != h and sum(1 for e in examples if e[5] == e[11]) < 10)) and len(examples) < 60:
                examples.append((tag, pos, side, i, 'gate', g, 'isTilting', j, 'exact', x, 'Hom', h, deg, fmt(E, rels)))
    return rows

def fmt(E, rels):
    return 'arrows ' + ' '.join('%d:%s>%s' % (k, e[0], e[1]) for k, e in enumerate(E)) + ' | zero ' + ' '.join('.'.join(map(str, r)) for r in rels)

def build(n, arrows, zeros):
    """library PathAlgebra from vertex lists: arrows [(s,t)], zeros [vertex sequence]"""
    q = nx.MultiDiGraph(); q.add_nodes_from(range(1, n + 1))
    E = [(s, t, q.add_edge(s, t)) for s, t in arrows]
    rels = []
    for z in zeros:
        p = tuple(next(e for e in E if e[0] == z[k] and e[1] == z[k + 1]) for k in range(len(z) - 1)); rels.append({p: 1})
    return procedure.toPathAlgebra(q, rels)

# ---------------------------------------------------------------- controls
def controls():
    T = collections.Counter(); ex = []
    print('[PLANTED] cases where the library gate must disagree with Hom (non-simple nonzero path killed into i):')
    # E-329 bag: x1:1>2 x2:2>3 a:4>5 b:5>4 g:5>1; zero a.b and b.a.g.x1 ; i = 1 (refusing path b.a.g: 5>4>5>1, non-simple)
    bag = build(5, [(1, 2), (2, 3), (4, 5), (5, 4), (5, 1)], [[4, 5, 4], [5, 4, 5, 1, 2]])
    # isTilting docstring case: 1<->2 with 2-cycle killed at 1 (1.2.1 = 0), 1>3>4 with (1 3 4) = 0, at 3:
    #   non-simple path 2>1>2>1>3?  we let the engine decide; library claims True, the docstring says the answer is False
    doc = build(4, [(1, 2), (2, 1), (1, 3), (3, 4)], [[1, 2, 1], [1, 3, 4]])
    # [SIMPLE] seed of toolsmith_cycgate (x1:1>2 x2:2>1 x3:2>3 ...): planted here as 1<->2, 2>3>4, zero 2.1.2.3? rebuilt:
    simple = build(4, [(1, 2), (2, 1), (2, 3), (3, 4)], [[1, 2, 1], [2, 1, 2, 3, 4]])
    for nm, A in (('bag (E-329)', bag), ('isTiltingAtVertex docstring', doc), ('simple-seed variant', simple)):
        compare(A, T, ex, 'planted:' + nm)
    print('[ACYCLIC controls] every LNA at n = 5, every vertex, both sides (gate = Hom expected, J != 0 cases included):')
    for row in lna_rows(5):
        compare(nk.LinearNakayamaAlgebra(5, row), T, ex, 'LNA n=5')
    report(T, ex)

def report(T, ex):
    keys = sorted({k[:2] for k in T if isinstance(k, tuple)})
    cols = ('rows', 'gate admits', 'Hom tilting', 'gate != Hom', 'gate over', 'isTilting != Hom', 'isTilting over', 'exact != Hom', 'gate&isTilting both != Hom')
    print('| sample | cycle vs i | ' + ' | '.join(cols) + ' |')
    print('|' + '---|' * (len(cols) + 2))
    for k in keys: print('| %s | %s | ' % k + ' | '.join(str(T[k + (c,)]) for c in cols) + ' |')
    for k in T:
        if not isinstance(k, tuple): print('  ', k, T[k])
    for e in ex: print('  [DISAGREE]', e)

# ---------------------------------------------------------------- walk
def lna_rows(n):
    cand = [(s, l) for s in range(1, n - 1) for l in range(2, n - s + 1)]
    for r in range(0, n):
        for rl in itertools.combinations(cand, r):
            if len({s for s, _ in rl}) < r: continue
            if any(a != b and a[0] <= b[0] and b[0] + b[1] <= a[0] + a[1] for a in rl for b in rl): continue
            row = [0] * (n - 2)
            for s, l in rl: row[s - 1] = l
            yield row

def cyclic(alg): return not nx.is_directed_acyclic_graph(alg.quiver)

def walk(N, D, MAXN=10 ** 9, roots=None, guard=True, label='LNA'):
    t0 = time.time(); T = collections.Counter(); ex = []; seen_cyc = {}; errmsg = {}
    if roots is None: roots = [nk.LinearNakayamaAlgebra(N, row) for row in lna_rows(N)]
    for root in roots:
        key0 = search._coxeterKeyOrNone(root); row = None
        seen = {lnaMoves._stateKey(root): 0}; frontier = [(root, False)]
        for d in range(D):
            nxt = []
            for alg, deep in frontier:
                isc = cyclic(alg)
                for v in sorted(alg.quiver.nodes):
                    for signed in (v, -v):
                        tgt = alg if signed > 0 else pathAlgebra.dualPathAlgebra(alg)
                        if not mutation.mutationIsPossibleAtVertex(tgt, v): continue
                        try: ch = mutation.quiverMutationAtVertices(lnaMoves._copy(alg), [signed])
                        except Exception as e:
                            # swallowed here; lnaMoves.movesFrom / localMutationSequences would raise instead
                            T['mutation error'] += 1; T['mutation error %s' % type(e).__name__] += 1
                            T['mutation error on %s parent' % ('cyclic' if isc else 'acyclic')] += 1
                            if type(e).__name__ not in errmsg: errmsg[type(e).__name__] = str(e)[:160]
                            continue
                        if ch is None: continue
                        chc = cyclic(ch)
                        if not chc and guard:
                            k = search._coxeterKeyOrNone(ch)
                            if k is not None and k != key0: continue
                        s = lnaMoves._stateKey(ch)
                        if s in seen: continue
                        seen[s] = d + 1
                        # LIB reaches ch iff every node on the way is acyclic (ch may itself be cyclic: a leaf)
                        libreach = not isc and not deep
                        if chc:
                            prev = seen_cyc.get(s)
                            if prev is None: seen_cyc[s] = (ch, libreach, row, d + 1)
                            elif libreach and not prev[1]: seen_cyc[s] = (ch, True, row, d + 1)
                        nxt.append((ch, deep or isc))
            frontier = nxt
            if len(seen) > MAXN: T['capped roots'] += 1; break
        T['nodes'] += len(seen)
    T['roots'] = len(roots)
    T['cyclic nodes'] = len(seen_cyc)
    print(label + (' guard' if guard else ' NO guard') + ' walk n=%d depth %d: %d roots, %d nodes, %d distinct cyclic nodes (%.0f s); comparing ...' % (N, D, len(roots), T['nodes'], len(seen_cyc), time.time() - t0), flush=True)
    for s, (ch, lib, row, d) in seen_cyc.items():
        T['cyclic LIB-reachable (leaf)'] += lib
        if model(ch) is None and model(pathAlgebra.dualPathAlgebra(ch)) is None: T['cyclic non-monomial'] += 1; continue
        T['cyclic monomial'] += 1
        compare(ch, T, ex, label + (' LIB leaf' if lib else ' LNAM only'))
    print('done %.0f s' % (time.time() - t0))
    for k, m in errmsg.items(): print('  [MUTATION ERROR sample] %s: %s' % (k, m))
    report(T, ex)

def random_root(n, rng):
    """random connected acyclic quiver on 1..n (arrows a < b, no parallel), monomial zero relations of length 2-3"""
    while True:
        arrows = sorted({(a, b) for a, b in (sorted(rng.sample(range(1, n + 1), 2)) for _ in range(rng.randint(n - 1, n + 2)))})
        g = nx.Graph(arrows); g.add_nodes_from(range(1, n + 1))
        if not nx.is_connected(g): continue
        succ = collections.defaultdict(list)
        for a, b in arrows: succ[a].append(b)
        paths = [[a, b, c] for a, b in arrows for c in succ[b]] + [[a, b, c, d] for a, b in arrows for c in succ[b] for d in succ[c]]
        rng.shuffle(paths); zeros = []
        for p in paths[:rng.randint(0, 3)]:
            if not any(len(z) <= len(p) and any(p[o:o + len(z)] == z for o in range(len(p) - len(z) + 1)) for z in zeros) and \
               not any(len(p) < len(z) and any(z[o:o + len(p)] == p for o in range(len(z) - len(p) + 1)) for z in zeros):
                zeros.append(p)
        return build(n, arrows, zeros)

def randwalk(n, count, D, seed):
    import random
    rng = random.Random(seed)
    walk(n, D, roots=[random_root(n, rng) for _ in range(count)], guard=True, label='RAND n=%d seed %d' % (n, seed))

def libcount(N, D):
    """instrumented library walk: does `mutationSearchDepthFirst` ever call the gate on a cyclic parent?"""
    T = collections.Counter(); orig = mutation.mutationIsPossibleAtVertex
    def spy(a, v, **kw):
        T['gate calls'] += 1; T['gate calls on cyclic parent'] += cyclic(a); return orig(a, v, **kw)
    mutation.mutationIsPossibleAtVertex = spy
    def visit(a, path): T['visits'] += 1; T['cyclic visits'] += cyclic(a)
    t0 = time.time()
    for row in lna_rows(N):
        r = nk.LinearNakayamaAlgebra(N, row)
        for start in (r, pathAlgebra.dualPathAlgebra(r)):
            search.mutationSearchDepthFirst(start, D, printOutput=False, visitor=visit)
    mutation.mutationIsPossibleAtVertex = orig
    print('LIBCOUNT n=%d depth %d: %s  %.0f s' % (N, D, dict(T), time.time() - t0))

if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'controls': controls()
    elif m == 'walk': walk(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9)
    elif m == 'nog': walk(int(sys.argv[2]), int(sys.argv[3]), guard=False)
    elif m == 'rand': randwalk(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
    elif m == 'libcount': libcount(int(sys.argv[2]), int(sys.argv[3]))
