"""Round 134 (toolsmith), response to referee item 3: are the two refused children at 303500030 one algebra?

Rebuilds the default depth-6 F-only ball at ROOT exactly as toolsmith_downstream.py step 1 does, collects the refused
(gate-admitted, key-kept, non-tilting) children, and for every pair tests isomorphism two ways:
  (a) key test: for each vertex bijection that is an isomorphism of the underlying quivers (arrow multiplicities kept),
      relabel X by it (shapeKeys.relabel) and compare fingerprint.canonicalKey with Y (the key quotients parallel-arrow
      renumbering and the sign gauge, E-164);
  (b) ideal test: some quiver isomorphism carries one ideal onto the other (arrowPaths.isInIdeal both ways; generating
      sets may differ), as rounds/129 maverick_inv.sameIdeal.
Also prints cheap invariants (vertex in/out degrees, #relations, SNF/signature) and, as a CONTROL, runs (a) on a child
against a random relabelling of itself (must say True) and on a refused child against a step-6 tilting child (must say False).
usage: timeout 10m .venv/bin/python -W ignore workshop/rounds/134/toolsmith_iso.py 303500030 root
"""
import sys, json, copy, random, itertools, collections
sys.path.insert(0, '.')
import importlib.util
import networkx as nx
spec = importlib.util.spec_from_file_location('dcount', 'workshop/rounds/133/toolsmith_dcount.py')
dc = importlib.util.module_from_spec(spec); spec.loader.exec_module(dc)
from quivermutation import search, mutation, fingerprint, procedure, arrowPaths, shapeKeys

def simple(a):
    g = nx.DiGraph(); g.add_nodes_from(a.quiver.nodes)
    for x, y in set(a.quiver.edges()): g.add_edge(x, y, m=a.quiver.number_of_edges(x, y))
    return g

def vertexIsos(X, Y):
    gm = nx.isomorphism.DiGraphMatcher(simple(X), simple(Y), edge_match=lambda p, q: p['m'] == q['m'])
    return gm.isomorphisms_iter()

def keyIso(X, Y, limit=200000):
    kY = fingerprint.canonicalKey(Y)
    for t, m in enumerate(vertexIsos(X, Y)):
        if t >= limit: return None
        if fingerprint.canonicalKey(shapeKeys.relabel(X, m)) == kY: return m
    return False

def idealIso(X, Y, limit=20000):
    rX, rY = procedure.relationsFrom(X), procedure.relationsFrom(Y)
    for t, m in enumerate(vertexIsos(X, Y)):
        if t >= limit: return None
        Xm = shapeKeys.relabel(X, m)
        # parallel bundles: try every renumbering of keys within each bundle
        bund = collections.defaultdict(list)
        for a in arrowPaths.arrowsOf(Xm.quiver): bund[a[:2]].append(a)
        ends = sorted(bund)
        for perm in itertools.product(*[itertools.permutations(bund[e]) for e in ends]):
            ren = {a: b for e, p in zip(ends, perm) for a, b in zip(sorted(bund[e]), p)}
            rXm = [{tuple(ren[a] for a in path): c for path, c in rel.items()} for rel in procedure.relationsFrom(Xm)]
            if all(arrowPaths.isInIdeal(Y.quiver, rY, r) for r in rXm if r) and all(arrowPaths.isInIdeal(Y.quiver, rXm, r) for r in rY if r):
                return m
    return False

def degs(a): return sorted((a.quiver.in_degree(v), a.quiver.out_degree(v)) for v in a.quiver.nodes)

def main(root, side):
    N, D = 11, 6
    s = dc.startOf(N, root, side); base = search._coxeterKeyOrNone(s)
    gate = mutation.mutationIsPossibleAtVertex; depthOf = {}; refused = []; tilt6 = []
    def visit(A, path): depthOf[id(A)] = len(path)
    def spy(alg, v):
        ok = gate(alg, v)
        if ok and id(alg) in depthOf:
            t = dc._tilt(alg, v)
            if (not t) or (depthOf[id(alg)] + 1 == D and len(tilt6) < 3):
                legal, c = dc.childOf(alg, v)
                if legal and search._coxeterKeyOrNone(c) == base: (tilt6 if t else refused).append((depthOf[id(alg)] + 1, v, c))
        return ok
    mutation.mutationIsPossibleAtVertex = spy
    search.mutationSearchDepthFirst(copy.deepcopy(s), D, printOutput=False, visitor=visit, visited=fingerprint.Visited(), baseKey=base)
    mutation.mutationIsPossibleAtVertex = gate
    uniq = {}
    for st, v, c in refused: uniq.setdefault(fingerprint.canonicalKey(c), (st, v, c))
    R = list(uniq.values())
    for i, (st, v, c) in enumerate(R):
        print(json.dumps(dict(child=i, step=st, v=v, vertices=len(c.quiver.nodes), arrows=c.quiver.number_of_edges(),
                              relations=len(procedure.relationsFrom(c)), degs=degs(c), inv=str(dc.inv(c)))), flush=True)
    for i, j in itertools.combinations(range(len(R)), 2):
        X, Y = R[i][2], R[j][2]
        k = keyIso(X, Y); d = idealIso(X, Y)
        print('PAIR %d,%d: quiver isos %d; key-iso %s; ideal-iso %s' % (i, j, sum(1 for _ in vertexIsos(X, Y)),
              'none' if k is False else ('cap' if k is None else 'YES %s' % k), 'none' if d is False else ('cap' if d is None else 'YES %s' % d)), flush=True)
    # controls
    X = R[0][2]; vs = list(X.quiver.nodes); p = vs[:]; random.Random(134).shuffle(p); Xp = shapeKeys.relabel(X, dict(zip(vs, p)))
    print('CONTROL relabelled self: labelled keys equal %s; key-iso %s; ideal-iso %s' % (fingerprint.canonicalKey(Xp) == fingerprint.canonicalKey(X),
          keyIso(Xp, X) is not False, idealIso(Xp, X) is not False), flush=True)
    if tilt6:
        T = tilt6[0][2]
        print('CONTROL refused child vs tilting step-6 child (in-class, different SNF): key-iso %s; ideal-iso %s'
              % (keyIso(X, T) is not False, idealIso(X, T) is not False), flush=True)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'root')
