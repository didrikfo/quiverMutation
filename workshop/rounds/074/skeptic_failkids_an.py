"""Round 074 (skeptic), T10: distinct children of the 70 key-kept failing edges (E-204 balls), by canonicalKey and by a vertex-relabelling-aware
isomorphism of the structure graph, and overlap with the E-151 children (E-157/E-163/E-212 joined set).
  skeptic_failkids_an.py DIR   DIR holds c{1,2}_s{0..3}.pkl from skeptic_failkids.py and e151_c{1,2}.pkl from rounds/049/toolsmith_collect.py (7 c 20000 ...)
Structure graph: vertex nodes; one node per ARROW (parallel arrows distinct) with src/tgt edges; relation nodes -> term nodes (label |coef|, number of arrows)
-> arrow nodes with position labels.  Coefficient SIGNS are ignored (necessary test; reported: how many children have a relation with a mixed-sign multi-term)."""
import sys, pickle, glob, collections
sys.path.insert(0, '.')
import networkx as nx
from quivermutation import procedure, fingerprint
D = sys.argv[1]
def struct(a):
    g = nx.DiGraph()
    for v in a.quiver.nodes(): g.add_node(('v', v), lab='v')
    for (x, y, k) in a.quiver.edges(keys=True):
        g.add_node(('a', x, y, k), lab='a'); g.add_edge(('a', x, y, k), ('v', x), lab='s'); g.add_edge(('a', x, y, k), ('v', y), lab='t')
    for i, rel in enumerate(procedure.relationsFrom(a)):
        g.add_node(('r', i), lab='r')
        for j, (path, c) in enumerate(rel.items()):
            g.add_node(('p', i, j), lab='p%s_%d' % (abs(c), len(path))); g.add_edge(('r', i), ('p', i, j), lab='rp')
            for q, arr in enumerate(path): g.add_edge(('p', i, j), ('a',) + tuple(arr), lab='pos%d' % q)
    return g
def iso(g, h):
    return nx.is_isomorphic(g, h, node_match=lambda a, b: a['lab'] == b['lab'], edge_match=lambda a, b: a['lab'] == b['lab'])
def invariant(g):
    return nx.weisfeiler_lehman_graph_hash(g, node_attr='lab', edge_attr='lab')
def classes_of(items):
    """items: list of (tag, alg).  returns list of classes [ [indices...] ] by exact structure-graph iso"""
    buckets = collections.defaultdict(list); gs = [struct(a) for _, a in items]
    for i, g in enumerate(gs): buckets[invariant(g)].append(i)
    out = []
    for b in buckets.values():
        reps = []
        for i in b:
            for r in reps:
                if iso(gs[r[0]], gs[i]): r.append(i); break
            else: reps.append([i])
        out += reps
    return out, gs
def mixed(a):
    return sum(1 for rel in procedure.relationsFrom(a) if len(rel) > 1)
res = {}
for c in (1, 2):
    fails = []; T = collections.Counter(); stops = []
    for s in range(4):
        t, f, stop = pickle.load(open('%s/c%d_s%d.pkl' % (D, c, s), 'rb')); stops.append(stop)
        if s == 0:  # shard 0 carries the failures of levels before the last
            pass
        fails += f
        for k, v in dict(t).items(): T[k] += v
    # dedupe edges: (parent key, v); None keys cannot dedupe -> keep (reported)
    seen = set(); edges = []
    for (d, pk, v, pa, ch) in fails:
        if pk is not None:
            if (pk, v) in seen: continue
            seen.add((pk, v))
        edges.append((d, pk, v, pa, ch))
    print('class', c, 'shards stopped by budget', stops, 'raw fail records', len(fails), 'distinct (parent key,v) edges', len(edges), 'parent keys None', sum(1 for e in edges if e[1] is None))
    res[c] = edges
allE = [(c,) + e for c in (1, 2) for e in res[c]]
print('total failing edges', len(allE), 'by class', {c: len(res[c]) for c in res})
keys = [fingerprint.canonicalKey(e[5]) for e in allE]
print('child canonicalKey None:', sum(k is None for k in keys))
kc = collections.Counter(k for k in keys if k is not None)
print('distinct child canonicalKeys (non-None)', len(kc), '+ None children', sum(k is None for k in keys))
print(' per class:', {c: len({k for k, e in zip(keys, allE) if e[0] == c and k is not None}) for c in (1, 2)}))
