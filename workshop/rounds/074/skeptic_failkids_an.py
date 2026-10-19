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
print(' per class:', {c: len({k for k, e in zip(keys, allE) if e[0] == c and k is not None}) for c in (1, 2)})
items = [(i, e[5]) for i, e in enumerate(allE)]
cl, gs = classes_of(items)
print('distinct children up to vertex relabelling + arrow relabelling (sign-blind structure iso):', len(cl))
for c in (1, 2):
    sub = [(i, e[5]) for i, e in enumerate(allE) if e[0] == c]
    print('  class', c, 'edges', len(sub), 'iso classes', len(classes_of(sub)[0]))
print('cross-class iso merges:', sum(1 for r in cl if len({allE[i][0] for i in r}) > 1))
print('iso classes with >1 distinct canonicalKey (relabelled copies):', sum(1 for r in cl if len({keys[i] for i in r if keys[i] is not None}) > 1))
print('iso-class sizes (edges):', sorted((len(r) for r in cl), reverse=True))
print('children with a multi-term relation (sign blind caveat):', sum(1 for e in allE if mixed(e[5])), 'of', len(allE))
print('child dims: arrows', collections.Counter(e[5].quiver.number_of_edges() for e in allE), 'parallel-arrow children', sum(1 for e in allE if len(fingerprint.arrowBundles(e[5].quiver)) > 0))
# E-151 membership: the E-151/E-154/E-157 rebuild is the same key-guarded BFS cut at 20000 expansions (rounds/049/toolsmith_collect.py).
# A failing edge is in it iff its parent was expanded before expansion 20000: prefix levels (d < 8) always (<= 12033 / 16784 expansions); a level-8 parent iff
# prefix + its position in the frontier < 20000.  Positions come from the checkpoint written by skeptic_failkids.py.
CAP = 20000; PRE = {1: 12033, 2: 16784}
inE151 = []
for i, e in enumerate(allE):
    c, d, pk = e[0], e[1], e[2]
    if d < 8: inE151.append(True); continue
    t0_, f0_, cur_, seen_ = pickle.load(open('%s/c%d_s0.pkl.frontier' % (D, c), 'rb')) if ('fr', c) not in globals() else globals()[('fr', c)]
    globals()[('fr', c)] = (t0_, f0_, cur_, seen_)
    ks = [fingerprint.canonicalKey(a) for a in cur_]
    pos = [j for j, k in enumerate(ks) if k is not None and k == pk]
    inE151.append(bool(pos) and PRE[c] + pos[0] < CAP)
    if not pos: print('WARNING parent key not located (None key?)', c, d)
for c in (1, 2):
    idx = [i for i, e in enumerate(allE) if e[0] == c and inE151[i]]
    print('class', c, 'edges inside the 20000-expansion window:', len(idx), '(E-157: %d)' % {1: 16, 2: 9}[c], '; distinct (parent key, v):', len({(allE[i][2], allE[i][3]) for i in idx}))
ref = [i for i in range(len(allE)) if inE151[i]]
out = [i for i in range(len(allE)) if not inE151[i]]
rk = {keys[i] for i in ref if keys[i] is not None}
print('E-151-window edges', len(ref), 'distinct child keys', len(rk), 'None', sum(keys[i] is None for i in ref))
rcl, _ = classes_of([(i, allE[i][5]) for i in ref]); print('E-151-window children, structure iso classes:', len(rcl))
# cross-check class 2 against the independent rebuild by rounds/049/toolsmith_collect.py
recs = pickle.load(open('%s/e151_c2.pkl' % D, 'rb')); recs = [r for r in recs if r.get('kind') == 'fail']
ek = {fingerprint.canonicalKey(r['childObj']) for r in recs}
print('collect.py class 2: failing records', len(recs), 'distinct child keys', len(ek), '; equal to my window keys for class 2:', ek == {keys[i] for i in ref if allE[i][0] == 2})
# overlap by iso
allc = [(0, allE[i][5]) for i in range(len(allE))]
acl, _ = classes_of(allc)
refset = set(ref)
both = [r for r in acl if any(j in refset for j in r)]
onlyN = [r for r in acl if not any(j in refset for j in r)]
print('all-70 iso classes', len(acl), '; containing an E-151-window edge:', len(both), '; with NO window edge (new):', len(onlyN))
print('by canonicalKey: keys outside window keys:', len({keys[i] for i in out if keys[i] is not None and keys[i] not in rk}), '; new edges with None key:', sum(keys[i] is None for i in out))
for r in onlyN:
    e = allE[r[0]]
    print('  NEW iso class: class', e[0], 'edges', len(r), 'parent depth', sorted({allE[j][1] for j in r}), 'child arrows', e[5].quiver.number_of_edges(), 'rels', len(procedure.relationsFrom(e[5])), 'canonicalKey None' if keys[r[0]] is None else '')
pickle.dump(dict(allE=[(e[0], e[1], e[3]) for e in allE], inE151=inE151, onlyN=onlyN, acl=acl), open(D + '/an_summary.pkl', 'wb'))

# ---- sign-aware refinement: within each sign-blind class, is every member gauge-isomorphic to the first (exists a structure iso whose relabelled relations
# have the same gauge-canonical form)?  Uses fingerprint._relationsUnder (the canonicalKey internals) with an explicit full arrow renaming.
from quivermutation import arrowPaths
def signedIso(A, B, limit=3000):
    gA, gB = struct(A), struct(B)
    gm = nx.isomorphism.DiGraphMatcher(gB, gA, node_match=lambda a, b: a['lab'] == b['lab'], edge_match=lambda a, b: a['lab'] == b['lab'])
    arrowsA = arrowPaths.arrowsOf(A.quiver); relA = procedure.relationsFrom(A); relB = procedure.relationsFrom(B)
    target = fingerprint._relationsUnder(relA, {}, arrowsA, gauge=True)
    for t, m in enumerate(gm.isomorphisms_iter()):
        if t >= limit: return None
        ren = {k[1:]: v[1:] for k, v in m.items() if k[0] == 'a'}   # B arrow (x,y,k) -> A arrow
        if fingerprint._relationsUnder(relB, ren, arrowsA, gauge=True) == target: return True
    return False
bad = 0; undecided = 0
for r in acl:
    for j in r[1:]:
        res_ = signedIso(allE[r[0]][5], allE[j][5])
        if res_ is None: undecided += 1
        elif not res_: bad += 1
print('sign-aware check: members not gauge-isomorphic to the first of their class:', bad, '; undecided (iso enumeration cap):', undecided)
# details of the representative of each new class
def show(a):
    return 'arrows ' + str(sorted(a.quiver.edges(keys=True))) + '\n      rels ' + str([{tuple(x for x in p): c for p, c in rel.items()} for rel in procedure.relationsFrom(a)])
for r in onlyN:
    j = min(r, key=lambda j: allE[j][5].quiver.number_of_edges() + allE[j][4].quiver.number_of_edges())
    e = allE[j]
    print('NEW class rep: cls', e[0], 'v', e[3], '\n   parent', show(e[4]), '\n   child ', show(e[5]))
newedges = [j for r in onlyN for j in r]
print('new-class edges', len(newedges), 'distinct canonicalKeys among them', len({keys[j] for j in newedges}), '; keys outside window keys that belong to joined classes:', len({keys[j] for j in out if j not in newedges and keys[j] not in rk}))
