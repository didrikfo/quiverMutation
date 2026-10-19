"""Round 075 (toolsmith), T10 item 1: the class-1 child of E-213 missed the LABELLED-key join at 6 + 5 (toolsmith_join.py).  Same F/R J = 0 moves, but
nodes are met up to ISOMORPHISM: WL hash of the structure graph (vertices, arrows, relation terms; rounds/073 struct) then explicit networkx iso.
  toolsmith_isojoin.py REP DT DC CAP BUDGET_SEC     (repo root; REP as in toolsmith_join.py).  Prints a HIT with total depth or a MISS with level sizes / cap / budget."""
import sys, re, ast, time
REP, DT, DC, CAP = (int(x) for x in sys.argv[1:5]); BUD = float(sys.argv[5])
src = open('workshop/rounds/075/toolsmith_join.py').read().split("# sanity:")[0]
sys.argv = ['x', str(REP), str(DT), str(DC), str(CAP)]
exec(compile(src, 'join', 'exec'))
from quivermutation import arrowPaths
def struct(a):   # copied from rounds/074/skeptic_failkids_an.py: arrows are nodes (parallel arrows kept distinct), relation terms labelled |coef|, length
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
def H(a): return nx.weisfeiler_lehman_graph_hash(struct(a), node_attr='lab', edge_attr='lab', iterations=4)
def iso(A, B, limit=3000):
    """sign-AWARE: some structure isomorphism whose relabelled relations have the same gauge-canonical form (as rounds/074 signedIso)."""
    gA, gB = struct(A), struct(B)
    gm = nx.isomorphism.DiGraphMatcher(gB, gA, node_match=lambda a, b: a['lab'] == b['lab'], edge_match=lambda a, b: a['lab'] == b['lab'])
    arrowsA = arrowPaths.arrowsOf(A.quiver); relA = procedure.relationsFrom(A); relB = procedure.relationsFrom(B)
    target = fingerprint._relationsUnder(relA, {}, arrowsA, gauge=True)
    for t, m in enumerate(gm.isomorphisms_iter()):
        if t >= limit: return False
        ren = {k[1:]: v[1:] for k, v in m.items() if k[0] == 'a'}
        if fingerprint._relationsUnder(relB, ren, arrowsA, gauge=True) == target: return True
    return False
def ballH(starts, depth, tag, t0):
    seen = {}; objs = {}; frontier = []
    for s in starts:
        k = fingerprint.canonicalKey(s)
        if k in seen: continue
        seen[k] = 0; objs.setdefault(H(s), []).append((s, 0)); frontier.append(s)
    lv = [len(frontier)]
    for d in range(1, depth + 1):
        nxt = []
        for alg in frontier:
            if time.time() - t0 > BUD: print(tag, 'BUDGET stop at level', d, lv, flush=True); return seen, objs, lv, False
            for kind, v, c in moves(alg):
                k = fingerprint.canonicalKey(c)
                if k is None or k in seen: continue
                seen[k] = d; objs.setdefault(H(c), []).append((c, d)); nxt.append(c)
        lv.append(len(nxt)); frontier = nxt
    return seen, objs, lv, True
T0 = time.time()
tseen, tobjs, tlv, tok = ballH(classes[base], DT, 'target', T0)
print('target ball (labelled keys) levels', tlv, 'hash buckets', len(tobjs), 'complete', tok, '%.0fs' % (time.time() - T0), flush=True)
for name, alg in (('child', C), ('parent', P)):
    t1 = time.time(); k0 = fingerprint.canonicalKey(alg); seen = {k0}; frontier = [alg]; lv = [1]; best = None; state = 'uncapped'
    def check(a, d):
        global best
        for (o, td) in tobjs.get(H(a), []):
            if iso(a, o) and (best is None or d + td < best[0]): best = (d + td, d, td)
    check(alg, 0)
    for d in range(1, DC + 1):
        nxt = []
        for a in frontier:
            if len(seen) >= CAP: state = 'CAPPED'; break
            if time.time() - t1 > BUD: state = 'BUDGET'; break
            for kind, v, c in moves(a):
                k = fingerprint.canonicalKey(c)
                if k is None or k in seen: continue
                seen.add(k); nxt.append(c); check(c, d)
        lv.append(len(nxt)); frontier = nxt
        if best or state != 'uncapped': break
    print(name, ('ISO-HIT total %d (alg side %d + target side %d)' % best) if best else 'MISS', 'levels', lv, state, '%.0fs' % (time.time() - t1), flush=True)
