"""Round 054, T10 (ii): witness path for group-A n=10 merge 05040330 -> 33460000.
Meet in the middle, matching up to relabelling (WL invariant + exact isomorphism of the quiver-with-relations).
Stages:  fwd D   : reach set of 05040330 member [0,5,0,4,0,3,3,0] to depth D (with dual), pickled
         back D  : reach sets of the 84 members of orbit 33460000 to depth D (pool of 4), pickled
         join    : join, verify isomorphism, build full path, replay with gate / tiltingPlus / key
Usage (repo root): timeout 10m .venv/bin/python workshop/rounds/054/experimentalist_amerge.py fwd 4 OUT.pkl"""
import sys, json, copy, time, pickle, multiprocessing
import networkx as nx
sys.path.insert(0, '.')
from quivermutation import nakayama as nk, mutation, procedure, reduction, search, pathAlgebra, arrowPaths, lnaMoves
START = [0, 5, 0, 4, 0, 3, 3, 0]
CK = 'workshop/rounds/051/experimentalist_merges10.jsonl'

def struct(edges, rels):
    """Structure graph: vertex nodes, arrow edges, path chains for relations."""
    g = nx.DiGraph()
    for v in set(x for e in edges for x in e): g.add_node(('v', v), lab='v')
    for a, b in edges: g.add_edge(('v', a), ('v', b), lab='a')
    for i, rel in enumerate(rels):
        g.add_node(('r', i), lab='r%d' % len(rel))
        for j, path in enumerate(rel):
            g.add_node(('p', i, j), lab='p%d' % len(path)); g.add_edge(('r', i), ('p', i, j), lab='rp')
            for k, v in enumerate(path):
                g.add_node(('q', i, j, k), lab='q'); g.add_edge(('p', i, j), ('q', i, j, k), lab='pq%d' % k)
                g.add_edge(('q', i, j, k), ('v', v), lab='qv')
    return g

def inv(edges, rels):
    return nx.weisfeiler_lehman_graph_hash(struct(edges, rels), node_attr='lab', edge_attr='lab', iterations=4)

def reach(start, depth):
    """{(edges, rels): path} over all quivers reached (right and, via dual, left steps), skipping parallel arrows."""
    found = {}
    def visit(q, mv, dual=False):
        if dual: q = pathAlgebra.dualPathAlgebra(q)
        k = search.quiverKey(q)
        if k is None: return
        p = [(-s if dual else s) for s in mv]
        if k not in found or len(p) < len(found[k]): found[k] = p
    alg = nk.LinearNakayamaAlgebra(10, list(start))
    search.mutationSearchDepthFirst(copy.deepcopy(alg), depth, [], 'm', printOutput=False, visitor=visit)
    search.mutationSearchDepthFirst(pathAlgebra.dualPathAlgebra(alg), depth, [], 'm', printOutput=False,
                                    visitor=lambda q, mv: visit(q, mv, True))
    return found

def reach_t(a): return tuple(a[0]), reach(*a)

if __name__ == '__main__':
    stage, D, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    t0 = time.time()
    if stage == 'fwd':
        f = reach(START, D); pickle.dump(f, open(out, 'wb'))
        print('fwd depth', D, 'quivers', len(f), 'seconds %.0f' % (time.time() - t0))
    elif stage == 'back':
        tg = sorted({tuple(json.loads(l)['member']) for l in open(CK) if json.loads(l)['orbit'] == '33460000'})
        lim = int(sys.argv[4]) if len(sys.argv) > 4 else len(tg)
        with multiprocessing.Pool(4) as pool:
            res = dict(pool.map(reach_t, [(m, D) for m in tg[:lim]], chunksize=1))
        pickle.dump(res, open(out, 'wb'))
        print('back depth', D, 'targets', len(res), 'quivers/target', [len(v) for v in res.values()][:5], 'seconds %.0f' % (time.time() - t0))

# ---------------------------------------------------------------- join and replay
def graphOf(key): return struct(list(key[0]), [list(map(list, r)) for r in key[1]])

def iso(kx, ky):
    """vertex map Y -> X if the labelled quivers-with-relations are isomorphic, else None."""
    gm = nx.isomorphism.DiGraphMatcher(graphOf(ky), graphOf(kx),
                                       node_match=lambda a, b: a['lab'] == b['lab'], edge_match=lambda a, b: a['lab'] == b['lab'])
    if not gm.is_isomorphic(): return None
    return {k[1]: v[1] for k, v in gm.mapping.items() if k[0] == 'v'}

def join(F, B, maxTotal=99):
    fi = {}
    for k, p in F.items(): fi.setdefault(inv(*k), []).append((k, p))
    hits = []
    for m, R in B.items():
        for ky, pb in R.items():
            for kx, pf in fi.get(inv(*ky), []):
                if len(pf) + len(pb) > maxTotal: continue
                phi = iso(kx, ky)
                if phi is None: continue
                back = [(-s if s > 0 else -s) for s in []]  # placeholder, replaced below
                inv_b = [(-1 if s > 0 else 1) * phi[abs(s)] for s in reversed(pb)]
                hits.append((len(pf) + len(pb), m, pf, pb, inv_b, kx, ky, phi))
    hits.sort(key=lambda h: h[0]); return hits

def replay(path, tiltingPlus):
    alg = nk.LinearNakayamaAlgebra(10, list(START)); out = []
    base = search._coxeterKeyOrNone(alg)
    for s in path:
        v = abs(s); dual = s < 0
        a = pathAlgebra.dualPathAlgebra(alg) if dual else alg
        gate = mutation.mutationIsPossibleAtVertex(a, v)
        J0 = bool(tiltingPlus(a.quiver, procedure.relationsFrom(a), v))
        a2 = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(a), v))
        alg = pathAlgebra.dualPathAlgebra(a2) if dual else a2
        out.append((s, gate, J0, search._coxeterKeyOrNone(alg) in (None, base)))
    return alg, out

def loadTP():
    ns = {}
    exec(compile(open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0], 'h015', 'exec'), ns)
    return ns['tiltingPlus']
