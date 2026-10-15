"""Round 067 (toolsmith), answers to referee: (1) end check, (2) non-vacuity control, for the n = 10 cycle scan.
(1) each replayed end is compared with the LNA built from the `end relLengths` printed in 054/experimentalist_amerge_out.txt:
    quiver edges equal, and relation structure graphs isomorphic (labelled-equal also reported).
(2) the same loop is run with a 3-cycle injected into one parent (an extra arrow n -> 1 on a COPY at step 3 of each path): cyclic must be > 0;
    and on a bare 3-cycle quiver with a relation.
Usage (repo root): .venv/bin/python workshop/rounds/067/toolsmith_cycscan2.py"""
import sys, copy
import networkx as nx
sys.path.insert(0, '.')
from quivermutation import nakayama as nk, mutation, reduction, pathAlgebra, procedure, arrowPaths as ap
START = [0, 5, 0, 4, 0, 3, 3, 0]
PATHS = [[2, 1, 6, 5, -4, -10, -8], [-10, -8, -4, -5, 2, 1, 6], [4, 3, 2, 1, 7, 6, 9], [4, 3, 7, 2, 1, 6, 9], [4, 3, 7, 9, 2, 1, 6]]
ENDREL = [[4, 5, 0, 5, 5, 0, 0, 0], [5, 5, 5, 0, 4, 4, 0, 0], [6, 0, 5, 0, 4, 0, 3, 0], [6, 0, 5, 0, 4, 0, 3, 0], [6, 0, 5, 0, 4, 0, 3, 0]]  # copied from the out file

def struct(a):
    edges = sorted(set(a.quiver.edges()))
    g = nx.DiGraph()
    for v in a.quiver.nodes(): g.add_node(('v', v), lab='v')
    for x, y in edges: g.add_edge(('v', x), ('v', y), lab='a')
    for i, rel in enumerate(procedure.relationsFrom(a)):
        g.add_node(('r', i), lab='r')
        for j, term in enumerate(rel):
            g.add_node(('p', i, j), lab='p%d' % len(term)); g.add_edge(('r', i), ('p', i, j), lab='rp')
            for k, arr in enumerate(term):
                g.add_edge(('p', i, j), ('v', arr[0]), lab='s%d' % k)
            g.add_edge(('p', i, j), ('v', term[-1][1]), lab='e')
    return g
def same(g, h):
    return nx.is_isomorphic(g, h, node_match=lambda x, y: x['lab'] == y['lab'], edge_match=lambda x, y: x['lab'] == y['lab'])
def labelled_equal(g, h):
    return set(g.nodes(data='lab')) == set(h.nodes(data='lab')) and set(g.edges(data='lab')) == set(h.edges(data='lab'))

def scan(inject=False, check_ends=False):
    tot = dict(paths=0, parents=0, cyclic=0, ends=0, endcyclic=0, ends_checked=0, ends_edges_equal=0, ends_rel_iso=0, ends_labelled_equal=0)
    for pi, p in enumerate(PATHS):
        alg = nk.LinearNakayamaAlgebra(10, list(START)); tot['paths'] += 1
        for i, s in enumerate(p):
            v = abs(s); a = pathAlgebra.dualPathAlgebra(alg) if s < 0 else alg
            q = a.quiver
            if inject and i == 3:  # copy only, never mutated
                q = nx.MultiDiGraph(q); q.add_edge(10, 1)
            tot['parents'] += 1; tot['cyclic'] += ap.hasOrientedCycle(q)
            a2 = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(a), v))
            alg = pathAlgebra.dualPathAlgebra(a2) if s < 0 else a2
        tot['ends'] += 1; tot['endcyclic'] += ap.hasOrientedCycle(alg.quiver)
        if check_ends:
            tgt = nk.LinearNakayamaAlgebra(10, ENDREL[pi] + [0] * 0) if len(ENDREL[pi]) == 8 else None
            g, h = struct(alg), struct(tgt)
            tot['ends_checked'] += 1
            tot['ends_edges_equal'] += set(alg.quiver.edges()) == set(tgt.quiver.edges())
            tot['ends_rel_iso'] += same(g, h); tot['ends_labelled_equal'] += labelled_equal(g, h)
            print('path', pi + 1, 'end edges', sorted(set(alg.quiver.edges())) == sorted(set(tgt.quiver.edges())), 'rel-iso', same(g, h), 'labelled-equal', labelled_equal(g, h),
                  'nrel', len(procedure.relationsFrom(alg)), 'vs', len(procedure.relationsFrom(tgt)))
    return tot
print('scan', scan(check_ends=True))
print('control (3-cycle injected into the step-4 parent of each path, scan otherwise identical)', scan(inject=True))
q = nx.MultiDiGraph([(1, 2), (2, 3), (3, 1)])
print('control (bare 3-cycle quiver) hasOrientedCycle =', ap.hasOrientedCycle(q))
