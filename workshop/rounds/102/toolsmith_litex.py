"""Round 102 (toolsmith): run the worked examples of arXiv:2112.08129 (Ex. 3.1, 3.2, Sect. 6) and arXiv:1504.02617 (Ex. 2.3/4.1, 4.2)
through quivermutation: gate (mutationIsPossibleAtVertex), J (perI, dim ker g_i; J = 0 iff empty), Coxeter-key guard, and the printed result.
usage (repository root):  .venv/bin/python workshop/rounds/102/toolsmith_litex.py"""
import sys, io, contextlib, itertools
import sympy, networkx as nx
sys.path.insert(0, '.'); sys.path.insert(0, 'tests')
import quivermutation as qm
from quivermutation import mutation, procedure, reduction, pathAlgebra, search, arrowPaths as ap, invariants
from helpers import path_algebra, quiet, arrow_set, rel_set

def perI(alg, v):
    """dim J_i = dim ker g_i for each i != v (copy of workshop/rounds/033/experimentalist_bothdie.py perI)."""
    Q = alg.quiver; rels = procedure.relationsFrom(alg); outs = ap.arrowsOutOf(Q, v); res = {}
    for i in Q.nodes:
        if i == v: continue
        P = ap.allPathsBetween(Q, i, v)
        if not P: continue
        dimV = len(P) - len(ap.idealBasis(Q, rels, i, v)); rows = []
        for q in P:
            row = {}
            for b in outs:
                for kk, x in ap.reduceAgainstPivots(ap.combination([q + (b,)]), ap.idealBasis(Q, rels, i, b[1])).items(): row[(b, kk)] = x
            rows.append(row)
        keys = sorted({k for r in rows for k in r}, key=repr)
        rk = sympy.Matrix([[r.get(k, 0) for k in keys] for r in rows]).rank() if keys else 0
        d = dimV - rk
        if d: res[i] = d
    return res

def show(a):
    return 'arrows %s | rels %s' % (sorted(arrow_set(a)), sorted(rel_set(a)))
def key(a): return search._coxeterKeyOrNone(a)
def step(name, alg, v, left=False):
    print('---', name, 'vertex', v, 'left' if left else 'right')
    base = key(alg); print(' start:', show(alg)); print(' key', base)
    print(' gate', mutation.mutationIsPossibleAtVertex(pathAlgebra.dualPathAlgebra(alg) if left else alg, v), '| J', perI(pathAlgebra.dualPathAlgebra(alg) if left else alg, v))
    f = mutation.leftQuiverMutationAtVertex if left else mutation.quiverMutationAtVertex
    try:
        raw = quiet(f, alg, v)
    except Exception as e:
        print(' RAISED', type(e).__name__, str(e)[:100]); return None
    print(' raw:', show(raw))
    red = quiet(reduction.reducePathAlgebra, raw); print(' reduced:', show(red)); k = key(red)
    print(' key after', k, '| guard (kept)' if k == base else '| guard REFUSES' if k is not None else '| key None (cyclic)')
    return red

def a2n3(n):
    N = 2 * n
    return path_algebra([(i, i + 1) for i in range(1, N)], [[(i, i + 1, i + 2, i + 3)] for i in range(1, N - 2)])
def rect(n):
    """kA_2 (x) kA_n: commutative squares, vertices (r, c) -> 2c + r + 1."""
    v = lambda r, c: 2 * c + r + 1
    arrows = [(v(0, c), v(1, c)) for c in range(n)] + [(v(r, c), v(r, c + 1)) for r in range(2) for c in range(n - 1)]
    rels = [[(v(0, c), v(1, c), v(1, c + 1)), (v(0, c), v(0, c + 1), v(1, c + 1))] for c in range(n - 1)]
    return path_algebra(arrows, rels)
def sect6(n, seq):
    cur = a2n3(n); base = key(cur); print('=== Fosse Sect. 6, A_{%d,3}, sequence %s; key of start %s' % (2 * n, seq, base))
    for v in seq:
        gate = mutation.mutationIsPossibleAtVertex(cur, v); j = perI(cur, v) if gate else None
        if not gate: print(' vertex', v, 'gate REFUSES'); return cur
        nxt = quiet(reduction.reducePathAlgebra, quiet(mutation.quiverMutationAtVertex, cur, v)); k = key(nxt)
        print(' mutate', v, '| gate True | J', j, '| key kept', k == base, '|', show(nxt)); cur = nxt
    return cur

def isoTo(a, b):
    """quiver-with-relations isomorphism (vertex bijection mapping arrows to arrows and the relation sets onto each other)."""
    Ga, Gb = nx.DiGraph(sorted(arrow_set(a))), nx.DiGraph(sorted(arrow_set(b)))
    if len(Ga) != len(Gb) or Ga.number_of_edges() != Gb.number_of_edges(): return None
    rb = rel_set(b)
    for m in nx.isomorphism.DiGraphMatcher(Ga, Gb).isomorphisms_iter():
        if {tuple(sorted(tuple(m[x] for x in p) for p in r)) for r in rel_set(a)} == rb: return m
    return None
def bfs(start, target, maxd, base):
    """BFS over gate-admitted, key-kept right mutations from start; first sequence reaching an algebra iso to target."""
    front = [(start, [])]; seen = {frozenset(rel_set(start)) | frozenset(arrow_set(start))}
    for d in range(maxd):
        nf = []
        for cur, seq in front:
            for v in sorted(cur.quiver.nodes):
                if not mutation.mutationIsPossibleAtVertex(cur, v): continue
                nxt = quiet(reduction.reducePathAlgebra, quiet(mutation.quiverMutationAtVertex, cur, v))
                if key(nxt) != base: continue
                sg = frozenset(rel_set(nxt)) | frozenset(arrow_set(nxt))
                if sg in seen: continue
                seen.add(sg)
                if isoTo(nxt, target) is not None: return seq + [v], nxt
                nf.append((nxt, seq + [v]))
        front = nf
    return None, None


def replay(n, seq, base):
    cur = a2n3(n)
    for v in seq:
        g = mutation.mutationIsPossibleAtVertex(cur, v); jj = perI(cur, v) if g else None
        nxt = quiet(reduction.reducePathAlgebra, quiet(mutation.quiverMutationAtVertex, cur, v))
        yield v, g, jj, key(nxt) == base
        cur = nxt
    yield 'final iso to rect(%d): %s' % (n, isoTo(cur, rect(n)) is not None), None, None, None

if __name__ == '__main__':
    import time; T0 = time.time()
    warnings_off = __import__('warnings').simplefilter('ignore')
    print('##### arXiv:2112.08129 Ex. 3.1 (acyclic)')

    A = path_algebra([(1,2),(2,3),(3,4),(3,5),(4,6),(5,6),(6,7)], [[(1,2,3)], [(3,4,6),(3,5,6)], [(4,6,7)], [(5,6,7)]])
    r = step('Fosse Ex 3.1', A, 3)
    want_arr = {(1,2),(2,4),(2,5),(4,3),(5,3),(3,6),(6,7)}
    want_rel = {((2,4,3),(2,5,3)), ((1,2,4),), ((1,2,5),), ((3,6,7),)}
    print(' printed final == repo:', arrow_set(r) == want_arr and rel_set(r) == want_rel)
    print('##### arXiv:2112.08129 Ex. 3.2 (cyclic, relation 1 ~> 1)')
    B = path_algebra([(1,2),(2,1)], [[(1,2,1)]])
    print('--- Fosse Ex 3.2 mutateAtVertex default:')
    try: quiet(qm.quiverMutationAtVertex, B, 1); print(' no error')
    except Exception as e: print(' RAISED', type(e).__name__)
    try:
        raw = quiet(procedure.toPathAlgebra, *procedure.mutateAtVertex(B.quiver, procedure.relationsFrom(B), 1, allowCyclic=True)); print(' allowCyclic=True raw:', show(raw))
    except Exception as e: print(' allowCyclic=True RAISED', type(e).__name__, e)
    print(' gate', mutation.mutationIsPossibleAtVertex(B, 1), 'key', key(B), 'J', perI(B, 1) if not nx.is_directed_acyclic_graph(B.quiver) else '')


    print('##### arXiv:2112.08129 Sect. 6, A_{2n,3} -> kA_2 (x) kA_n')
    for n, first in ((3, [1, 3]), (4, [1, 3, 5])):
        c = sect6(n, first); base = key(a2n3(n)); R = rect(n)
        print(' key(start) == key(rect):', base == key(R))
        if n == 3:
            e = {(1,4),(2,1),(2,3),(3,4),(4,5),(5,6)}
        s, res = bfs(c, R, 9, base); print(' BFS (gate + key guard) from the pushed-square quiver to the rectangle:', s)
        full = first + (s or [])
        print(' full replay with gate/J/key per step:'); [print('  ', t) for t in replay(n, full, base)]
    print('##### arXiv:1504.02617 Ex. 4.1 and 4.2 (A3/rad^2 at 1; 5-vertex at 2)')
    A = path_algebra([(1,2),(2,3)], [[(1,2,3)]])
    r = step('Opp Ex 4.1', A, 1); print(' printed final (phi*: 2->1, omega: 1->3, no relation) == repo:', arrow_set(r) == {(2,1),(1,3)} and rel_set(r) == set())
    step('Opp Ex 4.1 (repo left)', A, 1, left=True)
    B = path_algebra([(1,2),(2,3),(3,4),(2,5)], [[(1,2,3)],[(2,3,4)]])
    r = step('Opp Ex 4.2', B, 2)
    print(' printed final degree-0 part (beta*:3->2, delta*:5->2, s:2->4, delta alpha:1->5; relations d(alpha)=delta* delta alpha, d(s delta^-1)=s delta*) == repo:',
          arrow_set(r) == {(3,2),(5,2),(2,4),(1,5)} and rel_set(r) == {((1,5,2),),((5,2,4),)})
    print(' printed "quasi-isomorphic to A^op": repo result iso to A^op:', isoTo(r, pathAlgebra.dualPathAlgebra(B)) is not None, '| iso to A:', isoTo(r, B) is not None)
    step('Opp Ex 4.2 (repo left)', B, 2, left=True)
    print('time %.0fs' % (time.time() - T0))
