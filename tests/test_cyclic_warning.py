"""A cyclic quiver is outside the setting where the Cartan matrix, ideal membership and
minimal generators are proved correct (paths are simple in the vertices, research E-193):
they warn, and stay silent on an acyclic quiver.  Round 067, thread T10/T5."""
import warnings

import networkx as nx
import pytest

from quivermutation import arrowPaths as ap
from quivermutation import procedure


def quiver(arrows):
    q = nx.MultiDiGraph()
    for a, b in arrows:
        q.add_edge(a, b)
    return q


def test_hasOrientedCycle():
    assert not ap.hasOrientedCycle(quiver([(1, 2), (2, 3), (1, 3), (1, 3)]))
    assert ap.hasOrientedCycle(quiver([(1, 2), (2, 1)]))
    assert ap.hasOrientedCycle(quiver([(1, 1), (1, 2)]))


def test_acyclic_is_silent():
    q = quiver([(1, 2), (2, 3)])
    rel = [{((1, 2, 0), (2, 3, 0)): 1}]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ap.cartanMatrix(q, rel)
        ap.isInIdeal(q, rel, rel[0])
        procedure._minimalGenerators(q, rel)


def test_cyclic_warns_everywhere_and_keeps_the_answer():
    q = quiver([(1, 2), (2, 3), (3, 1)])
    rel = [{((1, 2, 0), (2, 3, 0)): 1}]
    with pytest.warns(ap.CyclicQuiverWarning, match="oriented cycle"):
        m = ap.cartanMatrix(q, rel)
    with pytest.warns(ap.CyclicQuiverWarning):
        assert ap.isInIdeal(q, rel, rel[0])
    with pytest.warns(ap.CyclicQuiverWarning):
        assert procedure._minimalGenerators(q, rel) == rel
    with warnings.catch_warnings():  # the warning does not change the value
        warnings.simplefilter("ignore")
        assert ap.cartanMatrix(q, rel) == m


def test_mutateAtVertex_refuses_cyclic_descendants():
    """E-211: a:1->2, s:2->3, b:3->4, c:4->3 with as, sbc, cb, mutated at 1: D(1) contains the 3 <-> 4 cycle."""
    q = nx.MultiDiGraph()
    names = {}
    for n, (x, y) in {"a": (1, 2), "s": (2, 3), "b": (3, 4), "c": (4, 3)}.items():
        names[n] = (x, y, q.add_edge(x, y))
    A = lambda *ns: tuple(names[n] for n in ns)
    rel = [{A("a", "s"): 1}, {A("s", "b", "c"): 1}, {A("c", "b"): 1}]
    with pytest.raises(procedure.CyclicDescendantsError, match="E-211"):
        procedure.mutateAtVertex(q, rel, 1)
    procedure.mutateAtVertex(q, rel, 1, allowCyclic=True)        # old behaviour kept on request
    # a cycle NOT below the vertex is not refused: mutate at 3? D(3) = {3,4} is cyclic; at 2 -> D(2) = {2,3,4} cyclic;
    # a cycle above the vertex: 1 <-> 2 and 2 -> 3, mutate at 3 (D(3) = {3})
    q2 = quiver([(1, 2), (2, 1), (2, 3)])
    procedure.mutateAtVertex(q2, [], 3)


def _e211():
    q = nx.MultiDiGraph()
    names = {}
    for n, (x, y) in {"a": (1, 2), "s": (2, 3), "b": (3, 4), "c": (4, 3)}.items():
        names[n] = (x, y, q.add_edge(x, y))
    A = lambda *ns: tuple(names[n] for n in ns)
    return q, [{A("a", "s"): 1}, {A("s", "b", "c"): 1}, {A("c", "b"): 1}]


def test_left_and_multi_mutation_refuse_cyclic_descendants():
    """Left mutation applies the guard to D in the opposite quiver, i.e. the ancestors in the original (r078 review)."""
    q, rel = _e211()
    with pytest.raises(procedure.CyclicDescendantsError):
        procedure.mutateAtVertices(q, rel, [1])
    qo, relo = procedure._opposite(q, rel)
    with pytest.raises(procedure.CyclicDescendantsError):
        procedure.mutateLeftAtVertex(q, rel, 4)          # ancestors of 4 in q = {1, 2, 3, 4}: cyclic 3 <-> 4
    with pytest.raises(procedure.CyclicDescendantsError):
        procedure.mutateAtVertices(q, rel, [-4])
    with pytest.raises(procedure.CyclicDescendantsError):
        procedure.mutateAtVertices(qo, relo, [4])
    # nothing above 1 in q: left mutation at 1 is not refused
    procedure.mutateLeftAtVertex(q, rel, 1)


def test_loop_at_or_below_the_vertex_refused():
    q = nx.MultiDiGraph()
    q.add_edge(1, 2)
    q.add_edge(2, 2)
    with pytest.raises(procedure.CyclicDescendantsError):       # loop below 1
        procedure.mutateAtVertex(q, [], 1)
    q2 = nx.MultiDiGraph()
    q2.add_edge(1, 1)
    with pytest.raises(procedure.CyclicDescendantsError):       # loop at the vertex itself
        procedure.mutateAtVertex(q2, [], 1)
    q3 = nx.MultiDiGraph()
    q3.add_edge(1, 2)
    q3.add_edge(1, 1)
    procedure.mutateAtVertex(q3, [], 2)                         # loop above 2 only: accepted
