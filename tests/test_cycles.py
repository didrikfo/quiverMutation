"""A cycle in the quiver used to crash the search with a RecursionError.

The relation-enumerating functions recurse along the arrows out of a vertex.
Without tracking which vertices are already on the recursion path, a cycle makes
them descend forever.  mutationSearchDepthFirst computed all the relations of a
quiver at the top of every node, *before* testing that quiver for cycles, so the
first mutation that produced a cyclic quiver killed the whole search on the
following node -- which is what stopped the length-12 run.

Restricting to simple paths is also the right semantics: the rest of the module
works with simple paths throughout, and on an acyclic quiver, where a vertex
cannot repeat on a path anyway, nothing changes.
"""

import sys

import networkx as nx
import pytest

from quivermutation import pathAlgebra as pac
import quivermutation as qm
from helpers import line_algebra, quiet


def cyclic_algebra():
    """1 -> 2 -> 3 -> 1, with the path 1 -> 2 -> 3 set to zero."""
    pa = pac.PathAlgebra()
    pa.add_arrows_from([[1, 2], [2, 3], [3, 1]])
    pa.add_rels_from([[[1, 2, 3]]])
    return pa


@pytest.fixture
def shallow_recursion():
    """Lower the recursion limit so an unbounded recursion fails fast."""
    original = sys.getrecursionlimit()
    sys.setrecursionlimit(300)
    yield
    sys.setrecursionlimit(original)


def test_the_quiver_really_has_a_cycle():
    assert list(nx.simple_cycles(cyclic_algebra().quiver))


def test_all_rels_between_vertices_terminates_on_a_cycle(shallow_recursion):
    pa = cyclic_algebra()
    assert qm.allRelsBetweenVertices(pa, 1, 3) == [[[1, 2, 3]]]
    assert qm.allRelsBetweenVertices(pa, 2, 1) == []


def test_all_rels_in_path_algebra_terminates_on_a_cycle(shallow_recursion):
    rels = qm.allRelsInPathAlgebra(cyclic_algebra())
    assert [[1, 2, 3]] in rels


def test_extend_rel_terminates_on_a_cycle(shallow_recursion):
    pa = cyclic_algebra()
    extended = qm.extendRel(pa, [[1, 2, 3]])
    assert [[1, 2, 3]] in extended
    # It may extend along 3 -> 1, but must not go round again.
    for rel in extended:
        for path in rel:
            assert len(path) == len(set(path)), f"{path} revisits a vertex"


def test_the_search_ends_the_branch_at_a_cycle_instead_of_crashing(shallow_recursion):
    """The search never descends from a cyclic quiver, so it should just stop."""
    collected = []
    quiet(qm.mutationSearchDepthFirst, cyclic_algebra(), 3, [], "cyclic",
          printOutput=False, collected=collected)
    assert collected == []      # a 3-cycle is not a line, so nothing is recorded


@pytest.mark.parametrize("length, rels", [(5, "300"), (6, "3030"), (7, "22230")])
def test_bounding_the_recursion_changes_nothing_on_an_acyclic_quiver(length, rels):
    """No vertex can repeat on a path in an acyclic quiver, so the visited set
    never excludes anything there."""
    pa = line_algebra(length, rels)
    for source in pa.vertices():
        for target in pa.vertices():
            found = qm.allRelsBetweenVertices(pa, source, target)
            for rel in found:
                for path in rel:
                    assert path[0] == source and path[-1] == target
                    assert len(path) == len(set(path))
