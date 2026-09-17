"""Reorienting a relation-free tree is a sequence of mutations.

The repo has always used the fact that the underlying undirected tree of a
relation-free quiver settles its class.  That is a statement about *derived*
equivalence, justified by BGP reflections; the classification then merges
*mutation* classes on it.  These tests close the step in between: the reflection
is what the procedure does at a source, the orientations of a tree are one
mutation class, and a sequence joining any two of them can be written down.
"""

import itertools

import networkx as nx
import pytest

import quivermutation as qm
from quivermutation import invariants as inv
from quivermutation import nakayama as nk
from quivermutation import pathAlgebra as pa
from quivermutation import reflections as rf
from quivermutation import treeSearch as ts


def treeQuivers(order):
    """Every tree of the order, in every orientation."""
    for tree in ts.treesOfOrder(order):
        graph = nx.convert_node_labels_to_integers(tree, first_label = 1, ordering = "sorted")
        for arrows in rf.orientationsOfTree(graph):
            yield graph, rf.quiverWithArrows(graph.nodes, arrows)


def reversedAt(pathAlg, vertex):
    """The same quiver with every arrow at one vertex turned round."""
    return rf.quiverWithArrows(
        pathAlg.quiver.nodes,
        [(head, tail) if vertex in (tail, head) else (tail, head)
         for tail, head in pathAlg.quiver.edges()])


@pytest.mark.parametrize("order", [4, 5, 6])
def test_mutating_at_a_source_is_the_reflection(order):
    """Right mutation at a source reverses its arrows and leaves no relations.

    This is the step the whole thing rests on, and it is not obvious from the
    procedure: the mutation of a general quiver creates relations and renumbers.
    At a source of a relation-free tree it does neither.
    """
    checked = 0
    for _graph, quiver in treeQuivers(order):
        for vertex in quiver.quiver.nodes:
            if quiver.quiver.in_degree(vertex):
                continue
            assert quiver.canMutateAt(vertex)
            mutated = qm.quiverMutationAtVertex(quiver, vertex).reduce()
            assert not mutated.rels
            assert rf.arrowsOf(mutated) == rf.arrowsOf(reversedAt(quiver, vertex))
            checked += 1
    assert checked


@pytest.mark.parametrize("order", [4, 5, 6])
def test_left_mutating_at_a_sink_is_the_reflection_too(order):
    """The mirror, which is what runs a reorientation the other way round."""
    checked = 0
    for _graph, quiver in treeQuivers(order):
        for vertex in quiver.quiver.nodes:
            if quiver.quiver.out_degree(vertex):
                continue
            mutated = qm.leftQuiverMutationAtVertex(quiver, vertex).reduce()
            assert not mutated.rels
            assert rf.arrowsOf(mutated) == rf.arrowsOf(reversedAt(quiver, vertex))
            checked += 1
    assert checked


@pytest.mark.parametrize("order", [4, 5, 6])
def test_every_orientation_is_reached_and_the_sequence_is_legal(order):
    """The sequence lands on the orientation asked for, by legal mutations only.

    `reorient` checks each step against the procedure's own admissibility --
    right mutations directly, left ones on the opposite algebra -- and compares
    the result with the target, so this is the engine confirming the
    construction rather than the construction asserting itself.
    """
    for graph, quiver in itertools.islice(treeQuivers(order), 0, None, 3):
        for target in rf.orientationsOfTree(graph):
            reached, sequence = rf.reorient(quiver, target)
            assert rf.arrowsOf(reached) == target
            assert not reached.rels
            assert inv.coxeterKey(reached) == inv.coxeterKey(quiver)
            assert len(sequence) <= order * graph.number_of_edges()


def test_the_reflection_sequence_refuses_what_it_cannot_do():
    line = rf.quiverWithArrows([1, 2, 3], [(1, 2), (2, 3)])
    with pytest.raises(ValueError):
        rf.reflectionSequence(rf.arrowsOf(line), frozenset([(1, 2), (3, 1)]))
    withRelations = nk.LinearNakayamaAlgebra.fromClassName("300")
    assert not rf.isReflectable(withRelations)
    with pytest.raises(ValueError):
        rf.reorient(withRelations, rf.arrowsOf(withRelations))


@pytest.mark.parametrize("order", [4, 5, 6, 7])
def test_right_mutations_alone_already_connect_every_orientation(order):
    """Reflecting at sources is enough; sinks only make the path shorter.

    Worth pinning because the search walks right mutations only.  The distances
    are the reason the lemma is not just a tidier way to say what a search would
    find anyway: at order 7 two orientations can be 12 right mutations apart,
    twice the depth a classification run can afford.
    """
    for tree in ts.treesOfOrder(order):
        graph = nx.convert_node_labels_to_integers(tree, first_label = 1, ordering = "sorted")
        start = next(iter(rf.orientationsOfTree(graph)))
        seen = {start}
        frontier = [start]
        while frontier:
            arrows = frontier.pop()
            heads = {head for _tail, head in arrows}
            for vertex in graph.nodes:
                if vertex in heads:
                    continue
                flipped = frozenset((head, tail) if vertex in (tail, head) else (tail, head)
                                    for tail, head in arrows)
                if flipped not in seen:
                    seen.add(flipped)
                    frontier.append(flipped)
        assert len(seen) == 2 ** graph.number_of_edges()


def test_the_line_has_an_orientation_that_is_an_lna():
    """A path-shaped tree reorients to `kA_n`, which is an LNA of the length."""
    graph = nx.path_graph(range(1, 7))
    zigzag = rf.quiverWithArrows(graph.nodes, [(1, 2), (3, 2), (3, 4), (5, 4), (5, 6)])
    reached, sequence = rf.reorient(zigzag, rf.linearOrientation(graph))
    assert sequence
    assert rf.arrowsOf(reached) == rf.arrowsOf(nk.LinearNakayamaAlgebra(6, "0000"))
    assert rf.linearOrientation(nx.star_graph(3)) is None


def sameAlgebra(first, second):
    """Whether two path algebras are isomorphic, quiver and relations together."""
    left = nx.DiGraph(first.quiver.edges())
    right = nx.DiGraph(second.quiver.edges())
    for mapping in nx.algorithms.isomorphism.DiGraphMatcher(left, right).isomorphisms_iter():
        carried = sorted(tuple(sorted(tuple(mapping[vertex] for vertex in path)
                                      for path in rel)) for rel in first.rels)
        theirs = sorted(tuple(sorted(tuple(path) for path in rel)) for rel in second.rels)
        if carried == theirs:
            return True
    return False


@pytest.mark.parametrize("first, second", [
    ("00300", "05000"),
    ("00030", "33000"),
    ("00000", "00202"),
])
def test_the_bridge_lands_on_the_other_algebra(first, second):
    """Two LNAs reaching one tree are joined by an exhibited mutation sequence.

    This is what makes a merge on a shared hereditary form constructive.  The
    reflection leg of `00300` to `05000` is nine mutations on its own, and the
    whole sequence sixteen: no classification search is run anywhere near that
    deep, so the sequence is one that only the construction can supply.
    """
    left = nk.LinearNakayamaAlgebra.fromClassName(first)
    right = nk.LinearNakayamaAlgebra.fromClassName(second)
    bridge = rf.mutationBridge(left, right, 4)
    assert bridge is not None
    assert sameAlgebra(bridge['algebra'], right)
    assert bridge['sequence'] == bridge['toTree'] + bridge['reflect'] + bridge['back']
    assert all(1 <= abs(step) <= left.length for step in bridge['sequence'])


def test_the_bridge_says_nothing_when_no_tree_is_reached():
    """`3033030` is in no quipu class, so nothing hereditary is ever reached."""
    left = nk.LinearNakayamaAlgebra.fromClassName("3033030")
    right = nk.LinearNakayamaAlgebra.fromClassName("3345000")
    assert rf.mutationBridge(left, right, 3) is None
