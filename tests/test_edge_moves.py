"""The doubling at an end, which no window rule can state.

A relation at the source gains a copy at the next vertex under two left
mutations, and dually at the sink.  `lnaMoves`' encoding cannot express it --
its window is the arrows the new relation covers, and a spectator relation is
allowed to start on the window's last arrow and run out the far side, which
`matchesAt` refuses.  So these are checked directly against the mutation engine,
the same three ways `verifyMove` checks a rule: the predicted LNA, every
mutation admissible, and the Coxeter polynomial held fixed (research F-029).
"""

import pytest

from quivermutation import edgeMoves as em
from quivermutation import lines
from quivermutation import lnaMoves as lm
from quivermutation import mutation
from quivermutation import nakayama as nk


def _engineAgrees(length, relLengths, predicted, sequence):
    """Run the sequence and check all three conditions."""
    startAlg = nk.LinearNakayamaAlgebra(length, list(relLengths))
    if not lm.isLegalSequence(startAlg, list(sequence)):
        return 'illegal mutation'
    mutated = lm._quiet(mutation.quiverMutationAtVertices,
                        lm._copy(startAlg), list(sequence))
    actual = lm.asRelLengths(mutated, length)
    if actual is None:
        return 'not a line'
    if tuple(actual) != tuple(predicted):
        return 'got ' + lines.className(actual)
    if not lm._sameCoxeter(length, list(relLengths), list(predicted)):
        return 'Coxeter polynomial moved'
    return None


# ---------------------------------------------------------------------------
# What the moves do
# ---------------------------------------------------------------------------


def test_the_source_doubling_copies_the_relation_to_the_next_vertex():
    assert em.sourceDoubling(12, (3, 0, 0, 0, 0, 0, 0, 0, 0, 0)) == (
        (3, 3, 0, 0, 0, 0, 0, 0, 0, 0), [-4, -4])


def test_the_source_doubling_tolerates_a_relation_the_window_overlaps():
    """The point of the module: the companion starts on the window's last arrow."""
    assert em.sourceDoubling(12, (3, 0, 0, 3, 0, 0, 0, 0, 0, 0)) == (
        (3, 3, 0, 3, 0, 0, 0, 0, 0, 0), [-4, -4])
    assert lm.rewritesOf(12, [3, 0, 0, 3, 0, 0, 0, 0, 0, 0], lm.ALL_MOVES) == []


def test_the_source_doubling_refuses_a_relation_further_inside_the_window():
    """A relation starting at vertex 3 is inside the arrows the new one covers."""
    assert em.sourceDoubling(8, (3, 0, 3, 0, 0, 0)) is None


def test_the_sink_doubling_is_the_dual_of_the_source_one():
    assert em.sinkDoubling(12, (0, 0, 0, 0, 0, 0, 0, 0, 3, 0)) == (
        (0, 0, 0, 0, 0, 0, 0, 3, 3, 0), [9, 9])


def test_the_collapses_undo_the_doublings():
    for length, relLengths in [(12, (3, 0, 0, 3, 0, 0, 0, 0, 0, 0)),
                               (12, (5, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
                               (10, (4, 0, 0, 0, 0, 0, 0, 0))]:
        doubled, _ = em.sourceDoubling(length, relLengths)
        assert em.sourceCollapse(length, doubled)[0] == tuple(relLengths)
    for length, relLengths in [(12, (0, 0, 0, 0, 0, 0, 0, 0, 3, 0)),
                               (12, (0, 0, 0, 0, 0, 0, 4, 0, 0, 0))]:
        answer = em.sinkDoubling(length, relLengths)
        if answer is None:
            continue
        doubled, _ = answer
        assert em.sinkCollapse(length, doubled)[0] == tuple(relLengths)


# ---------------------------------------------------------------------------
# Against the engine, over every LNA that carries the pattern
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("length", [7, 8, 9, 10])
def test_every_edge_move_holds_against_the_mutation_engine(length):
    """No failures anywhere, which is what makes these usable without a search."""
    confirmed = 0
    failures = []
    for relLengths in nk.allRelationLengths(length):
        for predicted, sequence in em.rewritesOf(length, relLengths):
            reason = _engineAgrees(length, relLengths, predicted, sequence)
            if reason is None:
                confirmed += 1
            else:
                failures.append((lines.className(relLengths),
                                 lines.className(predicted), reason))
    assert failures == []
    assert confirmed > 0


@pytest.mark.parametrize("length, count", [(7, 88), (8, 256), (9, 784)])
def test_the_number_of_edge_moves_available_is_pinned(length, count):
    total = sum(len(em.rewritesOf(length, relLengths))
                for relLengths in nk.allRelationLengths(length))
    assert total == count


@pytest.mark.parametrize("length", [7, 8, 9])
def test_each_move_fires_exactly_as_often_as_its_dual(length):
    """The relation dual of F-026 exchanges source for sink, so the counts match."""
    fired = {move.__name__: sum(1 for relLengths in nk.allRelationLengths(length)
                                if move(length, list(relLengths)) is not None)
             for move in em.MOVES}
    assert fired['sourceDoubling'] == fired['sinkDoubling']
    assert fired['sourceCollapse'] == fired['sinkCollapse']
    assert fired['sourceDoubling'] == fired['sourceCollapse']
