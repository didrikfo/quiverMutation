"""`proposition:doubleMutation` of arXiv:2310.08346, against the engine.

The proposition is proved in the paper, so the point here is not to doubt it but
to check that `doubleMutation` reads it right -- the conditions, the direction of
each length change, the dual taken through F-026, the `s = 1` extension the paper
does not state -- and that the sequence `[-t, -t]` is what the engine performs.
Checked the same three ways `lnaMoves.verifyMove` checks a rule: the predicted
LNA, every mutation admissible, and the Coxeter polynomial held fixed (F-032).
"""

import pytest

from quivermutation import doubleMutation as dm
from quivermutation import edgeMoves as em
from quivermutation import freeMoves as fm
from quivermutation import lines
from quivermutation import lnaMoves as lm
from quivermutation import mutation
from quivermutation import nakayama as nk


def _engineAgrees(length, relLengths, predicted, sequence):
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


def _row(length, intervals):
    return dm.relLengthsOf(length, intervals)


# ---------------------------------------------------------------------------
# The paper's own worked examples
# ---------------------------------------------------------------------------


def test_example_mutationToA11_5_first_two_steps():
    """Lambda_1 = L_8(Lambda) and Lambda_2 = L_9(Lambda_1), as drawn in the paper."""
    start = _row(11, [(1, 5), (2, 8), (5, 11)])
    first, sequence = dm.leftDoubleMutation(11, start, 2)
    assert dm.intervalsOf(first) == [(1, 6), (2, 8), (3, 9), (6, 11)]
    assert sequence == [-8, -8]
    second, _ = dm.leftDoubleMutation(11, first, 3)
    assert dm.intervalsOf(second) == [(1, 7), (3, 9), (4, 10), (7, 11)]


def test_example_A13tworelations_right_mutations_at_the_source():
    """R_1 twice, the dual the paper names but does not state."""
    start = _row(13, [(1, 8), (3, 9), (5, 11), (6, 13)])
    once, sequence = dm.rightDoubleMutation(13, start, 8)
    assert dm.intervalsOf(once) == [(1, 8), (2, 9), (4, 11), (5, 13)]
    assert sequence == [1, 1]
    twice, _ = dm.rightDoubleMutation(13, once, 8)
    assert dm.intervalsOf(twice) == [(1, 8), (3, 11), (4, 13)]


# ---------------------------------------------------------------------------
# What it contains
# ---------------------------------------------------------------------------


def test_on_an_isolated_equal_pair_it_is_the_pair_slide():
    """The companion lengthens onto r, stops being minimal, and goes."""
    before = _row(12, [(4, 7), (5, 8)])
    after, sequence = dm.leftDoubleMutation(12, before, 5)
    assert dm.intervalsOf(after) == [(5, 8), (6, 9)]
    assert sequence == [-8, -8]


def test_at_the_sink_it_is_the_end_pair_collapse():
    before = _row(10, [(6, 9), (7, 10)])
    after, sequence = dm.leftDoubleMutation(10, before, 7)
    assert dm.intervalsOf(after) == [(7, 10)]
    assert sequence == [-10, -10]


@pytest.mark.parametrize("length", [7, 8, 9])
def test_every_source_doubling_is_a_double_mutation_at_s_equal_one(length):
    """F-029 is the `s = 1` extension with nothing crossing r."""
    for relLengths in nk.allRelationLengths(length):
        doubled = em.sourceDoubling(length, relLengths)
        if doubled is not None:
            assert doubled == dm.leftDoubleMutation(length, relLengths, 1)


def test_a_relation_starting_at_t_minus_one_blocks_it():
    assert dm.leftDoubleMutation(10, _row(10, [(1, 4), (2, 6), (5, 8)]), 2) is None


def test_the_paper_hypothesis_is_needed_away_from_the_source():
    """No relation at s - 1 and s > 1: not the proposition, and not attempted."""
    assert dm.leftDoubleMutation(10, _row(10, [(3, 6)]), 3) is None


# ---------------------------------------------------------------------------
# Against the engine, everywhere it applies
# ---------------------------------------------------------------------------


def _verify(length):
    confirmed, failures = 0, []
    for relLengths in nk.allRelationLengths(length):
        for predicted, sequence in dm.rewritesOf(length, relLengths):
            reason = _engineAgrees(length, relLengths, predicted, sequence)
            if reason is None:
                confirmed += 1
            else:
                failures.append((lines.className(relLengths),
                                 lines.className(predicted), sequence, reason))
    return confirmed, failures


@pytest.mark.parametrize("length, count", [(5, 18), (6, 68), (7, 250), (8, 922)])
def test_every_double_mutation_holds_against_the_engine(length, count):
    confirmed, failures = _verify(length)
    assert failures == []
    assert confirmed == count


@pytest.mark.slow
@pytest.mark.parametrize("length, count", [(9, 3430), (10, 12868)])
def test_every_double_mutation_holds_against_the_engine_further_out(length, count):
    confirmed, failures = _verify(length)
    assert failures == []
    assert confirmed == count


@pytest.mark.parametrize("length", [7, 8, 9, 10])
def test_left_and_right_fire_equally_often(length):
    """The relation dual exchanges them, so the counts must match."""
    left = right = 0
    for relLengths in nk.allRelationLengths(length):
        for a, b in dm.intervalsOf(relLengths):
            left += dm.leftDoubleMutation(length, relLengths, a) is not None
            right += dm.rightDoubleMutation(length, relLengths, b) is not None
    assert left == right > 0


# ---------------------------------------------------------------------------
# What it reaches
# ---------------------------------------------------------------------------


def test_with_the_free_move_n_9_leaves_exactly_the_two_classes_no_quipu_names():
    """3345000, C(2,4,4), eight members; and 3033030, not piecewise hereditary."""
    result = fm.coverage(9, rules = [], free = True, doubles = True)
    left = sorted(lines.className(lna) for lna in result['uncovered'])
    assert len(left) == 9
    assert '3033030' in left and '3345000' in left
    orbitsLeft = [members for members in result['orbits'].values()
                  if not set(members) & result['covered']]
    assert sorted(len(members) for members in orbitsLeft) == [1, 8]


def test_n_8_needs_neither_the_table_nor_the_free_move():
    result = fm.coverage(8, rules = [], free = False, doubles = True)
    assert result['uncovered'] == []
