"""The two guards on rule discovery, both of which exist because of a false rule.

`discover.py` searches for rewrites and then has to decide which to believe. Two
things make that decision safe, and neither is obvious:

* a candidate whose tight window fails is retried with the window **widened**,
  since `describeLink` returns the smallest window containing what the link
  touches and that is often too small to state the rule's precondition;
* widening narrows where a rule fires, so a padded form can fire once or twice
  and pass vacuously. A candidate counts as verified only with enough
  confirmations behind it. That is R-008, which cost a false rule.
"""

import pytest

import discover
import lnaMoves


def test_maximum_overlap_is_what_the_theorem_needs():
    """Zero or one is 'almost separate'; two or more is what blocks the table."""
    assert discover.maximumOverlap([0, 0, 0, 0]) == 0
    assert discover.maximumOverlap([3, 0, 0, 3]) == 0        # far apart
    assert discover.maximumOverlap([3, 0, 3, 0]) == 1        # meeting at a vertex
    assert discover.maximumOverlap([3, 3, 0, 0]) == 2        # maximally overlapping
    assert discover.maximumOverlap([4, 4, 0, 0]) == 3


def test_padding_a_rule_moves_its_window_without_changing_what_it_does():
    """A padded rule fired one arrow further left is the original rule.

    This is the index arithmetic that makes widening safe: a relation at relative
    start s sits at absolute arrow windowStart + s - 1, and a signed offset v
    names vertex windowStart + |v| - 1, so shifting the window left by one adds
    one to every start and to every offset's magnitude.
    """
    tight = (5, ((1, 3), (2, 3)), ((0, 3), (1, 3)), (2, 2))
    padded = discover.paddedRule(tight, 1, 1)
    assert padded[0] == tight[0] + 2

    length = 11
    relLengths = [0, 0, 3, 3, 0, 0, 0, 0, 0]
    for windowStart in range(2, length - tight[0]):
        if not lnaMoves.matchesAt(length, relLengths, tight, windowStart):
            continue
        tightResult = lnaMoves.applyAt(length, relLengths, tight, windowStart)
        paddedResult = lnaMoves.applyAt(length, relLengths, padded, windowStart - 1)
        assert paddedResult is not None, "the padded rule should still apply"
        assert tightResult == paddedResult


def test_padding_only_ever_narrows_where_a_rule_fires():
    """The wider window is a stronger precondition, never a weaker one."""
    tight = (5, ((1, 3), (2, 3)), ((0, 3), (1, 3)), (2, 2))
    padded = discover.paddedRule(tight, 1, 1)
    length = 12
    relLengths = [2, 0, 0, 3, 3, 0, 0, 0, 0, 0]

    tightPositions = {w for w in range(1, length)
                      if lnaMoves.matchesAt(length, relLengths, tight, w)}
    paddedPositions = {w + 1 for w in range(1, length)
                       if lnaMoves.matchesAt(length, relLengths, padded, w)}
    assert paddedPositions <= tightPositions


def test_the_false_rule_of_R_008_is_rejected():
    """The rewrite that passed on one confirmation, and should not have.

    It is admitted by zero-failures alone at lengths 8 and 9, and fails as soon
    as there is room for it to fire properly. Checking at length 10 alone is
    enough to see it.
    """
    false = (8, ((2, 3), (3, 3)), ((2, 2), (5, 3)), (4, -7, -8, 4))
    confirmed, failures = lnaMoves.verifyMove(false, [10], checkCoxeter=True)
    assert failures, "R-008's rewrite must not verify clean"
    assert confirmed < 8, "and it must not reach the confirmation threshold either"
