"""Relation overlap: the coordinate the quipu theorem's reach is measured in.

Two consecutive relations of an LNA share some number of arrows, and the
theorem of arXiv:2305.06642 covers exactly the algebras where that number never
exceeds one.  These tests pin what overlap the code computes, that the rules
anchored to an end of the quiver hold there and nowhere else, and the two
measurements F-021 rests on: nothing below overlap two is ever left over, and
nothing at overlap two or more is reached without an end of the quiver or a
third overlapping relation.
"""

import pytest

from quivermutation import lnaMoves as lm
from quivermutation import nakayama as nk
from quivermutation import overlap as ov


# ---------------------------------------------------------------------------
# The measurement itself
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("relLengths, profile", [
    ([0, 0, 0, 0], ()),
    ([2, 0, 0, 0], ()),
    ([3, 3, 0, 0], (2,)),               # arrows 1-3 and 2-4, sharing 2 and 3
    ([3, 0, 2, 0], (1,)),               # arrows 1-3 and 3-4, sharing 3
    ([3, 0, 0, 2], (0,)),               # arrows 1-3 and 4-5, sharing nothing
    ([4, 4, 0, 0], (3,)),
    ([3, 3, 3, 0], (2, 2)),
])
def test_the_overlap_profile_counts_the_arrows_two_relations_share(relLengths, profile):
    assert ov.overlapProfile(relLengths) == profile
    assert ov.maxOverlap(relLengths) == (max(profile) if profile else 0)


def test_almost_separate_is_exactly_overlap_at_most_one():
    """The paper's condition and the overlap bound are the same condition.

    This is what makes overlap the right coordinate: the set the theorem covers
    is a sublevel set of `maxOverlap`, so "what the theorem misses" and "what
    overlaps by two or more" are the same question.
    """
    for length in range(4, 9):
        for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
            relLengths = algebra.relLengths
            assert (ov.isAlmostSeparate(length, relLengths)
                    == (ov.maxOverlap(relLengths) <= 1)), relLengths


def test_overlap_runs_break_where_the_relations_stop_overlapping_heavily():
    runs = ov.overlapRuns([3, 3, 0, 0, 3, 3, 0])
    assert runs == [[(1, 3), (2, 3)], [(5, 3), (6, 3)]]
    # A relation sharing only one arrow starts a new run, which is the point:
    # it is not part of the heavily overlapping cluster.
    assert ov.overlapRuns([2, 3, 3, 0, 0]) == [[(1, 2)], [(2, 3), (3, 3)]]


# ---------------------------------------------------------------------------
# Rules anchored to an end of the quiver
# ---------------------------------------------------------------------------


def test_an_anchored_rule_is_only_ever_tried_at_its_end():
    """`windowStartsFor` is the single gate every caller slides a rule through,
    so an anchored rule cannot be applied in the interior by any route."""
    left = (4, ((0, 3), (1, 3)), ((0, 3),), (1, 1), 'left')
    right = (4, ((0, 3), (1, 3)), ((1, 3),), (-5, -5), 'right')
    floating = (4, ((0, 2),), ((2, 2),), (-3, -4))
    assert lm.anchorOf(left) == 'left'
    assert lm.anchorOf(floating) is None
    assert lm.windowStartsFor(9, left) == [1]
    assert lm.windowStartsFor(9, right) == [5]   # arrows 5-8, the last one
    assert lm.windowStartsFor(9, floating) == [1, 2, 3, 4, 5]
    # The window has to fit at all.
    assert lm.windowStartsFor(4, left) == []


@pytest.mark.parametrize("rule", lm.endPairCollapseRules(6), ids=str)
def test_the_end_pair_collapse_holds_at_its_end(rule):
    """A maximally overlapping pair flush against an end of the quiver loses one
    of its two relations under two mutations at the end vertex.

    This is the rule that gets a heavily overlapping LNA into the theorem's
    reach at all, and the only kind of rule that reduces an isolated overlap.
    """
    width = rule[0]
    confirmed, failures = lm.verifyMove(rule, range(width + 1, width + 5))
    assert failures == []
    assert confirmed > 0


def test_the_end_pair_collapse_is_false_in_the_interior():
    """The reason it has to be anchored, stated as a check.

    Taking the same rewrite as a floating rule -- the identical window, relations
    and mutation offsets, without the anchor -- makes it false, because away from
    the end the mutations are not the same mutations.  Stating a rule that needs
    an end as though it held everywhere is how R-009's false rules arose.
    """
    anchored = (4, ((0, 3), (1, 3)), ((0, 3),), (1, 1), 'left')
    floating = anchored[:4]
    confirmed, failures = lm.verifyMove(floating, range(5, 9))
    assert failures, "the unanchored form must fail somewhere"
    assert floating not in lm.VERIFIED_MOVES


def test_the_anchored_rules_are_in_the_table_and_the_floating_half_is_unchanged():
    assert set(lm.endPairCollapseRules()) <= set(lm.ANCHORED_MOVES)
    assert lm.ALL_MOVES == lm.VERIFIED_MOVES + lm.ANCHORED_MOVES
    assert all(lm.anchorOf(rule) is None for rule in lm.VERIFIED_MOVES)
    assert all(lm.anchorOf(rule) is not None for rule in lm.ANCHORED_MOVES)


# ---------------------------------------------------------------------------
# What the whole table reaches, by overlap
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("length", [5, 6, 7, 8])
def test_the_theorem_and_the_moves_place_everything_below_overlap_two(length):
    """Half of F-021: the gap is entirely above the almost separate line.

    Every LNA the seeding plus the move orbits fail to place has two consecutive
    relations sharing two or more arrows.  So "the rows a search still has to
    find" and "the heavily overlapping rows" are the same set, and there is no
    second phenomenon hiding among the rest.
    """
    result = ov.coverage(length)
    for relLengths in result['uncovered']:
        assert ov.maxOverlap(relLengths) >= 2, relLengths
    for overlapValue, (_total, _covered, left) in result['byMaxOverlap'].items():
        if overlapValue <= 1:
            assert left == 0


def test_a_rewrite_reaching_an_overlapping_pair_needs_an_end_or_a_third_relation():
    """The other half of F-021, as a property of the rule table.

    No rule in the table takes a window holding exactly two relations that
    overlap in two or more arrows to anything overlapping less -- unless it is
    anchored to an end.  A third relation overlapping the pair is what unlocks
    it in the interior, and those rules are in the table.
    """
    def maxOverlapOf(relations):
        return ov.maxOverlap(_asRelLengths(relations))

    reducing = []
    for rule in lm.ALL_MOVES:
        before, after = rule[1], rule[2]
        if len(before) != 2 or maxOverlapOf(before) < 2:
            continue
        if maxOverlapOf(after) < maxOverlapOf(before):
            reducing.append(rule)
    assert reducing, "the end collapse must be among them"
    assert all(lm.anchorOf(rule) is not None for rule in reducing), [
        lm.formatMove(rule) for rule in reducing if lm.anchorOf(rule) is None]


def _asRelLengths(relations):
    """A window's relations as a relation-length row long enough to hold them."""
    width = max(start + arrows for start, arrows in relations) + 1
    row = [0] * width
    for start, arrows in relations:
        row[start] = arrows
    return row


@pytest.mark.slow
@pytest.mark.parametrize("rule", lm.ANCHORED_MOVES, ids=str)
def test_each_anchored_rule_holds_at_its_end(rule):
    """Re-run the verification the anchored table's entries were admitted by.

    A rule here claims one window position per length, so the lengths still have
    to follow the window -- R-009's lesson does not stop applying because there
    is only one position to check.
    """
    width = rule[0]
    lengths = range(width + 1, width + 5) if width + 4 <= 11 else range(width + 1, width + 3)
    confirmed, failures = lm.verifyMove(rule, lengths)
    assert failures == []
    assert confirmed > 0, "checked at lengths {0} for a window of {1} arrows".format(
        list(lengths), width)


def test_the_end_moves_are_what_carries_the_overlapping_rows():
    """The anchored half of the table is worth more than the floating half.

    Not a tautology and not a close thing: at n = 8 the floating rules place 13
    LNAs above the almost separate line and the anchored ones place another 101,
    and at n = 6 the two together leave nothing for a search at all.
    """
    withEnds = ov.coverage(6)
    assert len(withEnds['uncovered']) == 0
    floating = ov.coverage(8, lm.VERIFIED_MOVES)
    both = ov.coverage(8)
    assert len(floating['covered']) == 246
    assert len(both['covered']) == 347
