"""Mutation sequence shortcuts between LNAs.

A rule here is a *local rewrite*: a window of the quiver, what the relations
inside it become, and the mutation sequence that does it.  Nothing goes in
lnaMoves.VERIFIED_MOVES without verifyMove checking it against the mutation
engine at every window position of every LNA over a range of lengths, because
plausible-looking rules are often wrong -- "a lone relation of length 2 may be
deleted", stated on a two-arrow window, holds 63 times and fails 130 times.
"""

import itertools

import pytest
import sympy

from quivermutation import lnaMoves as lm
from quivermutation import nakayama as nk
import quivermutation as qm
from helpers import quiet


# ---------------------------------------------------------------------------
# The two moves this started from
# ---------------------------------------------------------------------------

PAIR_SLIDE_RIGHT = (5, ((0, 3), (1, 3)), ((1, 3), (2, 3)), (-5, -5))
PAIR_SLIDE_LEFT = (5, ((1, 3), (2, 3)), ((0, 3), (1, 3)), (2, 2))
MEETING_POINT_LEFT = (5, ((0, 3), (1, 3), (3, 2)), ((0, 2), (1, 3), (2, 3)), (2, 2))


def test_the_pair_slide_is_in_the_table_both_ways():
    """Two relations of equal length starting at consecutive vertices slide
    along the quiver: two right mutations at the first relation's source move
    them one arrow left, two left mutations at the second relation's target move
    them one arrow right."""
    assert PAIR_SLIDE_RIGHT in lm.VERIFIED_MOVES
    assert PAIR_SLIDE_LEFT in lm.VERIFIED_MOVES


def test_the_meeting_point_move_is_in_the_table():
    """A relation stays put while the meeting point of the two relations around
    it moves one arrow left, under two right mutations.

    In the form it was described: a relation of length r from v to v+r, one from
    v-1 to v+r-k and one from v+r-k to v+r+1; afterwards the middle relation is
    unchanged and the other two meet at v+r-k-1.  With r = 3 and k = 1 that is
    this entry, and the two mutations are at v rather than v-1.
    """
    assert MEETING_POINT_LEFT in lm.VERIFIED_MOVES

    width, before, after, sequence = MEETING_POINT_LEFT
    # The middle relation is the one that does not move.
    assert (1, 3) in before and (1, 3) in after
    # The other two meet one arrow earlier: 0+3 = 3 becomes 0+2 = 2.
    assert {(0, 3), (3, 2)} <= set(before)
    assert {(0, 2), (2, 3)} <= set(after)


def test_the_pair_slide_walks_a_pair_along_the_quiver_and_off_each_end():
    """A_10 with a pair of length-3 relations at the far left reaches every
    position of that pair -- and, at either end, loses one of the two.

    The pair slide alone gives the six positions.  What the last two members are
    is the point of F-021: an isolated pair overlapping in two arrows cannot be
    pulled apart anywhere in the interior, but against an end it collapses under
    two mutations at the end vertex.  So the orbit of a heavily overlapping LNA
    reaches an almost separate one exactly by walking to an end, which is why
    those rules matter out of all proportion to their two positions.
    """
    positions = {
        "33000000", "03300000", "00330000", "00033000", "00003300", "00000330",
    }
    orbit = lm.closureUnderMoves(10, [3, 3, 0, 0, 0, 0, 0, 0])
    assert set(orbit) == positions | {"30000000", "00000030"}
    assert set(lm.closureUnderMoves(10, [3, 3, 0, 0, 0, 0, 0, 0],
                                    rules = lm.VERIFIED_MOVES)) == positions


# ---------------------------------------------------------------------------
# The rule table
# ---------------------------------------------------------------------------

def test_the_table_is_not_empty_and_has_no_duplicates():
    assert len(lm.VERIFIED_MOVES) >= 15
    assert len(set(lm.VERIFIED_MOVES)) == len(lm.VERIFIED_MOVES)


@pytest.mark.parametrize("description", lm.VERIFIED_MOVES, ids=str)
def test_each_rule_is_well_formed(description):
    width, before, after, sequence = description
    assert width >= 1
    assert before != after
    assert sequence
    for relations in (before, after):
        for start, arrows in relations:
            assert arrows >= 2
            assert 0 <= start
            assert start + arrows <= width + 1
    for vertex in sequence:
        assert 1 <= abs(vertex) <= width + 1


def lengthsToCheck(width):
    """The lengths a rule of this window width can be checked at.

    A window of `width` arrows needs a quiver of at least `width + 1` vertices
    to sit in at all, so a *fixed* range of lengths cannot check every rule in
    the table: the pair slide is one statement for every relation length and its
    window grows with that length (F-013), and the widest rules in the table --
    window 11, relation length 9 -- do not fit in a quiver of length 8.

    This test used to ask for lengths 5 to 8 whatever the rule, which is
    `width + 1 .. width + 4` for a window of 4 and nothing at all for a window
    of 8 or more. The eight widest rules therefore got zero confirmations and
    the test failed on them -- not because they are wrong, but because it was
    verifying them where they cannot occur.

    Four lengths where that is affordable and two where it is not: verification
    enumerates every admissible LNA of the length, and there are 208012 of them
    at length 13.
    """
    if width + 4 <= 10:
        return range(width + 1, width + 5)
    return range(width + 1, width + 3)


@pytest.mark.slow
@pytest.mark.parametrize("description", lm.VERIFIED_MOVES, ids=str)
def test_each_rule_holds_wherever_it_applies(description):
    """Re-run the verification the table's entries were admitted by."""
    width = description[0]
    confirmed, failures = lm.verifyMove(description, lengthsToCheck(width))
    assert confirmed > 0, "checked at lengths {0} where the window is {1} arrows".format(
        list(lengthsToCheck(width)), width)
    assert failures == []


def test_the_lengths_a_rule_is_checked_at_can_contain_its_window():
    """The bug the fix above is for, stated as a check of its own.

    Every rule in the table must be verified at lengths that can hold its
    window, or `confirmed > 0` is asserting about nothing.
    """
    for description in lm.VERIFIED_MOVES:
        width = description[0]
        lengths = list(lengthsToCheck(width))
        assert lengths, description
        assert min(lengths) >= width + 1, (description, lengths)


def test_a_plausible_rule_that_is_actually_false_is_rejected():
    """Deleting a lone length-2 relation, stated on a two-arrow window.

    Relations of length 2 do not change the derived equivalence class, so this
    looks safe, but on a two-arrow window the rewrite does not see enough of the
    quiver and is wrong more often than right.  The three-arrow version that
    slides such a relation *is* in the table.
    """
    tooNarrow = (2, ((0, 2),), (), (1,))
    confirmed, failures = lm.verifyMove(tooNarrow, range(5, 8))
    assert failures, "this rule is false and must not verify"
    assert tooNarrow not in lm.VERIFIED_MOVES


def test_a_rule_built_from_an_inadmissible_mutation_is_rejected():
    """Landing on the predicted LNA is not enough.

    A mutation outside the procedure's admissibility condition still produces a
    quiver, just not a derived equivalent one, so a rewrite can hit exactly the
    predicted relation lengths and still be false.  This shortening of a
    length-3 relation to a length-2 one does precisely that: it was discovered,
    it always produces the predicted result, and every application of it is an
    illegal mutation.
    """
    shortening = (3, ((0, 3),), ((0, 2),), (3, -3))
    confirmed, failures = lm.verifyMove(shortening, range(5, 8))
    assert confirmed == 0
    assert {reason for *_rest, reason in failures} == {"illegal mutation"}
    assert shortening not in lm.VERIFIED_MOVES

    # Without the admissibility and Coxeter checks it would look perfect.
    lenient, lenientFailures = lm.verifyMove(shortening, range(5, 8), checkCoxeter=False)
    assert lenientFailures  # still caught, by the admissibility half


def test_every_rule_uses_only_admissible_mutations():
    """Each entry's sequence must be a chain of genuine tilting mutations."""
    for description in lm.VERIFIED_MOVES:
        width, before, _after, _offsets = description
        placed = False
        for length in range(width + 2, width + 5):
            for windowStart in range(1, length):
                relLengths = [0] * (length - 2)
                fits = True
                for start, arrows in before:
                    position = windowStart + start - 1
                    if position < 0 or position >= len(relLengths) or windowStart + start + arrows > length:
                        fits = False
                        break
                    relLengths[position] = arrows
                if not fits or not lm.isAdmissible(length, relLengths):
                    continue
                if not lm.matchesAt(length, relLengths, description, windowStart):
                    continue
                applied = lm.applyAt(length, relLengths, description, windowStart)
                if applied is None:
                    continue
                _predicted, sequence = applied
                assert lm.isLegalSequence(
                    nk.LinearNakayamaAlgebra(length, relLengths), sequence), description
                placed = True
            if placed:
                break
        assert placed, f"could not place {description} anywhere"


# ---------------------------------------------------------------------------
# Orbits
# ---------------------------------------------------------------------------

def test_an_orbit_records_a_mutation_path_that_works():
    for length, rels in [(8, "033000"), (9, "0022022"), (10, "22003300")]:
        relLengths = [int(c) for c in rels]
        for name, (sequence, numbering) in lm.closureUnderMoves(length, relLengths).items():
            if not sequence:
                continue
            mutated = quiet(qm.quiverMutationAtVertices,
                            lm._copy(nk.LinearNakayamaAlgebra(length, relLengths)),
                            list(sequence))
            assert qm.className(lm.asRelLengths(mutated, length)) == name, (rels, name, sequence)


@pytest.mark.slow
@pytest.mark.parametrize("length", [5, 6, 7])
def test_every_orbit_is_inside_one_derived_equivalence_class(length):
    """The Coxeter polynomial must be constant on an orbit -- every member is
    reached by an actual mutation sequence, so they are all derived equivalent."""
    for relLengths in itertools.product(range(0, length), repeat=length - 2):
        relLengths = list(relLengths)
        if not lm.isAdmissible(length, relLengths):
            continue
        orbit = lm.closureUnderMoves(length, relLengths)
        if len(orbit) < 2:
            continue
        expected = sympy.expand(
            quiet(qm.coxeterPoly, nk.LinearNakayamaAlgebra(length, relLengths)).as_expr())
        for name in orbit:
            got = sympy.expand(quiet(
                qm.coxeterPoly,
                nk.LinearNakayamaAlgebra(length, [int(c) for c in name])).as_expr())
            assert got == expected, (qm.className(relLengths), name)


def test_admissibility_matches_the_enumeration():
    """isAdmissible must accept exactly the LNAs generateAllPossibleLineRelations
    produces."""
    for length in range(3, 8):
        enumerated = {
            qm.className(qm.relationStringToLineRelLengths(length, qm.relSetToString(relSet)))
            for relSet in qm.generateAllPossibleLineRelations(length)
        }
        byPredicate = {
            "".join(str(n) for n in candidate)
            for candidate in itertools.product(range(0, length), repeat=length - 2)
            if lm.isAdmissible(length, list(candidate))
        }
        assert enumerated == byPredicate, length


# -- is a move a *local* rewrite? (research H-009) -------------------------

def coversFromLeft(relLengths, windowStart):
    """Whether a relation reaches the window's first arrow from outside it."""
    for start, arrows in lm.relationsOf(relLengths):
        if start < windowStart and start + arrows - 1 >= windowStart:
            return True
    return False


def matchesLocally(length, relLengths, description, windowStart):
    """`lm.matchesAt` decided from a bounded neighbourhood.

    `matchesAt` scans the whole relation-length row. This reads only

    * the window's own cells -- a relation starting at an offset inside the
      window and staying inside it has at most `width` arrows, so what a cell
      can say is bounded by the width; and
    * one bit, whether any relation covers the window's first arrow having
      started strictly before it.

    That bit is the part the per-vertex encoding cannot supply locally, since a
    relation of unbounded length reaches arbitrarily far to the right, and the
    part a per-arrow encoding carries for free.
    """
    width, before, _after, _offsets = description
    windowEnd = windowStart + width - 1
    if windowStart < 1 or windowEnd > length - 1:
        return False
    if coversFromLeft(relLengths, windowStart):
        return False
    inside = []
    for offset in range(width):
        cell = windowStart + offset
        arrows = relLengths[cell - 1] if cell - 1 < len(relLengths) else 0
        if not arrows:
            continue
        if cell + arrows - 1 > windowEnd:
            return False
        inside.append((offset, arrows))
    return tuple(sorted(inside)) == before


@pytest.mark.parametrize("length", [5, 6, 7])
def test_whether_a_move_applies_is_a_local_condition(length):
    """Every rule in the table, every LNA, every window position.

    This is the check research H-009 turns on. The suspicion there is that the
    move table is a one-dimensional cellular automaton -- a local rewrite on the
    row of relation lengths -- and the caveat is that a rule whose applicability
    depended on the row far away would not be local at all, whatever it looked
    like.

    It does not: the answer comes out of the window's cells plus one bit about
    what reaches into its left edge. Which also says what the right encoding
    is -- per arrow, carrying "covered", rather than per vertex carrying a
    relation length -- because that bit is exactly what a per-arrow row has and
    a per-vertex row does not.
    """
    matched = 0
    for relLengths in itertools.product(range(0, length), repeat = max(0, length - 2)):
        relLengths = list(relLengths)
        if not lm.isAdmissible(length, relLengths):
            continue
        for description in lm.VERIFIED_MOVES:
            for windowStart in range(1, length):
                actual = lm.matchesAt(length, relLengths, description, windowStart)
                assert actual == matchesLocally(length, relLengths, description, windowStart), (
                    length, relLengths, description, windowStart)
                matched += bool(actual)
    assert matched > 0, "no rule matched anywhere, so this checked nothing"


@pytest.mark.slow
@pytest.mark.parametrize("length", [8, 9])
def test_whether_a_move_applies_is_a_local_condition_further_out(length):
    test_whether_a_move_applies_is_a_local_condition(length)


# -- families whose mutation count grows with the parameter ----------------

@pytest.mark.parametrize("distance", [1, 2, 3, 4, 5])
def test_the_short_relation_slide_takes_one_mutation_per_arrow(distance):
    """The family's shape, without running the engine.

    A lone relation of two arrows travels `d` arrows under `d` mutations, so the
    window is `d + 2` wide and the sequence has `d` entries -- which is why
    discovery bounded at three mutations found only d = 1, 2, 3.  The pair slide
    is the contrast: two mutations whatever its parameter.
    """
    rules = {(rule[1], rule[2]): rule for rule in lm.shortRelationSlideRules(5)}
    right = rules[(((0, 2),), ((distance, 2),))]
    left = rules[(((distance, 2),), ((0, 2),))]
    for rule in (right, left):
        assert rule[0] == distance + 2
        assert len(rule[3]) == distance
    assert all(vertex < 0 for vertex in right[3])
    assert all(vertex > 0 for vertex in left[3])
    # Every mutation is inside the window it claims.
    for rule in (right, left):
        assert all(2 <= abs(vertex) <= rule[0] + 1 for vertex in rule[3])


def test_discovery_found_an_initial_segment_of_the_family():
    """Generating the family must not contradict the table it extends.

    What discovery can find is an initial segment: `d` needs `d` mutations, so a
    search bounded at `maxSteps` sees `d <= maxSteps` and no more.  The listed
    members must therefore be the small ones, with no gaps and nothing beyond.
    """
    listed = set(lm.DISCOVERED_MOVES)
    found = sorted(distance for distance in range(1, 8)
                   if all(rule in listed
                          for rule in lm.shortRelationSlideRules(distance)[-2:]))
    assert found == list(range(1, len(found) + 1)), "an initial segment, no gaps"
    assert len(found) >= 2
    assert set(lm.shortRelationSlideRules()) <= set(lm.VERIFIED_MOVES)


def test_both_families_are_in_the_table_and_nothing_is_duplicated():
    assert len(lm.VERIFIED_MOVES) == len(set(lm.VERIFIED_MOVES))
    for rule in lm.pairSlideRules() + lm.shortRelationSlideRules():
        assert rule in lm.VERIFIED_MOVES


@pytest.mark.parametrize("generator, atDistanceOne", [
    (lm.shortRelationSlideRules, 1),
    (lm.trailingRelationWalkRules, 2),
    (lm.spreadingPairRules, 2),
])
def test_each_slide_family_pays_one_mutation_per_arrow(generator, atDistanceOne):
    """All three families have the same shape: the d-th member costs d mutations
    for the travel, plus one per companion relation that has to be displaced.

    That is what makes them invisible to a bounded search past their first
    members, and what makes generating them worth more than searching deeper.
    """
    rules = generator()
    perDistance = {}
    for rule in rules:
        # A family's rules are one per d, except the lone slide, which has both
        # directions at each d.
        perDistance.setdefault(rule[0], []).append(rule)
    widths = sorted(perDistance)
    for index, width in enumerate(widths):
        for rule in perDistance[width]:
            assert len(rule[3]) == atDistanceOne + index, (rule, index)
        assert width == widths[0] + index


def test_the_inverse_of_a_rule_is_not_free():
    """E-019: reverse, negate, shift toward zero inverts a slide and little else.

    Kept as a test because the transform is tempting -- it turns the lone slide's
    right rule into exactly the left rule the table lists -- and it is wrong for
    the table at large, so the counterexample is worth pinning.
    """
    def inverseOf(rule):
        width, before, after, sequence = rule
        flipped = [-vertex for vertex in reversed(sequence)]
        return (width, after, before,
                tuple((v - 1) if v > 0 else (v + 1) for v in flipped))

    slide = lm.shortRelationSlideRules(3)
    for right, left in zip(slide[::2], slide[1::2]):
        assert inverseOf(right) == left

    # The spreading pair mixes directions, and its "inverse" is not one.
    spreading = lm.spreadingPairRules(2)[0]
    candidate = inverseOf(spreading)
    assert candidate not in lm.VERIFIED_MOVES
    confirmed, failures = lm.verifyMove(candidate, [candidate[0] + 1])
    assert failures or confirmed == 0
