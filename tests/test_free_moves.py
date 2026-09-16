"""Relations of two arrows, which cost no mutation at all.

`corollary:lengthtworelations` of arXiv:2310.08346 says a relation of two arrows
leaves the derived equivalence class alone.  These tests pin the three things
that follow (research F-028): such a relation never overlaps a neighbour in more
than one arrow, so deleting them never changes whether the quipu theorem names
an LNA; the reduced space is the LNAs of one fewer vertex; and adding the move
to the rule table does far more for the coverage than the table itself does.

The Coxeter polynomial is a derived invariant, so it has to be blind to these
relations, and the quipu theorem has to be blind to them too.  Both are checked
here -- two independent routes to the corollary, neither of them its proof.
"""

import pytest

from quivermutation import freeMoves as fm
from quivermutation import invariants
from quivermutation import lnaMoves as lm
from quivermutation import nakayama as nk
from quivermutation import overlap as ov
from quivermutation import quipuForms as qf


# ---------------------------------------------------------------------------
# The move itself
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("relLengths, stripped", [
    ((0, 0, 0, 0), (0, 0, 0, 0)),
    ((2, 0, 0, 0), (0, 0, 0, 0)),
    ((3, 0, 2, 0), (3, 0, 0, 0)),
    ((2, 0, 2, 0), (0, 0, 0, 0)),
    ((4, 0, 0, 3), (4, 0, 0, 3)),
])
def test_stripping_deletes_exactly_the_two_arrow_relations(relLengths, stripped):
    assert fm.stripLengthTwo(relLengths) == stripped


def test_the_two_arrow_relations_are_reported_by_start_vertex():
    assert fm.lengthTwoRelations((2, 0, 3, 2, 0)) == (1, 4)
    assert fm.lengthTwoRelations((3, 0, 4, 0, 0)) == ()


def test_a_reduced_lna_is_one_with_no_two_arrow_relation():
    assert fm.isReduced((3, 0, 4, 0))
    assert not fm.isReduced((3, 0, 2, 0))


def test_shortening_refuses_an_lna_that_is_not_reduced():
    with pytest.raises(ValueError):
        fm.shortenByOneArrow((3, 0, 2, 0))


# ---------------------------------------------------------------------------
# The overlap lemma: a two-arrow relation is never heavily overlapping
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("length", [5, 6, 7, 8, 9, 10])
def test_a_two_arrow_relation_never_shares_more_than_one_arrow(length):
    """Admissibility alone forces it, so no LNA of any length can break it."""
    for relLengths in nk.allRelationLengths(length):
        relations = lm.relationsOf(list(relLengths))
        profile = ov.overlapProfile(list(relLengths))
        for overlapValue, (first, second) in zip(profile, zip(relations, relations[1:])):
            if first[1] == 2 or second[1] == 2:
                assert overlapValue <= 1, (relLengths, first, second)


@pytest.mark.parametrize("length", [5, 6, 7, 8, 9, 10])
def test_stripping_never_changes_whether_the_theorem_names_the_lna(length):
    """So the move never names anything on its own -- all it does is bridge."""
    for relLengths in nk.allRelationLengths(length):
        stripped = fm.stripLengthTwo(relLengths)
        assert (ov.isAlmostSeparate(length, relLengths)
                == ov.isAlmostSeparate(length, stripped)), relLengths


@pytest.mark.parametrize("length", [5, 6, 7, 8, 9, 10])
def test_the_quipu_theorem_gives_the_same_name_after_stripping(length):
    """The theorem the classification is seeded from is blind to these relations."""
    for relLengths in nk.allRelationLengths(length):
        stripped = fm.stripLengthTwo(relLengths)
        assert (qf.quipuForAlmostSeparateLNA(length, list(relLengths))
                == qf.quipuForAlmostSeparateLNA(length, list(stripped))), relLengths


# ---------------------------------------------------------------------------
# The reduced space is one vertex smaller
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("length", [4, 5, 6, 7, 8, 9, 10, 11])
def test_shortening_is_a_bijection_onto_the_lnas_of_one_fewer_vertex(length):
    reduced = fm.reducedForms(length)
    images = sorted(fm.shortenByOneArrow(form) for form in reduced)
    assert images == sorted(nk.allRelationLengths(length - 1))


@pytest.mark.parametrize("length, count", [
    (6, 14), (7, 42), (8, 132), (9, 429), (10, 1430), (11, 4862),
])
def test_the_reduced_forms_are_counted_by_the_previous_catalan_number(length, count):
    assert len(fm.reducedForms(length)) == count
    assert count == len(nk.allRelationLengths(length - 1))


# ---------------------------------------------------------------------------
# What it buys the search
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("length, orbitsWithout, orbitsWith, leftWithout, leftWith", [
    (7, 20, 10, 0, 0),
    (8, 69, 22, 10, 1),
    (9, 380, 91, 222, 53),
])
def test_the_free_move_merges_far_more_than_the_rule_table_does(
        length, orbitsWithout, orbitsWith, leftWithout, leftWith):
    """The measurement of F-028, pinned so a change to the table shows up here."""
    without = ov.coverage(length)
    with_ = fm.coverage(length)
    assert len(without['orbits']) == orbitsWithout
    assert len(with_['orbits']) == orbitsWith
    assert len(without['uncovered']) == leftWithout
    assert len(with_['uncovered']) == leftWith


def test_with_the_edge_moves_as_well_a_8_needs_no_search_at_all():
    """The headline of F-028 and F-029 together, and `overlaps.py 8 --free --edges`."""
    result = fm.coverage(8, free = True, edges = True)
    assert len(result['uncovered']) == 0
    assert len(result['covered']) == len(result['lnas']) == 429
    assert len(result['orbits']) == 21


@pytest.mark.parametrize("length, orbits, left", [(9, 77, 37), (8, 21, 0)])
def test_the_free_move_and_the_edge_moves_together_are_pinned(length, orbits, left):
    result = fm.coverage(length, free = True, edges = True)
    assert len(result['orbits']) == orbits
    assert len(result['uncovered']) == left


@pytest.mark.parametrize("length", [6, 7, 8])
def test_the_free_orbits_are_coarser_than_the_move_orbits(length):
    """Every move orbit lies inside a free orbit -- the move only ever joins."""
    _, moveOrbits = ov.moveOrbits(length)
    _, freeOrbits = fm.derivedOrbits(length)
    within = {member: root for root, members in freeOrbits.items() for member in members}
    for members in moveOrbits.values():
        assert len({within[member] for member in members}) == 1


def test_the_table_deletes_a_two_arrow_relation_only_against_an_end():
    """Which is why the free move is worth having: it does not need the end.

    A lone two-arrow relation at the sink of A_8 is deleted by the rules; the
    same relation four arrows in is only ever slid, and one with a longer
    relation for company is not deleted even at the sink.
    """
    atTheEnd = lm.rewritesOf(8, [0, 0, 0, 0, 0, 2], lm.ALL_MOVES)
    assert (0, 0, 0, 0, 0, 0) in atTheEnd
    inTheMiddle = lm.rewritesOf(8, [0, 0, 0, 2, 0, 0], lm.ALL_MOVES)
    assert (0, 0, 0, 0, 0, 0) not in inTheMiddle
    withCompany = lm.rewritesOf(8, [3, 0, 0, 0, 2, 0], lm.ALL_MOVES)
    assert (3, 0, 0, 0, 0, 0) not in withCompany
    assert fm.stripLengthTwo((3, 0, 0, 0, 2, 0)) == (3, 0, 0, 0, 0, 0)


# ---------------------------------------------------------------------------
# The derived invariant, which has to agree
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("length", [5, 6, 7, 8, 9])
def test_stripping_keeps_the_coxeter_polynomial(length):
    """Necessary for the corollary, since the polynomial is a derived invariant."""
    cache = {}

    def poly(relLengths):
        if relLengths not in cache:
            cache[relLengths] = invariants.coxeterPoly(
                nk.LinearNakayamaAlgebra(length, list(relLengths)))
        return cache[relLengths]

    for relLengths in nk.allRelationLengths(length):
        assert poly(relLengths) == poly(fm.stripLengthTwo(relLengths)), relLengths
