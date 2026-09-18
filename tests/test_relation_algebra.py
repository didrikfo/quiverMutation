"""Relations as linear combinations of paths, with exact ideal membership.

The set-of-paths model cannot express a coefficient, so it decides "is this path
zero" by looking for a zero relation sitting inside it, and propagates equalities
through commutativity relations by a partial closure.  relationAlgebra decides it
by linear algebra over the ideal instead, which is exact.
"""

import pytest

from quivermutation import pathAlgebra as pac
import quivermutation as qm
from quivermutation import nakayama as nk
from quivermutation import relationAlgebra as ra
from helpers import quiet


def algebra(arrows, rels=()):
    pa = pac.PathAlgebra()
    pa.add_arrows_from([list(a) for a in arrows])
    pa.add_rels_from([[list(p) for p in rel] for rel in rels])
    return pa


GRID_ARROWS = [(1, 2), (2, 3), (1, 4), (2, 5), (3, 6), (4, 5), (5, 6)]
GRID_SQUARES = [[(1, 2, 5), (1, 4, 5)], [(2, 3, 6), (2, 5, 6)]]


# -- the combination type ------------------------------------------------

def test_a_single_arrow_path_is_not_mistaken_for_a_coefficient_pair():
    """The path along one arrow is a two-element tuple, like a (path, coeff) pair.

    Telling them apart by length is wrong; the first element decides, since a
    pair's is a path and a bare path's is a vertex.  A relation containing a
    path of length one is exactly what the mutation procedure produces before
    reduction cancels it, so this comes up immediately in practice.
    """
    assert ra.combination([[1, 2]]) == {(1, 2): 1}
    assert ra.combination([(1, 2)]) == {(1, 2): 1}
    assert ra.combination([((1, 2), 3)]) == {(1, 2): 3}
    assert ra.combination({(1, 2): 3}) == {(1, 2): 3}
    assert ra.asSum([[1, 2], [1, 3, 2]]) == {(1, 2): 1, (1, 3, 2): 1}


def test_an_inadmissible_relation_reaches_the_cartan_matrix_intact():
    """A relation containing a single arrow says that arrow equals the rest.

    quiverMutationAtVertex produces these routinely and reducePathAlgebra
    cancels them, so the exact Cartan matrix has to cope with one in place.
    """
    pa = algebra([(1, 2), (1, 3), (3, 2)], [[(1, 2), (1, 3, 2)]])
    rels = [ra.fromPathSet(rel) for rel in pa.rels]
    # The two paths from 1 to 2 are identified, leaving one.
    assert ra.homDimension(pa.quiver, rels, 1, 2) == 1


def test_combination_sums_repeats_and_drops_cancellations():
    assert ra.combination([((1, 2), 1), ((1, 2), 2)]) == {(1, 2): 3}
    assert ra.combination([((1, 2), 1), ((1, 2), -1)]) == {}
    assert ra.combination([[1, 2, 3], [1, 4, 3]]) == {(1, 2, 3): 1, (1, 4, 3): 1}


def test_signs_alone_do_not_close_under_addition():
    """p + q + r = 0 and p - q = 0 give 2p + r = 0, which needs a 2.

    This is why the coefficients are integers rather than just signs: the
    arithmetic leaves {-1, 0, +1} on the first step that combines two relations.
    """
    p, q, r = (1, 2, 5), (1, 3, 5), (1, 4, 5)
    sum_of_three = ra.combination([(p, 1), (q, 1), (r, 1)])
    difference = ra.combination([(p, 1), (q, -1)])
    combined = ra.add(sum_of_three, difference)
    assert combined == {p: 2, r: 1}
    assert set(combined.values()) - {-1, 0, 1}


def test_left_and_right_divide_match_the_papers_notation():
    """r/alpha keeps the paths starting with alpha and drops it; the mirror for
    the last arrow."""
    r = ra.combination([((1, 2, 5), 1), ((1, 3, 5), -1)])
    assert ra.leftDivide(r, 2) == {(2, 5): 1}
    assert ra.leftDivide(r, 3) == {(3, 5): -1}
    assert ra.leftDivide(r, 9) == {}       # no path of r begins with that arrow
    assert ra.rightDivide(r, 2) == {(1, 2): 1}
    assert ra.rightDivide(r, 3) == {(1, 3): -1}


def test_pre_and_post_composition():
    r = ra.combination([((2, 5), 1), ((2, 6), 1)])
    assert ra.preCompose(r, [1, 2]) == {(1, 2, 5): 1, (1, 2, 6): 1}
    assert ra.postCompose(ra.combination([((1, 2), 1)]), [2, 7]) == {(1, 2, 7): 1}


def test_a_combination_must_have_one_source_and_one_target():
    mixed = ra.combination([((1, 2), 1), ((1, 3, 4), 1)])
    with pytest.raises(ValueError):
        ra.isInIdeal(algebra([(1, 2), (1, 3), (3, 4)]).quiver, [], mixed)


# -- exact ideal membership ----------------------------------------------

def test_a_two_path_relation_is_read_as_a_difference():
    """A relation written as two paths means commutativity, so p - q, not p + q.

    Reading it as a sum is not harmless.  Three commutativity relations among
    three parallel paths p, q, r give p = -q, r = -q and p + r = -2q, which
    forces q = 0 and collapses a Hom space that should be one-dimensional.  This
    turned up when checking that reducePathAlgebra preserves the Cartan matrix:
    nine of 38095 reductions appeared to change it, and every one was this sign
    reading rather than a fault in the reduction.
    """
    p, q, r = (1, 2, 3, 7), (1, 2, 6, 7), (1, 5, 6, 7)
    assert ra.fromPathSet([p, q]) == {p: 1, q: -1}
    assert ra.asSum([p, q]) == {p: 1, q: 1}
    assert ra.fromPathSet([p]) == {p: 1}
    assert ra.fromPathSet([p, q, r]) == {p: 1, q: 1, r: 1}

    pa = algebra(
        [(1, 2), (1, 5), (2, 3), (2, 6), (3, 7), (5, 6), (6, 7)],
        [[p, q], [q, r], [p, r]],
    )
    asDifferences = [ra.fromPathSet(rel) for rel in pa.rels]
    asSums = [ra.asSum(rel) for rel in pa.rels]
    assert ra.homDimension(pa.quiver, asDifferences, 1, 7) == 1
    assert ra.homDimension(pa.quiver, asSums, 1, 7) == 0


def test_a_zero_relation_propagates_through_chained_commutative_squares():
    """The case the set-of-paths model gets wrong.

    In the 2x2 grid the two commutativity relations make all three paths from 1
    to 6 equal, so killing one kills all three.  The old pathHasZeroRel only
    recognises the one that literally contains the zero relation.
    """
    pa = algebra(GRID_ARROWS, GRID_SQUARES + [[(1, 2, 3, 6)]])
    rels = [ra.fromPathSet(rel) for rel in pa.rels]

    for path in [(1, 2, 3, 6), (1, 2, 5, 6), (1, 4, 5, 6)]:
        assert ra.isInIdeal(pa.quiver, rels, ra.combination([path])), path

    # What the old predicate sees, for the record.
    assert qm.pathHasZeroRel([1, 2, 3, 6], pa.rels)
    assert not qm.pathHasZeroRel([1, 2, 5, 6], pa.rels)
    assert not qm.pathHasZeroRel([1, 4, 5, 6], pa.rels)


def test_without_the_zero_relation_the_grid_has_one_path():
    pa = algebra(GRID_ARROWS, GRID_SQUARES)
    rels = [ra.fromPathSet(rel) for rel in pa.rels]
    assert ra.homDimension(pa.quiver, rels, 1, 6) == 1
    assert not ra.isInIdeal(pa.quiver, rels, ra.combination([(1, 4, 5, 6)]))
    # The three paths are equal, so any difference of two of them is in the ideal.
    difference = ra.combination([((1, 4, 5, 6), 1), ((1, 2, 3, 6), -1)])
    assert ra.isInIdeal(pa.quiver, rels, difference)


def test_the_zero_relation_collapses_the_whole_hom_space():
    pa = algebra(GRID_ARROWS, GRID_SQUARES + [[(1, 2, 3, 6)]])
    rels = [ra.fromPathSet(rel) for rel in pa.rels]
    assert ra.homDimension(pa.quiver, rels, 1, 6) == 0


def test_hom_dimension_of_a_line_with_a_relation():
    """A_5 with 1->2->3->4 = 0: one path between any two vertices, except that
    the killed one and anything containing it are gone."""
    pa = algebra([(1, 2), (2, 3), (3, 4), (4, 5)], [[(1, 2, 3, 4)]])
    rels = [ra.fromPathSet(rel) for rel in pa.rels]
    assert ra.homDimension(pa.quiver, rels, 1, 3) == 1
    assert ra.homDimension(pa.quiver, rels, 1, 4) == 0
    assert ra.homDimension(pa.quiver, rels, 1, 5) == 0
    assert ra.homDimension(pa.quiver, rels, 2, 5) == 1
    assert ra.homDimension(pa.quiver, rels, 1, 1) == 1   # the trivial path
    assert ra.homDimension(pa.quiver, rels, 5, 1) == 0   # nothing goes backwards


def test_the_trivial_path_is_never_in_an_admissible_ideal():
    pa = algebra([(1, 2), (2, 3)], [[(1, 2, 3)]])
    rels = [ra.fromPathSet(rel) for rel in pa.rels]
    for vertex in (1, 2, 3):
        assert ra.homDimension(pa.quiver, rels, vertex, vertex) == 1


# -- the cheap count, and the chain the old one missed -------------------

def test_the_cheap_count_closes_the_commutativity_relations_transitively():
    """A path is zero when a *chain* of identifications reaches a zero path.

    The quiver reached from `40030` at n = 7 by `[1, 4, 2, 5, 2]`.  Its four
    paths `1 ~~> 7` are joined in one class by two commutativity relations, and
    one member of that class, `4 -> 2 -> 7` extended, is a zero relation -- so
    every one of them is zero and the entry is 0.  Reaching it needs three
    substitutions in a row:

        1,4,2,5,7  ~  1,6,2,5,7   (by 1,4,2 = 1,6,2)
                   ~  1,6,2,7     (by 6,2,5,7 = 6,2,7)
                   ~  1,4,2,7     (by 1,4,2 = 1,6,2)  = 0

    `paths.numberOfPathsUpToRels` applies each relation at most once per pass and
    compares canonical forms, so it stopped short of that and counted the class
    as nonzero -- which moved the Coxeter key and is one of the four "clean"
    wrong-key nodes F-038 recorded at n = 7 to depth 5.  It was a mis-count, not
    a bad mutation.  `arrowPaths.homDimensionByClosure` takes the full closure
    and agrees with the exact answer here.
    """
    from quivermutation import arrowPaths as ap
    from quivermutation import procedure as pr
    from quivermutation import invariants as inv
    from quivermutation import mutation, reduction

    algebra = nk.LinearNakayamaAlgebra(7, "40030")
    base = inv.coxeterKey(algebra)
    for vertex in [1, 4, 2, 5, 2]:
        algebra = quiet(reduction.reducePathAlgebra,
                        quiet(mutation.quiverMutationAtVertex, algebra, vertex))
    assert not algebra.hasParallelArrows()
    assert algebra.rels == [[[1, 4, 2], [1, 6, 2]], [[3, 6, 2]], [[4, 2, 7]],
                            [[6, 2, 5, 7], [6, 2, 7]]]

    relations = pr.relationsFrom(algebra)
    assert ap.homDimensionByClosure(algebra.quiver, relations, 1, 7) == 0
    assert ap.homDimension(algebra.quiver, relations, 1, 7) == 0
    assert quiet(qm.numberOfPathsUpToRels, algebra, 1, 7) == 1
    assert inv.coxeterKey(algebra) == base


# -- agreement with the existing Cartan matrix ---------------------------

@pytest.mark.parametrize("length", [3, 4, 5, 6])
def test_exact_and_heuristic_cartan_matrices_agree_on_every_lna(length):
    """Neither the published results nor the search depend on the difference.

    The heuristic numberOfPathsUpToRels is wrong in general -- the grid above
    shows it -- but on the LNAs themselves it is right, so replacing it does not
    move any Coxeter polynomial in the tables.
    """
    for relSet in qm.generateAllPossibleLineRelations(length):
        relLengths = [0] * (length - 2)
        for rel in relSet:
            relLengths[rel[0][0] - 1] = len(rel[0]) - 1
        pa = nk.LinearNakayamaAlgebra(length, relLengths)
        assert quiet(qm.cartanMatrix, pa) == quiet(ra.cartanMatrixExact, pa), relLengths


@pytest.mark.slow
@pytest.mark.parametrize("length, rels", [(5, "300"), (6, "3030"), (6, "2300"), (7, "22230")])
def test_the_two_cartan_matrices_agree_along_mutation_paths(length, rels):
    """Walk every legal mutation to depth 3 and compare at each step.

    `qm.cartanMatrix` counts arrow paths, exactly or by the cheap closure, and
    the two must agree everywhere an LNA walk goes.  `ra.cartanMatrixExact` is
    the *vertex-model* version and is compared as well, but only while no
    parallel pair has appeared -- past that it counts two paths as one, which is
    the whole point of `arrowPaths`.
    """
    from helpers import line_algebra
    from quivermutation import arrowPaths as ap
    from quivermutation import procedure as pr

    def walk(pa, depth):
        if depth == 0:
            return
        for vertex in pa.vertices():
            if not quiet(qm.mutationIsPossibleAtVertex, pa, vertex):
                continue
            mutated = quiet(qm.quiverMutationAtVertex, pa, vertex)
            if any(ap.isIllegalRelation(mutated.quiver, r)
                   for r in pr.relationsFrom(mutated)):
                continue
            mutated = quiet(qm.reducePathAlgebra, mutated)
            where = f"A_{length}_{rels}: {sorted(mutated.arrows())} {mutated.rels}"
            # The cheap count only where it is provably right: a relation of
            # three or more paths has no reading in it at all, and two-path
            # relations whose identifications close a cycle lose a dimension it
            # cannot see.  `isMonomial` is the condition `integerCartanMatrix`
            # itself uses.  F-039.
            if ap.isMonomial(pr.relationsFrom(mutated)):
                assert quiet(qm.cartanMatrix, mutated) == quiet(
                    qm.cartanMatrix, mutated, False), where
            if not mutated.hasParallelArrows():
                assert quiet(qm.cartanMatrix, mutated) == quiet(
                    ra.cartanMatrixExact, mutated), where
            walk(mutated, depth - 1)

    walk(line_algebra(length, rels), 3)


@pytest.mark.slow
def test_reduction_preserves_the_cartan_matrix():
    """reducePathAlgebra must present the same algebra it was given.

    The raw output of quiverMutationAtVertex contains inadmissible and redundant
    relations; reduction cancels and drops them, which changes the quiver but is
    supposed to leave the algebra alone.  The Cartan matrix is the sharpest cheap
    witness of that, since it is indexed by the vertices and those do not move.

    The full sweep -- every legal mutation of depth <= 3 out of every LNA of
    length 5 to 8, 38095 reductions -- comes back clean.  It is cut to lengths 5
    and 6 here to stay affordable.  It is also what turned up the sign reading
    fixed in fromPathSet: nine cases appeared to fail, all of them commutativity
    relations being read as sums.
    """
    from helpers import line_algebra
    from quivermutation import arrowPaths as ap
    from quivermutation import procedure as pr

    checked = 0

    def walk(pa, depth):
        nonlocal checked
        if depth == 0:
            return
        for vertex in pa.vertices():
            if not quiet(qm.mutationIsPossibleAtVertex, pa, vertex):
                continue
            raw = quiet(qm.quiverMutationAtVertex, pa, vertex)
            if any(ap.isIllegalRelation(raw.quiver, rel)
                   for rel in pr.relationsFrom(raw)):
                continue
            # The arrow-path Cartan matrix, not `ra.cartanMatrixExact`: the raw
            # output of a mutation is exactly where parallel arrows turn up, and
            # the vertex-model count is wrong there.
            before = quiet(qm.cartanMatrix, raw)
            reduced = quiet(qm.reducePathAlgebra, raw)
            after = quiet(qm.cartanMatrix, reduced)
            checked += 1
            assert before == after, (
                "reduction changed the algebra\n"
                f"  raw: {sorted((a[0], a[1]) for a in raw.arrows())} {raw.rels}\n"
                f"  red: {sorted((a[0], a[1]) for a in reduced.arrows())} {reduced.rels}"
            )
            # The cheap count must agree too, wherever it is taken -- which is
            # on a monomial ideal and nowhere else.  F-039.
            if ap.isMonomial(pr.relationsFrom(reduced)):
                assert after == quiet(qm.cartanMatrix, reduced, False)
            walk(reduced, depth - 1)

    for length in (5, 6):
        for relSet in qm.generateAllPossibleLineRelations(length):
            relLengths = [0] * (length - 2)
            for rel in relSet:
                relLengths[rel[0][0] - 1] = len(rel[0]) - 1
            walk(line_algebra(length, relLengths), 3)
    assert checked > 1000
