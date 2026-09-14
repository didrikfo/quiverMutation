"""The mutation procedure on linear combinations, against the model it replaced.

`procedure` is the procedure of arXiv:2112.08129 with relations as
`relationAlgebra` combinations: the coefficients come out of the steps rather
than being guessed back from a set of paths.  It is what
`mutation.quiverMutationAtVertex` and `reduction.reducePathAlgebra` now run.

Three things are worth pinning, and they are different things:

* that it gives what the set-of-paths implementation gave, on everything that
  implementation was trusted for -- which is what makes replacing it safe;
* that the coefficients it produces are the ones the paper's steps call for, on
  a worked example small enough to read;
* that its admissibility condition is the paper's, which is *not* what the
  search uses -- see the module docstring of `mutation`.
"""

import collections
import contextlib
import io

import pytest
import sympy

import quivermutation as qm
from quivermutation import nakayama as nk
from quivermutation import procedure as pr
from quivermutation import relationAlgebra as ra
from helpers import path_algebra, quiet


def shape(pathAlg):
    """A hashable form of a path algebra: its arrows and its relation sets."""
    arrows = tuple(sorted((a[0], a[1]) for a in pathAlg.quiver.edges()))
    rels = tuple(sorted(tuple(sorted(tuple(p) for p in rel)) for rel in pathAlg.rels))
    return arrows, rels


# -- the steps, read off a small example ----------------------------------

def test_the_steps_produce_the_coefficients_the_paper_asks_for():
    """A_4 with 1 -> 2 -> 3 -> 4 = 0, mutated at 3.

    Vertex 3 has one arrow out (to 4) and one in (from 2), and one relation
    starts at 3 -- none, here -- so the steps that fire are 1, 2, 4 and 6:

    * step 1 composes 2 -> 3 and 3 -> 4 into 2 -> 4;
    * step 2 flips 3 -> 4 to 4 -> 3;
    * step 4 turns the arrow 2 -> 3 into the relation (2, 4, 3) = 0, a sum over
      the single arrow out of 3;
    * step 6 extends the relation 1 -> 4 through the new composite, giving
      (1, 2, 4) = 0.

    The coefficients are all 1 here because each sum has one term; the point of
    the test is that the *paths* are right and no relation arrives with a
    coefficient the model cannot express.
    """
    algebra = nk.LinearNakayamaAlgebra(4, "30")
    quiver, relations = pr.mutateAtVertex(algebra.quiver, pr.relationsFrom(algebra), 3)

    assert sorted(quiver.edges()) == [(1, 2), (2, 4), (4, 3)]
    asWritten = {tuple(sorted(r.items())) for r in relations}
    assert asWritten == {(((1, 2, 4), 1),), (((2, 4, 3), 1),)}


def test_step_five_is_a_difference_not_a_sum():
    """`rbar alpha* = r / alpha` is a relation with a minus sign in it.

    A_5 with 1 -> 2 -> 3 = 0 and 3 -> 4 -> 5 = 0, mutated at 3: the relation
    out of 3 becomes the arrow 3 -> 5, the arrow 3 -> 4 flips to 4 -> 3, and
    step 5 says (4, 3, 5) equals the part of the relation beginning with
    3 -> 4, which is the path (4, 5).  So the relation is (4, 3, 5) - (4, 5),
    and reading it as a sum would be a different relation.
    """
    algebra = nk.LinearNakayamaAlgebra(5, "202")   # relations 1 -> 3 and 3 -> 5
    quiver, relations = pr.mutateAtVertex(algebra.quiver, pr.relationsFrom(algebra), 3)
    stepFive = [r for r in relations if sorted(r) == [(4, 3, 5), (4, 5)]]
    assert stepFive, [sorted(r.items()) for r in relations]
    assert stepFive[0][(4, 3, 5)] == -stepFive[0][(4, 5)]


# -- agreement with the implementation it replaced ------------------------

@pytest.mark.parametrize("length", [4, 5, 6])
def test_it_agrees_with_the_set_of_paths_procedure_on_every_lna(length):
    """One mutation at every admissible vertex of every LNA of the length.

    Gated on `mutationIsPossibleAtVertex`, the stricter condition, so both are
    compared on exactly the mutations the published classification used.
    """
    checked = 0
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        relations = pr.relationsFrom(algebra)
        for vertex in range(1, length + 1):
            if not quiet(qm.mutationIsPossibleAtVertex, algebra, vertex):
                continue
            quiver, mutated = pr.reduce(
                *pr.mutateAtVertex(algebra.quiver, relations, vertex))
            assert shape(pr.toPathAlgebra(quiver, mutated)) == shape(
                quiet(qm.reducePathAlgebra, quiet(qm.quiverMutationAtVertex, algebra, vertex))
            ), (algebra, vertex)
            checked += 1
    assert checked > 0


@pytest.mark.slow
@pytest.mark.parametrize("length", [7, 8])
def test_it_agrees_with_the_set_of_paths_procedure_further_out(length):
    test_it_agrees_with_the_set_of_paths_procedure_on_every_lna(length)


def test_the_coefficients_survive_a_chain_of_mutations():
    """A walk keeps its coefficients rather than guessing them again each step.

    `toPathAlgebra` records them on the algebra and `relationsFrom` reads them
    back, so three mutations in a row are exact throughout.  The cache is
    checked against the relation sets, so editing `rels` behind its back falls
    back to the guess rather than going stale.
    """
    algebra = nk.LinearNakayamaAlgebra(6, "3030")
    walked = quiet(qm.quiverMutationAtVertices, algebra, [4, 1, 2])
    assert walked.relCombinations is not None
    assert [ra.toPathSet(c) for c in walked.relCombinations] == walked.rels

    walked.rels = [[[1, 2, 3]]]
    assert pr.relationsFrom(walked) == [ra.combination([(1, 2, 3)])]


def test_the_opposite_algebra_keeps_the_coefficients():
    """Left mutation goes through the opposite algebra, which must not lose signs."""
    relations = [ra.combination([((1, 2, 4), 1), ((1, 3, 4), -1)])]
    quiver, opposed = pr._opposite(
        path_algebra([(1, 2), (1, 3), (2, 4), (3, 4)]).quiver, relations)
    assert sorted(quiver.edges()) == [(2, 1), (3, 1), (4, 2), (4, 3)]
    assert opposed == [ra.combination([((4, 2, 1), 1), ((4, 3, 1), -1)])]


# -- the reduction --------------------------------------------------------

def test_the_note_after_step_seven_substitutes_rather_than_deletes():
    """A relation with a path of length one removes that arrow from the quiver.

    `alpha + p = 0` says `alpha = -p`, so the arrow goes and every path through
    it is rewritten.  Here 1 -> 4 is a single arrow in a relation with the path
    1 -> 2 -> 3 -> 4, so the arrow goes and the relation on 1 -> 4 -> 5 becomes
    one on 1 -> 2 -> 3 -> 4 -> 5.
    """
    quiver = path_algebra([(1, 2), (2, 3), (3, 4), (1, 4), (4, 5)]).quiver
    relations = [
        ra.combination([((1, 4), 1), ((1, 2, 3, 4), 1)]),
        ra.combination([((1, 4, 5), 1)]),
    ]
    reducedQuiver, reduced = pr.reduce(quiver, relations)
    assert (1, 4) not in list(reducedQuiver.edges())
    assert [ra.toPathSet(r) for r in reduced] == [[[1, 2, 3, 4, 5]]]


def test_minimality_is_decided_over_the_ideal_not_by_containment():
    """A relation the others imply is dropped even with no path in common.

    `1 -> 2 -> 4 = 1 -> 3 -> 4` and `1 -> 2 -> 4 = 0` together imply
    `1 -> 3 -> 4 = 0`, which shares no path with either.  A syntactic test for
    a subrelation cannot see that; reducing against a basis of the ideal can.
    """
    quiver = path_algebra([(1, 2), (1, 3), (2, 4), (3, 4)]).quiver
    commutativity = ra.combination([((1, 2, 4), 1), ((1, 3, 4), -1)])
    oneIsZero = ra.combination([((1, 2, 4), 1)])
    theOther = ra.combination([((1, 3, 4), 1)])
    _, reduced = pr.reduce(quiver, [commutativity, oneIsZero, theOther])
    assert len(reduced) == 2, [ra.toPathSet(r) for r in reduced]
    assert ra.isInIdeal(quiver, reduced, theOther)


# -- admissibility --------------------------------------------------------

@pytest.mark.parametrize("length", [4, 5, 6])
def test_on_an_lna_the_two_admissibility_conditions_agree(length):
    """Every vertex of a line has one arrow out, which is where they coincide."""
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        relations = pr.relationsFrom(algebra)
        for vertex in range(1, length + 1):
            assert pr.isMutable(algebra.quiver, relations, vertex) is bool(
                quiet(qm.mutationIsPossibleAtVertex, algebra, vertex)), (algebra, vertex)


def test_where_they_differ_the_paper_allows_more_and_the_class_is_kept():
    """The stricter condition refuses mutations the paper permits.

    After one mutation a vertex can have two arrows out, and there
    `mutationIsPossibleAtVertex` rejects as soon as *one* of them kills a
    nonzero path where the paper only asks that *one keeps it*.  Every mutation
    that difference refuses is legitimate, so the Coxeter polynomial must not
    move across it -- which is the check R-005 exists to insist on.

    Research F-015 is the exhaustive version of this.
    """
    tally = collections.Counter()
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(5):
        relations = pr.relationsFrom(algebra)
        for first in range(1, 6):
            if not quiet(qm.mutationIsPossibleAtVertex, algebra, first):
                continue
            quiver, mutated = pr.reduce(*pr.mutateAtVertex(algebra.quiver, relations, first))
            here = pr.toPathAlgebra(quiver, mutated)
            base = sympy.expand(quiet(qm.coxeterPoly, here).as_expr())
            for second in sorted(quiver.nodes):
                strict = bool(quiet(qm.mutationIsPossibleAtVertex, here, second))
                exact = pr.isMutable(quiver, mutated, second)
                if strict == exact:
                    continue
                assert (strict, exact) == (False, True), (algebra, first, second)
                assert len(set(quiver.successors(second))) > 1, (algebra, first, second)
                onward = pr.toPathAlgebra(*pr.reduce(
                    *pr.mutateAtVertex(quiver, mutated, second)))
                assert sympy.expand(quiet(qm.coxeterPoly, onward).as_expr()) == base
                tally['refused by the strict condition'] += 1
    assert tally['refused by the strict condition'] > 0


def test_step_seven_finds_the_relation_the_old_implementation_missed():
    """The two cases at n = 7 where the set-of-paths procedure lost a relation.

    Step 7 is an iff over coefficients that are themselves paths, so the
    relations out of `i*` are the whole kernel -- not only the combinations of
    the arrows `r-bar`.  Reading it the narrow way loses relations, and these
    are the smallest cases found where it does: three mutations out of an LNA of
    length 7.  Research R-007 and F-015.

    `2;5;6;7 = 0` is the relation in question, and it is *not* in the ideal the
    old answer generated -- so the two were different algebras, not two
    presentations of one.  Nothing that was being checked caught it: both
    answers have the same Coxeter polynomial, both come back under a left
    mutation at the same vertex, and both reach the quipu the theorem names.
    """
    for name in ("40030", "44030"):
        algebra = nk.LinearNakayamaAlgebra(7, name)
        quiver, relations = algebra.quiver, pr.relationsFrom(algebra)
        for vertex in (1, 4, 2):
            quiver, relations = pr.reduce(*pr.mutateAtVertex(quiver, relations, vertex))
        reached = pr.toPathAlgebra(quiver, relations)
        assert [2, 5, 6, 7] in reached.rels[0] or any(
            [2, 5, 6, 7] in rel for rel in reached.rels), (name, reached.rels)
        # and it is a relation the paper's step 7 requires, not a consequence of
        # the others: dropping it leaves an ideal it is not in
        others = [r for r in relations if sorted(r) != [(2, 5, 6, 7)]]
        assert not ra.isInIdeal(quiver, others, ra.combination([(2, 5, 6, 7)])), name
