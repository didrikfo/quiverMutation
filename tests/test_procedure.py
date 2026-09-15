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
* that its admissibility condition -- now the search's gate -- allows
  everything its stricter predecessor did and more, and that nothing it newly
  allows moves the Coxeter polynomial.
"""

import collections
import contextlib
import io

import pytest
import sympy

import networkx as nx

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

def strictlyMutable(pathAlg, vertex):
    """The admissibility gate the search used before `procedure.isMutable`.

    Kept here, in a test, because it is the thing that was replaced and the
    only way to say what replacing it changed.  It was two tests of one idea,
    in two places, and both were the *strict* reading -- a vertex is refused as
    soon as **one** arrow out of it kills a nonzero path, where the paper's
    theorem refuses it only when a path dies against **every** arrow:

    * `mutation.mutationIsPossibleAtVertex` walked the relations and refused the
      vertex on any minimal zero relation whose last arrow left it and whose
      truncation was not itself written as a relation;
    * `search.mutationSearchDepthFirst` then counted, for every predecessor `v`
      and every arrow `i -> w`, the paths `v -> i` against the paths `v -> w`,
      and refused the vertex if any arrow lost one.

    Both also decided "nonzero" syntactically -- a zero relation written inside
    the path, or a path count up to the commutativity relations -- where the
    replacement decides it over the ideal.
    """
    if not bool(pathAlg.out_arrows(vertex)):
        return False
    if any(arrow[2] > 0 for arrow in pathAlg.arrows()):
        return False
    allRels = quiet(qm.allRelsInPathAlgebra, pathAlg)
    for rel in allRels:
        if len(rel) == 1 and rel[0][-2] == vertex and [rel[0][:-1]] not in allRels:
            return False
    successors = list(pathAlg.quiver.successors(vertex))
    for source in nx.dfs_preorder_nodes(nx.reverse(pathAlg.quiver), vertex):
        intoVertex = quiet(qm.numberOfPathsUpToRels, pathAlg, source, vertex)
        for target in successors:
            if intoVertex > quiet(qm.numberOfPathsUpToRels, pathAlg, source, target):
                return False
    return True


@pytest.mark.parametrize("length", [4, 5, 6])
def test_on_an_lna_the_gate_and_its_predecessor_agree(length):
    """Every vertex of a line has one arrow out, which is where they coincide.

    Which is why the switch does not touch the classification's first step: the
    LNAs it starts from are lines, and the two criteria part company only once a
    mutation has made a vertex branch.
    """
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        relations = pr.relationsFrom(algebra)
        for vertex in range(1, length + 1):
            assert pr.isMutable(algebra.quiver, relations, vertex) is strictlyMutable(
                algebra, vertex), (algebra, vertex)


@pytest.mark.parametrize("length", [5, 6])
def test_the_gate_allows_more_than_its_predecessor_and_keeps_the_class(length):
    """The switch only ever permits more, and never permits something wrong.

    The first half is a statement about the two criteria: after a mutation a
    vertex can have several arrows out, and there the predecessor refused as
    soon as one of them killed a nonzero path.  Every disagreement must be in
    that direction, and at a branching vertex.

    The second half is the part that matters, and is why R-005 exists. The
    paper's criterion rules mutation *out*; it does not rule it *in*, since the
    real condition is on the algebra and is not equivalent to any condition on
    the quiver. So a mutation this gate newly allows could in principle fail to
    be a derived equivalence, and the Coxeter polynomial would move across it.
    It does not, anywhere here.
    """
    tally = collections.Counter()
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        relations = pr.relationsFrom(algebra)
        for first in range(1, length + 1):
            if not pr.isMutable(algebra.quiver, relations, first):
                continue
            quiver, mutated = pr.reduce(*pr.mutateAtVertex(algebra.quiver, relations, first))
            here = pr.toPathAlgebra(quiver, mutated)
            base = sympy.expand(quiet(qm.coxeterPoly, here).as_expr())
            for second in sorted(quiver.nodes):
                wasAllowed = strictlyMutable(here, second)
                isAllowed = pr.isMutable(quiver, mutated, second)
                if wasAllowed == isAllowed:
                    continue
                assert (wasAllowed, isAllowed) == (False, True), (algebra, first, second)
                assert len(set(quiver.successors(second))) > 1, (algebra, first, second)
                onward = pr.toPathAlgebra(*pr.reduce(
                    *pr.mutateAtVertex(quiver, mutated, second)))
                assert sympy.expand(quiet(qm.coxeterPoly, onward).as_expr()) == base, (
                    algebra, first, second)
                tally['newly allowed'] += 1
    assert tally['newly allowed'] > 0


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
