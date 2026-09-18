"""Quivers with parallel arrows: the model, the procedure, and the invariant.

The mutation procedure of arXiv:2112.08129 produces parallel arrows -- step 1
adds a composite `alpha beta: h -> j` whether or not `h -> j` is an arrow
already, and step 3 adds one arrow `i* -> k` per relation `i ~~> k`.  Until the
relations named their arrows, the repo could not state such a quiver: a path was
a sequence of vertices, the Cartan matrix counted one path where there were two,
and the gate refused mutation at every vertex of a quiver that had a parallel
pair anywhere.  Research F-038 recorded the consequence and called it harmless;
these tests are what says it is now correct instead of avoided.

What each part pins:

* `arrowPaths` counts two parallel arrows as two paths, and lifting a vertex
  sequence is refused exactly where it is ambiguous;
* the Coxeter polynomial of the Kronecker algebra, which no other test in the
  repo can even build;
* the mutation that produced the smallest wrong key of F-038 keeps its key now,
  and the mutations out of the quiver it lands on keep it too;
* step 3 gives two arrows when two relations share both ends, and step 7 reads
  each of them back as the relation it came from rather than giving up;
* a relation carried past the mutated vertex tells the new composite arrow from
  an arrow that was there already.
"""

import networkx as nx
import pytest
import sympy

import quivermutation as qm
from quivermutation import arrowPaths as ap
from quivermutation import invariants as inv
from quivermutation import mutation
from quivermutation import nakayama as nk
from quivermutation import pathAlgebra as pa
from quivermutation import procedure as pr
from quivermutation import relationAlgebra as ra
from quivermutation import reduction
from quivermutation import search

from helpers import quiet


def quiver(arrows):
    """A MultiDiGraph from (tail, head) pairs, keys assigned in order."""
    q = nx.MultiDiGraph()
    for tail, head in arrows:
        q.add_edge(tail, head)
    return q


def walk(algebra, vertices):
    """Mutate at each vertex in turn, reducing after each."""
    for vertex in vertices:
        algebra = quiet(reduction.reducePathAlgebra,
                        quiet(mutation.quiverMutationAtVertex, algebra, vertex))
    return algebra


# -- the model ------------------------------------------------------------

def test_two_parallel_arrows_are_two_paths():
    q = quiver([(1, 2), (1, 2)])
    assert ap.arrowsOf(q) == [(1, 2, 0), (1, 2, 1)]
    assert ap.hasParallelArrows(q)
    assert ap.allPathsBetween(q, 1, 2) == [((1, 2, 0),), ((1, 2, 1),)]
    assert ap.homDimension(q, [], 1, 2) == 2
    assert ap.homDimensionByClosure(q, [], 1, 2) == 2


def test_a_vertex_sequence_cannot_be_lifted_along_a_parallel_pair():
    q = quiver([(1, 2), (1, 2), (2, 3)])
    with pytest.raises(ValueError, match = "parallel"):
        ap.liftPath(q, [1, 2, 3])
    assert ap.liftPath(quiver([(1, 2), (2, 3)]), [1, 2, 3]) == ((1, 2, 0), (2, 3, 0))
    with pytest.raises(ValueError, match = "not an arrow"):
        ap.liftPath(quiver([(1, 2)]), [1, 3])


def test_lifting_and_projecting_are_inverse_where_nothing_is_parallel():
    algebra = nk.LinearNakayamaAlgebra(6, "3030")
    relations = pr.relationsFrom(algebra)
    assert ap.projectToPathSets(relations) == algebra.rels
    assert ap.describesRels(relations, algebra.rels)
    assert relations == ap.lift(algebra.quiver, algebra.rels)


def test_a_combination_of_two_arrow_paths_is_not_read_as_a_pair():
    """The obvious test for `(path, coefficient)` is wrong here.

    An arrow path of two arrows is itself a two-element tuple, so a pair has to
    be told by its *second* element being a number rather than an arrow.
    """
    twoArrows = ((1, 2, 0), (2, 3, 0))
    assert ap.combination([twoArrows]) == {twoArrows: 1}
    assert ap.combination([(twoArrows, -3)]) == {twoArrows: -3}


def test_dividing_and_composing_are_by_arrow_not_by_endpoint():
    """`r / alpha` is division by an arrow, which is the point of naming them.

    On two arrows `1 -> 2`, the relation `alpha_0 beta - alpha_1 beta` divides to
    `beta` on one and `-beta` on the other; dividing by "the arrow to 2" could
    only give one answer for both.
    """
    q = quiver([(1, 2), (1, 2), (2, 3)])
    first, second, onward = (1, 2, 0), (1, 2, 1), (2, 3, 0)
    relation = ap.combination([((first, onward), 1), ((second, onward), -1)])
    assert ap.source(relation) == 1 and ap.target(relation) == 3
    assert ap.leftDivide(relation, first) == {(onward,): 1}
    assert ap.leftDivide(relation, second) == {(onward,): -1}
    assert ap.rightDivide(relation, onward) == {(first,): 1, (second,): -1}
    assert ap.preCompose(ap.combination([(onward,)]), (first,)) == {(first, onward): 1}
    assert ap.postCompose(ap.combination([(first,)]), (onward,)) == {(first, onward): 1}
    # the trivial path is the identity of the composition
    assert ap.preCompose(relation, ()) == relation
    assert ap.postCompose(relation, ()) == relation


# -- the invariant --------------------------------------------------------

def test_the_kronecker_algebra_has_the_coxeter_polynomial_it_should():
    """Two arrows 1 -> 2 and no relations: the path algebra of the Kronecker quiver.

    Its Cartan matrix is [[1, 0], [2, 1]] -- two paths from 1 to 2 -- and its
    Coxeter polynomial is `(lambda - 1)^2`.  In the vertex model the matrix read
    [[1, 0], [1, 1]] and the polynomial came out as that of `A_2`, so the
    smallest quiver with a parallel pair was already mis-measured.
    """
    algebra = pa.PathAlgebra()
    algebra.quiver = quiver([(1, 2), (1, 2)])
    assert inv.integerCartanMatrix(algebra) == [[1, 0], [2, 1]]
    lam = sympy.Symbol("lambda")
    assert sympy.expand(inv.coxeterPoly(algebra).as_expr()) == sympy.expand((lam - 1) ** 2)
    assert inv.coxeterKey(algebra) == (1, -2, 1)


def test_the_cheap_and_exact_counts_agree_on_a_quiver_with_parallel_arrows():
    algebra = walk(nk.LinearNakayamaAlgebra(6, "3030"), [1, 3, 4, 1, 4])
    assert algebra.hasParallelArrows()
    relations = pr.relationsFrom(algebra)
    assert (ap.cartanMatrix(algebra.quiver, relations, exact = True)
            == ap.cartanMatrix(algebra.quiver, relations, exact = False))


# -- the procedure --------------------------------------------------------

def test_the_smallest_wrong_key_of_f038_was_the_invariant_not_the_mutation():
    """`3030` at n = 6, mutated at [1, 3, 4, 1, 4], is the first parallel pair.

    E-033 found this by walking the search tree with the Coxeter key in hand: the
    fifth step produces two arrows `1 -> 6` and the key moved from
    `(1, 1, -1, -2, -1, 1, 1)` to `(1, 1, 0, -1, 0, 1, 1)`.  The mutation was
    right and the *measurement* was wrong -- counting the parallel pair once.
    Every step is admissible and the key holds the whole way now.
    """
    algebra = nk.LinearNakayamaAlgebra(6, "3030")
    base = inv.coxeterKey(algebra)
    assert base == (1, 1, -1, -2, -1, 1, 1)
    here = algebra
    for vertex in [1, 3, 4, 1, 4]:
        assert quiet(mutation.mutationIsPossibleAtVertex, here, vertex)
        here = walk(here, [vertex])
        assert inv.coxeterKey(here) == base, sorted(here.quiver.edges(keys = True))
    assert here.hasParallelArrows()
    assert sorted(here.quiver.edges(keys = True)) == [
        (1, 6, 0), (1, 6, 1), (2, 5, 0), (3, 1, 0), (5, 1, 0), (6, 4, 0)]


def test_a_quiver_with_parallel_arrows_is_mutable_and_stays_in_its_class():
    """The gate used to refuse every vertex of it, so the branch died there.

    The refusal was honest about being a restriction of the model rather than of
    the procedure, and naming the arrows removes it.  Every vertex the paper's
    criterion now admits keeps the Coxeter key, including the one that leads back
    out to a quiver with no parallel arrows at all -- which is a region of the
    mutation graph the search could not reach before.
    """
    algebra = walk(nk.LinearNakayamaAlgebra(6, "3030"), [1, 3, 4, 1, 4])
    base = inv.coxeterKey(nk.LinearNakayamaAlgebra(6, "3030"))
    assert algebra.hasParallelArrows()
    admissible = [v for v in algebra.vertices()
                  if quiet(mutation.mutationIsPossibleAtVertex, algebra, v)]
    assert admissible == [1, 2, 3]
    assert [v for v in algebra.vertices()
            if quiet(mutation.mutationIsPossibleAtVertex, algebra, v,
                     allowParallelArrows = False)] == []
    for vertex in admissible:
        assert inv.coxeterKey(walk(algebra, [vertex])) == base, vertex
    assert not walk(algebra, [3]).hasParallelArrows()


def test_a_relation_past_the_mutated_vertex_keeps_the_arrow_it_ran_along():
    """Splicing the vertex out of a path gives the *composite* arrow, not any arrow.

    At the fifth step above, the relation `5 -> 1 -> 4 -> 6 = 5 -> 1 -> 6` runs
    through the mutated vertex 4 on one side only.  Afterwards its two paths use
    the two different arrows `1 -> 6`: the composite of `1 -> 4` and `4 -> 6`,
    and the arrow `1 -> 6` that was there before.  As vertex sequences both read
    `5, 1, 6` and the relation collapsed to nothing.
    """
    algebra = walk(nk.LinearNakayamaAlgebra(6, "3030"), [1, 3, 4, 1, 4])
    relations = pr.relationsFrom(algebra)
    parallel = [r for r in relations if len(r) == 2]
    assert len(parallel) == 1, relations
    paths = sorted(parallel[0])
    assert [ap.projectPath(p) for p in paths] == [[5, 1, 6], [5, 1, 6]]
    assert paths[0] != paths[1]
    assert {p[-1] for p in paths} == {(1, 6, 0), (1, 6, 1)}


def test_step_three_gives_one_arrow_per_relation_and_step_seven_reads_it_back():
    """Two relations `1 ~~> 4` become two arrows `1* -> 4`.

    A vertex sequence cannot tell them apart, so `procedure` used to skip step 7
    for that target outright -- the old `_relationOfFirstArrow` returned None as
    soon as a target had two relation arrows.  Named, each candidate path's first
    arrow is the relation it came from, and step 5 divides by the arrow it means.
    """
    algebra = pa.PathAlgebra()
    algebra.add_arrows_from([[1, 2], [1, 3], [2, 4], [3, 4]])
    algebra.add_rels_from([[[1, 2, 4]], [[1, 3, 4]]])
    relations = pr.relationsFrom(algebra)
    assert pr.isMutable(algebra.quiver, relations, 1)

    mutatedQuiver, mutated = pr.mutateAtVertex(algebra.quiver, relations, 1)
    assert sorted(mutatedQuiver.out_edges(1, keys = True)) == [(1, 4, 0), (1, 4, 1)]

    reducedQuiver, reduced = pr.reduce(mutatedQuiver, mutated)
    reachedAlgebra = pr.toPathAlgebra(reducedQuiver, reduced)
    assert reachedAlgebra.hasParallelArrows()
    assert inv.coxeterKey(reachedAlgebra) == inv.coxeterKey(algebra)
    # each of the two relations out of 1 gave a relation on its own arrow
    assert sorted(sorted(r)[0][-1] for r in reduced) == [(1, 4, 0), (1, 4, 1)]


def test_the_reduction_removes_the_arrow_the_relation_names():
    """The paper's Note deletes the arrow a relation writes down, not a namesake.

    `1 -> 2` twice, with the *second* of them equal to the path `1 -> 4 -> 2`.
    That relation has one path of length one, so the Note removes that arrow and
    substitutes it out -- and the arrow it removes has to be the one the relation
    names.  `remove_edge` on a pair of endpoints drops whichever of a parallel
    pair networkx reaches first.
    """
    q = quiver([(1, 2), (1, 2), (1, 4), (4, 2)])
    doomed = (1, 2, 1)
    relations = [ap.combination([((doomed,), 1),
                                 (((1, 4, 0), (4, 2, 0)), -1)])]
    reducedQuiver, reduced = pr.reduce(q, relations)
    assert sorted(reducedQuiver.edges(keys = True)) == [(1, 2, 0), (1, 4, 0), (4, 2, 0)]
    assert reduced == []


# -- what it means for the walk ------------------------------------------

def test_the_dual_of_a_parallel_pair_is_a_parallel_pair():
    """Left mutation goes through the opposite algebra, so the dual has to carry
    the arrow names or a parallel-arrow algebra cannot be left-mutated at all."""
    algebra = walk(nk.LinearNakayamaAlgebra(6, "3030"), [1, 3, 4, 1, 4])
    dual = pa.dualPathAlgebra(algebra)
    assert dual.hasParallelArrows()
    assert dual.arrowRels is not None
    assert pr.relationsFrom(dual) == [dict(r) for r in dual.arrowRels]
    assert pa.dualPathAlgebra(dual).rels == algebra.rels


def test_the_projection_of_a_parallel_relation_cannot_be_read_back():
    """Which is why `arrowRels` is the record and `rels` only the table's key.

    The relation between the two parallel paths projects to the vertex sequence
    `5, 1, 6` written twice, and the repo reads a two-path relation as a
    difference -- so reading the projection back gives `p - p`, which is the
    *zero* combination and no relation at all.  Nothing about the vertex model
    can recover it, and `relationsFrom` must not try.
    """
    algebra = walk(nk.LinearNakayamaAlgebra(6, "3030"), [1, 3, 4, 1, 4])
    doubled = [rel for rel in algebra.rels if len(rel) == 2 and rel[0] == rel[1]]
    assert doubled, algebra.rels
    assert ra.fromPathSet(doubled[0]) == {}

    for relation in pr.relationsFrom(algebra):
        assert not ap.isIllegalRelation(algebra.quiver, relation)
    # `relationsFrom` reads the arrow relations, not the projection
    assert algebra.arrowRels is not None
    assert pr.relationsFrom(algebra) == [dict(r) for r in algebra.arrowRels]
    with pytest.raises(ValueError, match = "parallel"):
        ap.lift(algebra.quiver, algebra.rels)


def test_the_cheap_count_has_no_reading_of_a_three_path_relation():
    """Which is why the search's key is exact off a monomial ideal.

    `34400` at `n = 7` by `[1, 3, 4, 2, 2, 1]` reaches a quiver carrying
    `-(1,2,7) + (1,4,7) + (1,6,7) = 0` -- step 4 at a vertex with three arrows
    out.  That relation cuts the span of the paths `1 ~~> 7` from three
    dimensions to two.  The cheap count collects one-path relations as zeros and
    two-path relations as identifications, and has nothing to do with a sum of
    three, so it ignores it and reads 3; the key then moves, which is one of the
    nodes F-038 called "clean" at `n = 7` depth 6.
    """
    algebra = walk(nk.LinearNakayamaAlgebra(7, "34400"), [1, 3, 4, 2, 2, 1])
    assert algebra.rels == [[[1, 2, 7], [1, 4, 7], [1, 6, 7]], [[3, 1, 6]], [[5, 1, 2]]]
    relations = pr.relationsFrom(algebra)
    assert not ap.isMonomial(relations)
    assert ap.homDimensionByClosure(algebra.quiver, relations, 1, 7) == 3
    assert ap.homDimension(algebra.quiver, relations, 1, 7) == 2
    assert inv.coxeterKey(algebra) == inv.coxeterKey(nk.LinearNakayamaAlgebra(7, "34400"))


@pytest.mark.slow
@pytest.mark.parametrize("length, depth", [(6, 4), (7, 3)])
def test_the_cheap_count_is_exact_where_the_key_takes_it(length, depth):
    """The cheap route is taken exactly on a monomial ideal, where it is right.

    `isMonomial` is the condition, and the point of the sweep is that it is the
    *right* condition: on every node of an LNA walk whose relations are all
    single paths, the cheap count and the rank over the ideal agree.
    """
    monomial = [0]

    def visitor(pathAlg, path):
        relations = pr.relationsFrom(pathAlg)
        if not ap.isMonomial(relations):
            return
        monomial[0] += 1
        assert (ap.cartanMatrix(pathAlg.quiver, relations, exact = False)
                == ap.cartanMatrix(pathAlg.quiver, relations, exact = True)), (
                    sorted(pathAlg.quiver.edges(keys = True)), pathAlg.rels, path)

    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        search.mutationSearchDepthFirst(algebra, depth, printOutput = False,
                                        visitor = visitor)
    assert monomial[0] > 100


def test_the_search_walks_through_a_parallel_arrow_quiver():
    """The node used to be terminal, so everything below it was unreachable."""
    algebra = nk.LinearNakayamaAlgebra(6, "3030")
    seen = []

    def visitor(pathAlg, path):
        if pathAlg.hasParallelArrows():
            seen.append(list(path))

    search.mutationSearchDepthFirst(algebra, 6, printOutput = False, visitor = visitor)
    assert seen, "no parallel-arrow quiver was reached at all"
    # and the search went on from at least one of them
    assert any(len(path) > min(len(p) for p in seen) for path in seen)
