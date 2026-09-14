"""The mutation procedure itself, checked against the worked examples of

    D. Fosse, "A combinatorial procedure for tilting mutation",
    arXiv:2112.08129, section "Statement of main theorem".
"""

import pytest

import quivermutation as qm
from helpers import arrow_set, path_algebra, quiet, rel_set


def test_paper_example_acyclic():
    """The first worked example of the paper: mutate the 7-vertex quiver at 3.

    Before:  1 -a-> 2 -b-> 3 -{c,e}-> {4,5} -{d,z}-> 6 -h-> 7
             relations  ba = 0,  dc + ze = 0,  hd = 0,  hz = 0.
    """
    before = path_algebra(
        arrows=[(1, 2), (2, 3), (3, 4), (3, 5), (4, 6), (5, 6), (6, 7)],
        rels=[[(1, 2, 3)], [(3, 4, 6), (3, 5, 6)], [(4, 6, 7)], [(5, 6, 7)]],
    )

    after = quiet(qm.quiverMutationAtVertex, before, 3)

    # Step 1 adds the compositions through 3, step 2 flips the arrows out of 3,
    # step 3 turns the relation out of 3 into an arrow.  The arrows d: 4 -> 6
    # and z: 5 -> 6 survive mutation but are cancelled during reduction.
    assert arrow_set(after) == {
        (1, 2),
        (2, 4), (2, 5),          # step 1: cb and eb
        (4, 3), (5, 3),          # step 2: c* and e*
        (3, 6),                  # step 3: the relation dc + ze becomes an arrow
        (4, 6), (5, 6),          # d and z, cancelled below
        (6, 7),
    }

    reduced = quiet(qm.reducePathAlgebra, after)

    assert arrow_set(reduced) == {
        (1, 2), (2, 4), (2, 5), (4, 3), (5, 3), (3, 6), (6, 7),
    }
    assert rel_set(reduced) == {
        ((2, 4, 3), (2, 5, 3)),  # step 4: c*cb + e*eb = 0
        ((1, 2, 4),),            # step 6: extend ba = 0 along c
        ((1, 2, 5),),            # step 6: extend ba = 0 along e
        ((3, 6, 7),),            # step 7: h * (the new arrow) = 0
    }


def test_mutation_at_a_vertex_killed_by_a_relation_is_forbidden():
    """A_5 with the relation 1->2->3->4 = 0 cannot be mutated at vertex 3.

    Vertex 3 has the single arrow 3 -> 4 out of it, and the nonzero path
    1 -> 2 -> 3 becomes zero when composed with it, so Hom(P_3*[1], Lambda) is
    nonzero and the result would not be a tilting complex.  This is the case
    the paper singles out: "if there is only one arrow alpha starting in i,
    then mutation is not possible in i if there is a minimal zero relation
    whose last arrow is alpha".
    """
    pa = path_algebra([(1, 2), (2, 3), (3, 4), (4, 5)], [[(1, 2, 3, 4)]])

    assert not quiet(qm.mutationIsPossibleAtVertex, pa, 3)
    # 5 is a sink, so there is nothing to approximate P_5 by.  1, 2 and 4 are
    # each fine: the arrow out of them kills no nonzero path.
    assert not quiet(qm.mutationIsPossibleAtVertex, pa, 5)
    for vertex in (1, 2, 4):
        assert quiet(qm.mutationIsPossibleAtVertex, pa, vertex)


def test_every_vertex_but_the_sink_is_mutable_without_relations():
    pa = path_algebra([(1, 2), (2, 3), (3, 4)])
    assert [quiet(qm.mutationIsPossibleAtVertex, pa, v) for v in (1, 2, 3, 4)] == [
        True, True, True, False,
    ]


@pytest.mark.parametrize(
    "arrows, rels, vertex",
    [
        ([(1, 2), (2, 3), (3, 4), (4, 5)], [], 2),
        ([(1, 2), (2, 3), (3, 4), (4, 5)], [[(1, 2, 3, 4)]], 4),
        ([(1, 2), (2, 3), (3, 4), (4, 5)], [[(2, 3, 4, 5)]], 1),
        ([(1, 2), (2, 3), (3, 4), (4, 5)], [[(1, 2, 3)], [(2, 3, 4)]], 4),
    ],
)
def test_left_mutation_undoes_right_mutation(arrows, rels, vertex):
    """mu_i^L(mu_i^R(Lambda)) = Lambda, where mutation at i is legal.

    quiverMutationAtVertices reads a negative vertex as a left mutation, and
    left mutation is the same procedure with every arrow reversed.
    """
    start = path_algebra(arrows, rels)
    assert quiet(qm.mutationIsPossibleAtVertex, start, vertex)

    there = quiet(qm.quiverMutationAtVertices, path_algebra(arrows, rels), [vertex])
    back = quiet(qm.quiverMutationAtVertices, there, [-vertex])

    assert arrow_set(back) == arrow_set(start)
    assert rel_set(back) == rel_set(start)


@pytest.mark.xfail(
    reason="Step 3's cyclic case is not implemented: for a minimal relation "
           "r: i --> i the procedure calls for one arrow (alpha r-bar): i* -> "
           "t(alpha) per arrow alpha out of i, but the code adds a single "
           "arrow from r's source to r's target, which is a loop i* -> i*. "
           "The LNA search never reaches this case because it stops "
           "descending as soon as a cycle appears.",
    strict=True,
)
def test_paper_example_with_a_cycle():
    """The second worked example of the paper: 1 <=> 2 with (ba): 1 --> 1."""
    before = path_algebra(arrows=[(1, 2), (2, 1)], rels=[[(1, 2, 1)]])

    after = quiet(qm.quiverMutationAtVertex, before, 1)

    assert arrow_set(after) == {
        (2, 2),   # step 1: ab
        (2, 1),   # step 2: a*
        (1, 2),   # step 3: a * (ba)-bar, one per arrow out of 1
    }
