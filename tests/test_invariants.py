"""Derived invariants must survive mutation.

The Coxeter polynomial of an algebra is invariant under derived equivalence, and
the mutation procedure produces a derived equivalent algebra at every step, so
it must be constant along any legal mutation sequence.  That makes it the
sharpest cheap check on the mutation engine: any step that mangles the arrows or
the relations will almost always move the polynomial.
"""

import pytest

import quivermutation as qm
from helpers import coxeter_poly, line_algebra, quiet

# (length, per-vertex relation lengths) covering: no relations, a single
# relation of each length, two separate relations, two overlapping relations,
# and a chain of length-2 relations.
LNAS = [
    (5, "000"),
    (5, "300"),
    (5, "030"),
    (5, "230"),
    (5, "400"),
    (6, "0000"),
    (6, "3030"),
    (6, "2300"),
    (6, "0400"),
    (7, "00000"),
    (7, "30300"),
    (7, "22230"),
    (7, "03030"),
    (7, "40030"),
]


@pytest.mark.parametrize("length, rels", LNAS)
def test_coxeter_polynomial_is_constant_along_every_mutation_path(length, rels):
    """Walk every legal mutation of depth <= 3 and check the polynomial holds."""
    start = line_algebra(length, rels)
    expected = coxeter_poly(start)

    def walk(path_alg, depth, history):
        if depth == 0:
            return
        all_rels = quiet(qm.allRelsInPathAlgebra, path_alg)
        for vertex in path_alg.vertices():
            if not quiet(qm.mutationIsPossibleAtVertex, path_alg, vertex):
                continue
            mutated = quiet(qm.quiverMutationAtVertex, path_alg, vertex)
            if any(quiet(qm.isIllegalRelation, mutated, rel) for rel in mutated.rels):
                continue
            mutated = quiet(qm.reducePathAlgebra, mutated)
            got = coxeter_poly(mutated)
            assert got == expected, (
                f"A_{length}_{rels} mutated at {history + [vertex]}: "
                f"got {got}, expected {expected}"
            )
            walk(mutated, depth - 1, history + [vertex])

    walk(start, 3, [])


@pytest.mark.parametrize("length, rels", LNAS)
def test_the_relation_dual_is_derived_equivalent(length, rels):
    """Reversing every arrow of an LNA keeps it in the same class.

    This is the third of the class-preserving operations of arXiv:2305.06642,
    and here it says the opposite algebra of an LNA has the same Coxeter
    polynomial.
    """
    start = line_algebra(length, rels)
    opposite = quiet(qm.dualPathAlgebra, start)
    # Relabel the reversed line back to 1 -> 2 -> ... -> n so the Cartan matrix
    # is built on the standard numbering.
    relabelled, _ = quiet(qm.relabelLineAlgebra, opposite, {})
    assert coxeter_poly(relabelled) == coxeter_poly(start)
