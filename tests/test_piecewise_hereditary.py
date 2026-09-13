"""Certificates that an LNA is not piecewise hereditary, and canonical types.

A piecewise hereditary algebra's derived category is that of a hereditary abelian
category, and by Happel's classification such a category is either the module
category of a hereditary algebra or derived equivalent to a canonical algebra.
So an LNA that is not piecewise hereditary lies in no quipu class at all, and one
that is piecewise hereditary but not of tree type should be of canonical type
with the Coxeter polynomial naming the weights.

The two criteria are from arXiv:2310.08346, whose own claims are what these check
against: every LNA of length at most 8 is piecewise hereditary, and exactly one of
length 9 is not.
"""

import pytest
import sympy

import nakayama as nk
import piecewiseHereditary as pwh


# The paper's (**): A_9 with relations 1 -> 4, 3 -> 6, 4 -> 7, 6 -> 9.
PAPER_A9 = "3033030"
# The smallest algebra the big-overlap criterion applies to: A_13, 1 -> 10, 4 -> 13.
PAPER_A13 = [9, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0]


def test_the_papers_A9_example_is_the_expected_quiver():
    algebra = nk.LinearNakayamaAlgebra(9, PAPER_A9)
    assert [(start, start + arrows) for start, arrows in algebra.relations()] == [
        (1, 4), (3, 6), (4, 7), (6, 9),
    ]


def test_the_papers_A9_example_is_certified():
    name, (alpha, beta) = pwh.certificate(9, [int(c) for c in PAPER_A9])
    assert "Proposition A9" in name
    assert alpha == (3, 3) and beta == (4, 3)     # the overlapping pair
    assert pwh.overlapInArrows(alpha, beta) == 2


def test_the_papers_A13_example_is_certified():
    name, (alpha, beta) = pwh.certificate(13, PAPER_A13)
    assert "Proposition A13" in name
    assert pwh.overlapInArrows(alpha, beta) >= 6


@pytest.mark.parametrize("length", range(4, 9))
def test_nothing_of_length_at_most_8_is_certified(length):
    """The paper says every LNA of length <= 8 is piecewise hereditary."""
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        assert not pwh.isNotPiecewiseHereditary(length, algebra.relLengths), algebra


def test_exactly_one_LNA_of_length_9_is_certified():
    """And the paper says exactly one of length 9 is not -- this one."""
    certified = [
        algebra.className()
        for algebra in nk.LinearNakayamaAlgebra.allOfLength(9)
        if pwh.isNotPiecewiseHereditary(9, algebra.relLengths)
    ]
    assert certified == [PAPER_A9]


def test_the_flanked_pair_criterion_needs_nine_vertices():
    """Its statement assumes the quiver has at least 9 vertices."""
    assert pwh.failsFlankedPairCriterion(8, [3, 0, 3, 3, 0, 3]) is None


def test_overlap_in_arrows():
    assert pwh.overlapInArrows((1, 3), (3, 3)) == 1     # arrows 1,2,3 and 3,4,5
    assert pwh.overlapInArrows((3, 3), (4, 3)) == 2
    assert pwh.overlapInArrows((1, 9), (4, 9)) == 6
    assert pwh.overlapInArrows((1, 2), (5, 2)) == 0


# -- canonical algebras ---------------------------------------------------

def test_the_canonical_coxeter_polynomial_and_vertex_count():
    x = sympy.Symbol("x")
    assert pwh.canonicalCoxeterPolynomial((2,), x) == sympy.expand((x - 1) ** 2 * (1 + x))
    for weights in [(2, 4, 4), (3, 3, 3), (2, 2, 2, 2), (2, 3, 6)]:
        polynomial = pwh.canonicalCoxeterPolynomial(weights, x)
        assert sympy.Poly(polynomial, x).degree() == pwh.canonicalVertexCount(weights)


def test_the_four_tubular_types():
    for weights in pwh.TUBULAR_TYPES:
        assert pwh.isTubular(weights)
        # A tubular Coxeter polynomial has (x-1)^2 as a factor.
        x = sympy.Symbol("x")
        polynomial = pwh.canonicalCoxeterPolynomial(weights, x)
        assert sympy.rem(polynomial, (x - 1) ** 2, x) == 0
    assert not pwh.isTubular((2, 3, 7))


def test_the_non_quipu_class_of_length_9_is_tubular_of_type_2_4_4():
    """A_9 class 3345000 is piecewise hereditary but of canonical type.

    It is not certified non-piecewise-hereditary -- and could not be, since the
    paper says only 3033030 is -- and it is not a quipu class, since its Coxeter
    polynomial is not that of any quipu of order 9.  Its polynomial is exactly
    that of the canonical algebra of tubular weight type (2,4,4).
    """
    algebra = nk.LinearNakayamaAlgebra(9, "3345000")
    assert algebra.quipu() is None
    assert not pwh.isNotPiecewiseHereditary(9, algebra.relLengths)
    weights = pwh.canonicalWeightType(9, algebra.coxeterPolynomial())
    assert weights == (2, 4, 4)
    assert pwh.isTubular(weights)


def test_the_extended_dynkin_class_is_both_a_quipu_and_canonical():
    """D~_5 is tame hereditary, so it is a quipu *and* of canonical type (2,2,3).

    The two identifications agree rather than compete, which is the point: a
    canonical weight type is only used to name a class the quipu form cannot.
    """
    algebra = nk.LinearNakayamaAlgebra(6, "3030")
    assert algebra.quipuName() == "P^(1,1)_(1,0,1)"
    assert pwh.canonicalWeightType(6, algebra.coxeterPolynomial()) == (2, 2, 3)


def test_weight_types_of_an_order_all_have_that_many_vertices():
    for order in range(4, 11):
        types = pwh.weightTypesOfOrder(order)
        assert types
        for weights in types:
            assert pwh.canonicalVertexCount(weights) == order
            assert all(weight >= 2 for weight in weights)
            assert list(weights) == sorted(weights)


def test_the_stored_polynomial_string_parses():
    """The table stores them printed with the variable called 'lambda', which
    Python's tokenizer will not accept as an identifier."""
    text = "lambda**9 + lambda**8 - 2*lambda**5 - 2*lambda**4 + lambda + 1"
    parsed = pwh.parseCoxeterPolynomial(text)
    assert pwh.canonicalWeightType(9, text) == (2, 4, 4)
    assert parsed.free_symbols == {sympy.Symbol("x")}
