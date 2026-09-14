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


# -- propagating certificates by deleting vertices ------------------------

def test_removing_a_vertex_shortens_the_quiver_by_one():
    result = pwh.removeVertex(9, [int(c) for c in PAPER_A9], 5)
    assert result is not None
    length, relLengths = result
    assert length == 8
    assert len(relLengths) == 6


def test_deleting_an_interior_vertex_merges_two_arrows():
    """A relation spanning both arrows through the deleted vertex loses one;
    one spanning just one of them keeps its length and reaches further.

    A_6 with the single relation 2 -> 5 (three arrows, 2, 3 and 4).  Deleting
    vertex 3 merges arrows 2 and 3, so the relation is left with two arrows.
    Deleting vertex 6 instead touches none of them, so it is unchanged.
    """
    assert pwh.removeVertex(6, [0, 3, 0, 0], 3) == (5, [0, 2, 0])
    assert pwh.removeVertex(6, [0, 3, 0, 0], 6) == (5, [0, 3, 0])


def test_deleting_an_end_vertex_drops_the_relation_that_reaches_it():
    """There is nowhere for such a relation to be extended to."""
    # A_5, relations 1 -> 3 (arrows 1,2) and 2 -> 5 (arrows 2,3,4).
    # Deleting vertex 5 loses arrow 4, so 2 -> 5 keeps two of its three arrows,
    # while 1 -> 3 is untouched.
    assert pwh.removeVertex(5, [2, 3, 0], 5) == (4, [2, 2])
    # Deleting vertex 1 loses arrow 1.  That leaves 1 -> 3 with a single arrow, so
    # it goes; 2 -> 5 never used arrow 1, so it keeps all three and simply shifts
    # down to start at the new vertex 1.
    assert pwh.removeVertex(5, [2, 3, 0], 1) == (4, [3, 0])


def test_removing_a_vertex_always_gives_an_admissible_algebra_or_nothing():
    for length in range(4, 9):
        for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
            for vertex in range(1, length + 1):
                result = pwh.removeVertex(length, algebra.relLengths, vertex)
                if result is None:
                    continue
                shorter, relLengths = result
                assert shorter == length - 1
                # It must be constructible as an LNA, which validates it.
                nk.LinearNakayamaAlgebra(shorter, relLengths)


@pytest.mark.parametrize("length", range(4, 9))
def test_deletion_certifies_nothing_below_length_9(length):
    """The strongest check on the deletion rule.

    Every LNA of length at most 8 is piecewise hereditary, so any certificate
    here would mean the deletion construction is wrong.  A rule that dropped or
    extended the wrong relation would show up at once.
    """
    cache = {}
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        chain = pwh.notPiecewiseHereditaryByDeletion(length, algebra.relLengths, cache)
        assert chain is None, (algebra, chain)


def test_deletion_still_certifies_exactly_one_algebra_of_length_9():
    cache = {}
    certified = [
        algebra.className()
        for algebra in nk.LinearNakayamaAlgebra.allOfLength(9)
        if pwh.notPiecewiseHereditaryByDeletion(9, algebra.relLengths, cache) is not None
    ]
    assert certified == [PAPER_A9]


def test_a_certificate_chain_ends_at_a_direct_criterion():
    cache = {}
    chain = pwh.notPiecewiseHereditaryByDeletion(9, [int(c) for c in PAPER_A9], cache)
    assert chain is not None
    length, relLengths, reason = chain[-1]
    assert "Proposition" in reason
    assert pwh.certificate(length, relLengths) is not None


@pytest.mark.slow
def test_adding_a_vertex_at_either_end_keeps_the_certificate():
    """Corollary introducevertexatend, checked in the direction it is used.

    Extending the paper's A_9 example by a vertex at either end must stay
    certified -- and the deletion recursion is what has to see it, since it is
    the contrapositive of the same corollary.
    """
    cache = {}
    for name in ("30330300", "03033030", "23033030"):
        relLengths = [int(c) for c in name]
        assert pwh._isAdmissible(10, relLengths), name
        assert pwh.notPiecewiseHereditaryByDeletion(10, relLengths, cache) is not None, name


@pytest.mark.slow
def test_deletion_reaches_more_than_the_direct_criteria():
    """It should be a strict improvement, or it is not worth the recursion."""
    cache = {}
    algebras = nk.LinearNakayamaAlgebra.allOfLength(10)
    direct = {a.className() for a in algebras
              if pwh.isNotPiecewiseHereditary(10, a.relLengths)}
    byDeletion = {a.className() for a in algebras
                  if pwh.notPiecewiseHereditaryByDeletion(10, a.relLengths, cache) is not None}
    assert direct < byDeletion


# -- comparing class names ------------------------------------------------

def test_a_quipu_and_its_own_canonical_type_are_not_two_classes():
    """A tame hereditary quipu is also derived equivalent to a canonical algebra.

    So 'P^(1,1)_(1,4,1)' and 'C(2,2,7)' can be two names for one class -- they
    are, at order 10 -- and treating differing name strings as proof of
    distinctness reports such a pair as separated when it is a merge candidate.
    """
    import mutationClassTable as mct

    quipu = nk.QuipuAlgebra((1, 4, 1), (1, 1))
    assert len(quipu.vertices()) == 10
    assert pwh.canonicalWeightType(10, quipu.coxeterPolynomial()) == (2, 2, 7)

    assert mct.formsAreCompatible("P^(1,1)_(1,4,1)", "C(2,2,7)")
    assert not mct.formsAreCompatible("P^(1,1)_(1,4,1)", "C(2,4,5)")
    assert not mct.formsAreCompatible("P^(1,4)_(1,0,2)", "P^(3,3)_(1,0,1)")
    # A quipu of finite representation type has no canonical weight type at all.
    assert not mct.formsAreCompatible("P^(0)_(0,9)", "C(2,2,7)")


def test_the_negative_certificate_never_settles_a_comparison():
    import mutationClassTable as mct

    assert mct.formsAreCompatible("P^(0)_(0,9)", mct.NOT_PIECEWISE_HEREDITARY)
    assert mct.formsAreCompatible(mct.NOT_PIECEWISE_HEREDITARY,
                                  mct.NOT_PIECEWISE_HEREDITARY)


def test_quipu_names_round_trip():
    import quipuForms as qf

    for k, m in [((1, 0, 1), (1, 4)), ((1, 1), (3,)), ((0, 9), (0,))]:
        name = qf.formatQuipu((k, m))
        assert qf.parseQuipuName(name) == (k, m)
    assert qf.parseQuipuName("C(2,2,7)") is None
    assert qf.parseQuipuName("not piecewise hereditary") is None
