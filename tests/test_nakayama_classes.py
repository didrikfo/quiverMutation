"""The LNA and quipu algebra classes, and the theorem that connects them."""

import pytest

from quivermutation import nakayama as nk
from quivermutation import quipuForms as qf
from helpers import coxeter_poly, dynkin_A_coxeter, dynkin_D_coxeter, quiet
from paper_classification import PAPER_CLASSES, rel_lengths


# -- LinearNakayamaAlgebra ------------------------------------------------

def test_the_quiver_and_relations_are_built_from_the_relation_lengths():
    a = nk.LinearNakayamaAlgebra(6, "2030")
    assert sorted(a.vertices()) == [1, 2, 3, 4, 5, 6]
    assert sorted((x[0], x[1]) for x in a.arrows()) == [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]
    assert a.rels == [[[1, 2, 3]], [[3, 4, 5, 6]]]
    assert a.relations() == [(1, 2), (3, 3)]


def test_the_three_names_agree():
    a = nk.LinearNakayamaAlgebra(7, "22300")
    assert a.className() == "22300"
    assert a.relationString() == "1;2;3|2;3;4|3;4;5;6"
    assert nk.LinearNakayamaAlgebra.fromClassName("22300") == a
    assert nk.LinearNakayamaAlgebra.fromRelationString(7, a.relationString()) == a


def test_the_one_and_two_vertex_lines_take_no_relation_lengths():
    """A_1 and A_2 admit no relation, so n - 2 is not the count for them.

    This was the one place the refactor of the line helpers found something
    wrong rather than just moving it: `lineQuiverExample(1, [])` printed an
    error and returned an algebra with no vertices at all, and the n = 1 row of
    the paper's table was being checked against that.
    """
    for length in (1, 2):
        a = nk.LinearNakayamaAlgebra(length, [])
        assert sorted(a.vertices()) == list(range(1, length + 1))
        assert a.rels == []
        assert a.relations() == []
        assert a.className() == ""
        assert a.kupischSeries() == tuple(range(length, 0, -1))
    with pytest.raises(ValueError):
        nk.LinearNakayamaAlgebra(0, [])
    with pytest.raises(ValueError):
        nk.LinearNakayamaAlgebra(1, [2])


@pytest.mark.parametrize(
    "length, rels",
    [(5, "900"), (5, "10"), (4, "1"), (6, "0005")],
)
def test_a_relation_that_does_not_fit_is_rejected(length, rels):
    with pytest.raises(ValueError):
        nk.LinearNakayamaAlgebra(length, rels)


@pytest.mark.parametrize("length", [3, 4, 5, 6, 7])
def test_there_are_catalan_many_of_each_length(length):
    import math
    catalan = math.comb(2 * (length - 1), length - 1) // length
    algebras = nk.LinearNakayamaAlgebra.allOfLength(length)
    assert len(algebras) == catalan
    assert len(set(algebras)) == catalan       # hashable and distinct


@pytest.mark.parametrize("length", [4, 5, 6, 7])
def test_the_kupisch_series_is_the_column_sums_of_the_cartan_matrix(length):
    """dim P_i is the i-th column sum, since column i counts the paths out of i."""
    for a in nk.LinearNakayamaAlgebra.allOfLength(length):
        cartan = a.cartanMatrix()
        columnSums = tuple(sum(cartan[row, col] for row in range(length)) for col in range(length))
        assert a.kupischSeries() == columnSums, a


@pytest.mark.parametrize("length", [4, 5, 6, 7])
def test_the_kupisch_series_is_admissible(length):
    """c_i <= c_{i+1} + 1, and c_n = 1: the classical condition."""
    for a in nk.LinearNakayamaAlgebra.allOfLength(length):
        series = a.kupischSeries()
        assert series[-1] == 1
        assert all(c >= 1 for c in series)
        assert all(series[i] <= series[i + 1] + 1 for i in range(length - 1))


@pytest.mark.parametrize("length", [4, 5, 6, 7])
def test_the_relation_dual_is_an_involution_and_preserves_the_class(length):
    for a in nk.LinearNakayamaAlgebra.allOfLength(length):
        dual = a.relationDual()
        assert dual.relationDual() == a
        assert dual.coxeterPolynomial() == a.coxeterPolynomial()


# -- QuipuAlgebra and the theorem -----------------------------------------

def test_a_quipu_algebra_has_no_relations_and_the_right_shape():
    q = nk.QuipuAlgebra((1, 1), (2,))
    assert q.rels == []
    assert len(q.vertices()) == 5
    assert q.canonicalForm() == qf.canonicalTreeForm(
        qf.graphFromQuipuParameters((1, 1), (2,)))


@pytest.mark.parametrize("length", [n for n in sorted(PAPER_CLASSES) if n >= 2])
def test_an_lna_and_its_quipu_have_the_same_coxeter_polynomial(length):
    """The whole theorem chain, checked end to end.

    The LNA is derived equivalent to the quipu algebra, so the two must have the
    same Coxeter polynomial.  The LNA's is computed from its own Cartan matrix
    and the quipu's from the tree's, by entirely separate routes.
    """
    for label, entries in PAPER_CLASSES[length].items():
        for starts, lengths in entries:
            a = nk.LinearNakayamaAlgebra(length, rel_lengths(length, starts, lengths))
            q = nk.QuipuAlgebra.fromLNA(a)
            assert q.coxeterPolynomial() == a.coxeterPolynomial(), (label, a)


@pytest.mark.parametrize("length", [4, 5, 6, 7])
def test_every_lna_with_almost_separate_relations_matches_its_quipu(length):
    """Not just the ones the paper tabulates: every LNA the theorem covers."""
    checked = 0
    for a in nk.LinearNakayamaAlgebra.allOfLength(length):
        if not a.hasAlmostSeparateRelations():
            continue
        q = nk.QuipuAlgebra.fromLNA(a)
        assert len(q.vertices()) == length, a
        assert q.coxeterPolynomial() == a.coxeterPolynomial(), a
        checked += 1
    assert checked > 0


@pytest.mark.parametrize("length", [4, 5, 6, 7])
def test_the_quipu_round_trips_back_to_an_equivalent_lna(length):
    """quipu -> LNA -> quipu is the identity on the quipu."""
    for a in nk.LinearNakayamaAlgebra.allOfLength(length):
        if not a.hasAlmostSeparateRelations():
            continue
        q = nk.QuipuAlgebra.fromLNA(a)
        back = q.correspondingLNA()
        assert back.quipu() == q.quipuParameters(), (a, back)
        assert back.coxeterPolynomial() == a.coxeterPolynomial(), (a, back)


def test_an_lna_outside_the_theorem_refuses_to_name_a_quipu():
    overlapping = nk.LinearNakayamaAlgebra(6, "3300")   # two 3-arrow relations sharing two arrows
    assert not overlapping.hasAlmostSeparateRelations()
    assert overlapping.quipu() is None
    with pytest.raises(ValueError):
        nk.QuipuAlgebra.fromLNA(overlapping)


@pytest.mark.parametrize("length", [4, 5, 6, 7, 8])
def test_the_line_without_relations_is_the_dynkin_A_quipu(length):
    a = nk.LinearNakayamaAlgebra(length, [0] * (length - 2))
    assert a.quipuName() == "P^(0)_(0,{0})".format(length - 1)
    assert a.coxeterPolynomial() == dynkin_A_coxeter(length)


@pytest.mark.parametrize("length", [4, 5, 6, 7, 8])
def test_the_D_class_representative(length):
    a = nk.LinearNakayamaAlgebra(length, [3] + [0] * (length - 3))
    assert a.quipuName() == "P^({0})_(1,1)".format(length - 3)
    assert a.coxeterPolynomial() == dynkin_D_coxeter(length)
