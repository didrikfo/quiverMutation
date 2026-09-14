"""Cross-checks of this code against the published classification."""

import itertools

import pytest

from helpers import coxeter_poly, dynkin_A_coxeter, dynkin_D_coxeter, line_algebra
from paper_classification import PAPER_CLASSES, rel_lengths


def _algebra(length, entry):
    starts, lengths = entry
    return line_algebra(length, rel_lengths(length, starts, lengths))


@pytest.mark.parametrize("length", sorted(PAPER_CLASSES))
def test_members_of_a_paper_class_share_a_coxeter_polynomial(length):
    """The Coxeter polynomial is a derived invariant, so it must be constant on
    each of the paper's derived equivalence classes."""
    for label, members in PAPER_CLASSES[length].items():
        polys = {coxeter_poly(_algebra(length, entry)) for entry in members}
        assert len(polys) == 1, f"n={length} class {label} split: {polys}"


@pytest.mark.parametrize("length", sorted(PAPER_CLASSES))
def test_distinct_paper_classes_have_distinct_coxeter_polynomials(length):
    """For n <= 8 the Coxeter polynomial happens to separate every class.

    This is not true in general -- two LNAs that are not derived equivalent can
    share a Coxeter polynomial -- which is why merging classes in the generated
    CSV needs more than the polynomial.  But up to n = 8 it does separate them,
    so a failure here means a computed polynomial is wrong.
    """
    polys = {
        label: coxeter_poly(_algebra(length, members[0]))
        for label, members in PAPER_CLASSES[length].items()
    }
    for (a, pa), (b, pb) in itertools.combinations(polys.items(), 2):
        assert pa != pb, f"n={length}: classes {a} and {b} share {pa}"


@pytest.mark.parametrize("length", range(2, 9))
def test_linear_quiver_has_the_dynkin_A_coxeter_polynomial(length):
    assert coxeter_poly(line_algebra(length, [0] * (length - 2))) == dynkin_A_coxeter(length)


@pytest.mark.parametrize("length", range(4, 9))
def test_relations_of_length_two_do_not_change_the_class(length):
    """Adding relations of length 2 to A_n keeps it derived equivalent to A_n.

    This is the first of the paper's three class-preserving operations, and the
    reason the published table only lists relations of length >= 3.
    """
    all_two = [2] * (length - 2)
    assert coxeter_poly(line_algebra(length, all_two)) == dynkin_A_coxeter(length)


@pytest.mark.parametrize("length", range(4, 9))
def test_the_D_class_has_the_dynkin_D_coxeter_polynomial(length):
    """A_{n,(1)}^{(3)} is derived equivalent to the path algebra of D_n."""
    starts, lengths = PAPER_CLASSES[length][f"D_{length}"][0]
    pa = line_algebra(length, rel_lengths(length, starts, lengths))
    assert coxeter_poly(pa) == dynkin_D_coxeter(length)
