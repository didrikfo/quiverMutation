"""What the hereditary-form column is allowed to decide.

Three values live in that column and they carry different weight.  A quipu name
or a tree encoding is a *proof*: a mutation path reached a relation-free quiver,
or the quipu theorem handed the name over.  The negative certificate is a proof
too, of the opposite thing.  A `C(...)` canonical weight type is neither -- it is
read off the class' own Coxeter polynomial, so within a group of classes sharing
one polynomial it says the same thing about every member and separates nothing.

Confusing the third kind with the first is what put the n = 9 classification at
22 classes instead of 20 (F-018), so this pins the distinction at every level it
matters: the predicates, the merge report, and one real pair from that run.
"""

import pytest

from quivermutation import classification as cl
from quivermutation import mutationClassTable as mct
from quivermutation import piecewiseHereditary as ph
from quivermutation import quipuForms as qf

POLY = "shared"
NUMBERING = "1;2;3;4;5;6;7;8;9"


def tableOf(rows, length = 9):
    """A table from (relations, class, coxeter, form) tuples."""
    return mct.MutationClassTable(
        [[relations, className, "", coxeter, NUMBERING, form]
         for relations, className, coxeter, form in rows], length)


# -- which values are proofs ----------------------------------------------

@pytest.mark.parametrize("form, identifying, coxeterDerived, proved", [
    ("P^(5)_(1,2)", True, False, True),
    ("(((()())())())", True, False, True),
    ("C(2,4,4)", True, True, False),
    ("C(2,3,5)", True, True, False),
    (mct.NOT_PIECEWISE_HEREDITARY, False, False, False),
    ("", False, False, False),
])
def test_what_a_form_value_proves(form, identifying, coxeterDerived, proved):
    assert mct.isIdentifyingForm(form) == identifying
    assert mct.isCoxeterDerivedForm(form) == coxeterDerived
    assert mct.isProvedForm(form) == proved


def test_two_classes_are_not_merged_on_a_shared_canonical_type():
    """Both read the same weight type off the same polynomial, so this is
    circular: all it says is that the two share a Coxeter polynomial."""
    table = tableOf([
        ("1;2;3", "A", POLY, "C(2,3,5)"),
        ("2;3;4", "B", POLY, "C(2,3,5)"),
    ])
    assert table.classesByHereditaryForm() == {}
    assert cl.mergeReport(table)["certain"] == {}


def test_two_classes_are_merged_on_a_shared_quipu():
    table = tableOf([
        ("1;2;3", "A", POLY, "P^(5)_(1,2)"),
        ("2;3;4", "B", POLY, "P^(5)_(1,2)"),
    ])
    assert cl.mergeReport(table)["certain"] == {"P^(5)_(1,2)": {"A", "B"}}


# -- what the merge report may separate -----------------------------------

def test_a_canonical_type_never_separates_a_class_from_a_quipu_class():
    """The bug of F-018, at the level it lived at.

    Class B has no name of its own, so the pipeline reads a weight type off its
    Coxeter polynomial -- which is A's polynomial.  The two strings then differ,
    and the old report concluded the classes were distinct.  The right answer is
    that nothing has been decided and the pair needs a search.
    """
    table = tableOf([
        ("1;2;3", "P^(5)_(1,2)", POLY, "P^(5)_(1,2)"),
        ("2;3;4", "2334400", POLY, "C(2,3,5)"),
    ])
    report = cl.mergeReport(table)
    assert report["separated"] == {}
    assert report["candidate"] == {POLY: {"P^(5)_(1,2)", "2334400"}}


def test_two_different_quipus_are_separated():
    """The F-010 situation: cospectral quipus, so one polynomial and two classes."""
    table = tableOf([
        ("1;2;3", "A", POLY, "P^(1,4)_(1,0,1)"),
        ("2;3;4", "B", POLY, "P^(1,2)_(1,1,2)"),
    ])
    report = cl.mergeReport(table)
    assert report["separated"] == {POLY: {"A", "B"}}
    assert report["candidate"] == {}


def test_the_negative_certificate_separates_a_class_from_a_quipu_class():
    """Unlike a weight type this one is proved, and it is a real separation: a
    class that is not piecewise hereditary is in no quipu class."""
    table = tableOf([
        ("1;2;3", "A", POLY, "P^(5)_(1,2)"),
        ("2;3;4", "B", POLY, mct.NOT_PIECEWISE_HEREDITARY),
    ])
    assert cl.mergeReport(table)["separated"] == {POLY: {"A", "B"}}


def test_two_classes_that_are_only_known_not_to_be_quipus_stay_candidates():
    table = tableOf([
        ("1;2;3", "A", POLY, mct.NOT_PIECEWISE_HEREDITARY),
        ("2;3;4", "B", POLY, mct.NOT_PIECEWISE_HEREDITARY),
    ])
    report = cl.mergeReport(table)
    assert report["separated"] == {}
    assert report["candidate"] == {POLY: {"A", "B"}}


# -- the affine tree a domestic weight type names -------------------------

@pytest.mark.parametrize("weights, affine, quipu", [
    ((2, 3, 3), "E~6", "P^(2)_(2,2)"),
    ((2, 3, 4), "E~7", "P^(3)_(1,3)"),
    ((2, 3, 5), "E~8", "P^(5)_(1,2)"),
    ((2, 2, 6), "D~8", "P^(1,1)_(1,3,1)"),
    ((2, 2, 3), "D~5", "P^(1,1)_(1,0,1)"),
])
def test_a_domestic_weight_type_is_a_quipu_class(weights, affine, quipu):
    """Why a domestic `C(...)` is always an unfound merge, computed rather than
    asserted.

    A canonical algebra of domestic weight type is derived equivalent to the
    path algebra of the corresponding extended Dynkin diagram, and for these
    types that diagram is a tree -- in fact a quipu.  So the class is one of the
    quipu classes, which the quipu theorem has already named under another name
    in the same table.  The two n = 9 cases are (2, 3, 5) and (2, 2, 6), and the
    quipus they come out at are exactly the two classes the run wrongly split
    them from.
    """
    assert ph.isDomestic(weights)
    assert not ph.isTubular(weights)
    assert ph.affineTypeOfDomesticWeightType(weights) == affine
    tree = ph.affineTreeOfDomesticWeightType(weights)
    assert tree.number_of_nodes() == ph.canonicalVertexCount(weights)
    assert tree.number_of_edges() == tree.number_of_nodes() - 1      # a tree
    assert qf.formatQuipu(qf.quipuParameters(tree)) == quipu


@pytest.mark.parametrize("weights", [(2, 4, 4), (3, 3, 3), (2, 3, 6), (2, 2, 2, 2)])
def test_a_tubular_weight_type_is_not_a_quipu_class(weights):
    """The other side: tubular is the boundary past which the canonical algebra
    is derived equivalent to no hereditary algebra, so `C(2,4,4)` at n = 9 is a
    class of its own and not a merge waiting to happen."""
    assert ph.isTubular(weights)
    assert not ph.isDomestic(weights)
    assert ph.affineTypeOfDomesticWeightType(weights) is None
    assert ph.affineTreeOfDomesticWeightType(weights) is None


def test_two_weights_are_domestic_but_not_a_tree():
    """A~ is a cycle, so `C(p,q)` names a hereditary class that is not a quipu
    one.  Naming the type without offering a tree is the honest answer."""
    assert ph.isDomestic((2, 5))
    assert ph.affineTypeOfDomesticWeightType((2, 5)) == "A~6"
    assert ph.affineTreeOfDomesticWeightType((2, 5)) is None


def test_the_n_9_weight_types_are_read_off_the_right_polynomials():
    """`canonicalWeightType` really does return these two for those classes, so
    the pair above is the pair the pipeline meets."""
    import sympy

    x = sympy.Symbol("x")
    for weights in [(2, 3, 5), (2, 2, 6), (2, 4, 4)]:
        polynomial = str(ph.canonicalCoxeterPolynomial(weights, x)).replace("x", "lambda")
        assert ph.canonicalWeightType(9, polynomial) == weights


# -- the real pair, end to end --------------------------------------------

def test_a_class_named_only_by_its_polynomial_is_still_searched():
    """The n = 9 merge the run missed, at the depth it was always available at.

    `1;2;3;4;5|2;3;4;5;6;7` is a member of the class the run called 2334400 and
    named C(2,3,5); two mutations take it to a member of P^(5)_(1,2).  The old
    resolution step never looked, because it only searched classes with *no*
    form at all and this one had a form that proved nothing.
    """
    table = tableOf([
        ("1;2;3;4;5|2;3;4;5;6;7", "2334400", POLY, "C(2,3,5)"),
        ("1;2;3;4;5|2;3;4;5;6;7|6;7;8", "P^(5)_(1,2)", POLY, "P^(5)_(1,2)"),
    ])
    merges = cl.resolveMergeCandidates(table, 9, 2, printOutput = False)
    assert merges == [("2334400", "P^(5)_(1,2)")]
    assert table.classNames() == {"P^(5)_(1,2)"}
    # The merged class keeps the form that proves something, on every row.
    assert {row[5] for row in table.rows()} == {"P^(5)_(1,2)"}
    assert cl.mergeReport(table)["candidate"] == {}


def test_the_other_n_9_pair_merges_too():
    """`2;3;4;5|3;4;5;6|6;7;8;9`, of the class the run called 2233030 and named
    C(2,2,6), reaches a member of P^(1,1)_(1,3,1) in two mutations."""
    table = tableOf([
        ("2;3;4;5|3;4;5;6|6;7;8;9", "2233030", POLY, "C(2,2,6)"),
        ("1;2;3;4|2;3;4;5|6;7;8;9", "P^(1,1)_(1,3,1)", POLY, "P^(1,1)_(1,3,1)"),
    ])
    assert cl.resolveMergeCandidates(table, 9, 2, printOutput = False) == \
        [("2233030", "P^(1,1)_(1,3,1)")]
    assert table.classNames() == {"P^(1,1)_(1,3,1)"}
