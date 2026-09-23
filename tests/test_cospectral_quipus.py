"""Where the Coxeter polynomial stops separating derived equivalence classes.

The Coxeter polynomial of the path algebra of a tree is determined by the tree's
adjacency spectrum, so two *cospectral* non-isomorphic quipus give algebras that
are not derived equivalent yet share a Coxeter polynomial.  That is the whole
reason the classification cannot merge classes by Coxeter polynomial alone -- and
it can be mapped out completely, for any order, without running a single
mutation.
"""

import collections

import networkx as nx
import pytest
import sympy

from quivermutation import nakayama as nk
from quivermutation import quipuForms as qf

# Orders 1 to 12.  The paper's table gives 4, 6 and 11 for orders 6, 7 and 8.
QUIPU_COUNTS = [1, 1, 1, 2, 2, 4, 6, 11, 18, 36, 64, 127]


@pytest.mark.parametrize("order, expected", [
    (order, expected) if order < 12 else
    pytest.param(order, expected, marks = pytest.mark.slow)
    for order, expected in enumerate(QUIPU_COUNTS, start=1)])
def test_the_number_of_quipus_of_each_order(order, expected):
    assert len(qf.allQuipusOfOrder(order)) == expected


@pytest.mark.parametrize("order", range(1, 11))
def test_every_quipu_of_an_order_has_that_many_vertices_and_is_canonical(order):
    quipus = qf.allQuipusOfOrder(order)
    for k, m in quipus:
        graph = qf.graphFromQuipuParameters(k, m)
        assert graph.number_of_nodes() == order
        assert qf.quipuParameters(graph) == (k, m)     # already canonical
    # No two of them are the same tree.
    forms = {qf.canonicalTreeForm(qf.graphFromQuipuParameters(k, m)) for k, m in quipus}
    assert len(forms) == len(quipus)


@pytest.mark.parametrize("order", range(1, 10))
def test_the_two_routes_to_the_quipus_of_an_order_agree(order):
    """Enumerating parameters and enumerating trees must find the same quipus.

    `allQuipusOfOrder` enumerates the P^(m)_(k) parameter pairs and
    canonicalises them; `quipusByTreeEnumeration` enumerates the non-isomorphic
    trees of the order and keeps the ones whose degrees say they are quipus.
    Nothing about the notation enters the second one's choice of trees, so
    agreement checks the parameter enumeration and both readings of the
    definition of a quipu at once.
    """
    assert qf.allQuipusOfOrder(order) == qf.quipusByTreeEnumeration(order)


@pytest.mark.slow
@pytest.mark.parametrize("order", [10, 11, 12])
def test_the_two_routes_agree_further_out(order):
    assert qf.allQuipusOfOrder(order) == qf.quipusByTreeEnumeration(order)


def test_the_degree_test_and_the_main_string_test_agree_on_non_quipus():
    """Both readings must reject the same trees, not just accept the same ones."""
    star = nx.Graph([(1, 2), (1, 3), (1, 4), (1, 5)])                # degree 4
    # A degree-3 vertex with a degree-3 vertex down each of its three branches.
    # Four branch vertices is the smallest number that can fail to lie on one
    # path: with three they always do, since a path through two of them passes
    # through the third.
    offMainString = nx.Graph([(1, 2), (1, 3), (1, 4), (2, 5), (2, 6),
                              (3, 7), (3, 8), (4, 9), (4, 10)])
    for graph in (star, offMainString):
        assert not qf.isQuipuByDegrees(graph)
        assert not qf.isQuipu(graph)


@pytest.mark.parametrize("order", range(4, 9))
def test_no_quipu_of_order_at_most_8_is_cospectral_with_another(order):
    """This is why the published n <= 8 table is clean.

    Up to order 8 the Coxeter polynomial separates every class, so grouping by it
    gives the right answer there and only there.
    """
    assert qf.cospectralQuipuGroups(order) == {}


def test_the_first_collision_is_at_order_9():
    """P^(1,4)_(1,0,1) and P^(1,2)_(1,1,2) are cospectral but not isomorphic.

    They correspond to the LNAs A_{9,(1,3)}^{(3,6)} and A_{9,(1,4)}^{(3,4)},
    which are therefore *not* derived equivalent despite having the same Coxeter
    polynomial.  The classification's merge report flags this pair as
    'separated', and it is the smallest such pair.
    """
    groups = qf.cospectralQuipuGroups(9)
    assert len(groups) == 1
    (quipus,) = groups.values()
    assert {qf.formatQuipu(q) for q in quipus} == {"P^(1,4)_(1,0,1)", "P^(1,2)_(1,1,2)"}

    first, second = (nk.QuipuAlgebra(*q) for q in quipus)
    assert not nx.is_isomorphic(first.underlyingGraph(), second.underlyingGraph())
    assert first.canonicalForm() != second.canonicalForm()
    assert sympy.expand(first.coxeterPolynomial()) == sympy.expand(second.coxeterPolynomial())

    lnas = {q.correspondingLNA().className() for q in (first, second)}
    assert lnas == {"3060000", "3004000"}


@pytest.mark.parametrize("order, expected", [(9, 1), (10, 2), (11, 4)])
def test_how_many_collisions_each_order_has(order, expected):
    assert len(qf.cospectralQuipuGroups(order)) == expected


@pytest.mark.slow
@pytest.mark.parametrize("order", range(4, 12))
def test_cospectral_is_exactly_the_same_coxeter_polynomial(order):
    """The two ways of finding the collisions must give the same groups.

    Cospectrality is a statement about the undirected tree and costs one
    characteristic polynomial; equal Coxeter polynomials is a statement about the
    algebra and costs a Cartan matrix inverse.  They agree for every order
    checked, which means the cheap computation maps the failure of the Coxeter
    polynomial completely.
    """
    quipus = qf.allQuipusOfOrder(order)
    byCoxeter = collections.defaultdict(list)
    for parameters in quipus:
        polynomial = sympy.expand(nk.QuipuAlgebra(*parameters).coxeterPolynomial())
        byCoxeter[polynomial].append(parameters)

    coxeterGroups = {tuple(sorted(g)) for g in byCoxeter.values() if len(g) > 1}
    cospectralGroups = {tuple(sorted(g)) for g in qf.cospectralQuipuGroups(order).values()}
    assert coxeterGroups == cospectralGroups
