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

import nakayama as nk
import quipuForms as qf

# Orders 1 to 12.  The paper's table gives 4, 6 and 11 for orders 6, 7 and 8.
QUIPU_COUNTS = [1, 1, 1, 2, 2, 4, 6, 11, 18, 36, 64, 127]


@pytest.mark.parametrize("order, expected", list(enumerate(QUIPU_COUNTS, start=1)))
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
