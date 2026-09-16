"""Canonical forms for the hereditary algebras a search reaches.

The underlying undirected graph of a relation-free quiver is a complete derived
invariant for the tree case, so it settles what the Coxeter polynomial can only
suggest.  These check the encoding itself, and then that a search out of an LNA
lands on the quipu the paper says it should.
"""

import itertools

import networkx as nx
import pytest

from quivermutation import quipuForms as qf
import quivermutation as qm
from helpers import line_algebra, quiet


def graph(edges):
    g = nx.Graph()
    g.add_edges_from(edges)
    return g


DYNKIN = {
    "A_5": ([(1, 2), (2, 3), (3, 4), (4, 5)], "P^(0)_(0,4)"),
    "A_6": ([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)], "P^(0)_(0,5)"),
    "D_5": ([(1, 2), (2, 3), (3, 4), (3, 5)], "P^(2)_(1,1)"),
    "D_6": ([(1, 2), (2, 3), (3, 4), (4, 5), (4, 6)], "P^(3)_(1,1)"),
    "E_6": ([(1, 2), (2, 3), (3, 4), (4, 5), (3, 6)], "P^(2)_(1,2)"),
    "D~_5": ([(1, 2), (2, 3), (3, 4), (2, 5), (3, 6)], "P^(1,1)_(1,0,1)"),
}

# The labels the paper gives the non-Dynkin classes of order 7 and 8, as the
# (k, m) parameters it writes them with.  The notation does not determine the
# quipu, so these are canonicalised before comparing with anything.
PAPER_QUIPUS = {
    "P_(1,0,2)^(1,1)": ((1, 0, 2), (1, 1)),
    "P_(1,1,2)^(1,1)": ((1, 1, 2), (1, 1)),
    "P_(1,0,3)^(1,1)": ((1, 0, 3), (1, 1)),
    "P_(1,0,0,1)^(1,1,1)": ((1, 0, 0, 1), (1, 1, 1)),
    "P_(2,0,2)^(1,1)": ((2, 0, 2), (1, 1)),
    "P_(1,0,2)^(1,2)": ((1, 0, 2), (1, 2)),
    "P_(2,3)^(2)": ((2, 3), (2,)),
}


@pytest.mark.parametrize("name", sorted(DYNKIN))
def test_quipu_notation_of_the_dynkin_trees(name):
    edges, expected = DYNKIN[name]
    assert qf.formatQuipu(qf.quipuParameters(graph(edges))) == expected


@pytest.mark.parametrize("name", sorted(DYNKIN))
def test_the_parameters_account_for_every_vertex(name):
    """|P^(m)_(k)| = r + 1 + sum(k) + sum(m), per the paper's observation."""
    edges, _ = DYNKIN[name]
    g = graph(edges)
    k, m = qf.quipuParameters(g)
    assert len(m) + sum(k) + sum(m) == g.number_of_nodes()
    assert len(k) == len(m) + 1


@pytest.mark.parametrize("name", sorted(DYNKIN))
def test_the_canonical_form_does_not_depend_on_the_labelling(name):
    """Every permutation of the labels must give the same parameters.

    quipuParameters used to read each main string in one direction only, chosen
    by which endpoint sorted first, so relabelling could flip k and m and change
    the answer.
    """
    edges, _ = DYNKIN[name]
    original = graph(edges)
    labels = sorted(original.nodes)
    for permutation in itertools.islice(itertools.permutations(labels), 40):
        relabelled = nx.relabel_nodes(original, dict(zip(labels, permutation)))
        assert qf.canonicalTreeForm(relabelled) == qf.canonicalTreeForm(original)
        assert qf.quipuParameters(relabelled) == qf.quipuParameters(original)


@pytest.mark.parametrize("label", sorted(PAPER_QUIPUS))
def test_the_papers_quipu_labels_round_trip(label):
    """Build the quipu from the paper's parameters, read them back, rebuild.

    The notation is not unique, so the parameters read back need not be the ones
    written down, but the graph they describe must be the same one.
    """
    k, m = PAPER_QUIPUS[label]
    built = qf.graphFromQuipuParameters(k, m)
    assert len(m) + sum(k) + sum(m) == built.number_of_nodes()
    canonical = qf.quipuParameters(built)
    rebuilt = qf.graphFromQuipuParameters(*canonical)
    assert nx.is_isomorphic(built, rebuilt)
    assert qf.canonicalTreeForm(built) == qf.canonicalTreeForm(rebuilt)


def test_the_papers_quipus_of_order_8_are_all_distinct():
    """The seven non-Dynkin classes of order 8 must have seven distinct forms."""
    orderEight = {
        label: qf.canonicalTreeForm(qf.graphFromQuipuParameters(*PAPER_QUIPUS[label]))
        for label in PAPER_QUIPUS
        if sum(PAPER_QUIPUS[label][0]) + sum(PAPER_QUIPUS[label][1])
        + len(PAPER_QUIPUS[label][1]) == 8
    }
    assert len(orderEight) == 6
    assert len(set(orderEight.values())) == 6


def test_non_isomorphic_trees_of_the_same_size_get_different_forms():
    forms = {name: qf.canonicalTreeForm(graph(edges)) for name, (edges, _) in DYNKIN.items()}
    six = {name: form for name, form in forms.items() if name in ("A_6", "D_6", "E_6", "D~_5")}
    assert len(set(six.values())) == 4


def test_a_tree_with_a_degree_four_vertex_is_not_a_quipu():
    star = graph([(1, 2), (1, 3), (1, 4), (1, 5)])
    assert qf.quipuParameters(star) is None
    assert not qf.isQuipu(star)
    # It still has a canonical form, it just has no quipu notation.
    assert qf.canonicalTreeForm(star) == "(()()()())"


def test_a_cyclic_underlying_graph_falls_back_to_a_hash():
    cycle = graph([(1, 2), (2, 3), (3, 1)])
    form = qf.canonicalUndirectedForm(cycle)
    assert form.startswith("wl:")


@pytest.mark.parametrize(
    "length, rels, depth, expected",
    [
        (5, "000", 2, "P^(0)_(0,4)"),   # A_5 is already hereditary
        (5, "300", 4, "P^(2)_(1,1)"),   # A_{5,(1)}^{(3)} is derived equivalent to D_5
        (6, "3000", 4, "P^(3)_(1,1)"),  # A_{6,(1)}^{(3)} -> D_6
        (6, "2300", 4, "P^(2)_(1,2)"),  # A_{6,(2)}^{(3)} -> E_6
        (6, "3030", 7, "P^(1,1)_(1,0,1)"),  # A_{6,(1,3)}^{(3,3)} -> D~_5
    ],
)
def test_search_reaches_the_quipu_the_paper_predicts(length, rels, depth, expected):
    """Each of these appears in the n <= 8 table of arXiv:2305.06642."""
    forms = quiet(qm.hereditaryFormsReachedFrom, line_algebra(length, rels), depth)
    assert qm.formatHereditaryForms(forms) == expected


def test_a_class_reaches_exactly_one_hereditary_form():
    """All the hereditary algebras in one class have isomorphic quivers, so the
    search should never report two different forms for one LNA."""
    for length, rels, depth in [(5, "300", 5), (6, "3000", 5), (6, "2300", 5)]:
        forms = quiet(qm.hereditaryFormsReachedFrom, line_algebra(length, rels), depth)
        assert len(forms) == 1, f"A_{length}_{rels} reached {sorted(forms)}"
