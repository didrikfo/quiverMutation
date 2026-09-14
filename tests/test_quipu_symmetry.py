"""Which quipus are the same quipu, and which LNAs are therefore one class.

The quipu notation P^(m_0,...,m_r)_(k_0,...,k_{r+1}) does not determine the
graph.  Two re-readings give the same tree:

* reading the main string from the other end, which reverses both tuples;
* exchanging an *end* segment of the main string with the cord at the outermost
  foot -- `k_0` with `m_0`, or `k_{r+1}` with `m_r` -- because those are the two
  branches hanging off that foot.

Both are isomorphisms of the tree, so the LNAs they name are derived equivalent,
and the classification must not separate them.  The exchange does **not** extend
to an interior gap `k_i`, `0 < i < r + 1`: there the foot has a third branch
running on along the main string, so the two are not interchangeable.

Getting this wrong in either direction is the expensive mistake in this project.
Missing the end exchange splits one class in two; applying it at an interior gap,
or merging on the Coxeter polynomial instead, merges two classes into one -- which
is what the hand-made n = 9 workbook did (R-006, F-014).

These tests pin the symmetry down from both sides and check that the two sides
agree: `quipuForms.quipuParameters`, which canonicalises the notation by
searching every reading of the graph, and the explicit operations of
`cor:EquivNakayamaAlgebras` on the LNA, which are implemented independently of it
in `nakayama.LinearNakayamaAlgebra`.
"""

import collections
import itertools

import networkx as nx
import pytest

import nakayama as nk
import quipuForms as qf


# -- every name of a given order ------------------------------------------

def quipuNames(order):
    """Every (k, m) parameter pair describing a quipu on `order` vertices.

    Not up to isomorphism: this is the raw notation, many names per tree, which
    is the point -- the tests below are about which of them name the same tree.
    """
    names = []
    for cords in range(1, order + 1):
        budget = order - cords
        if budget < 0:
            continue
        for m in qf._compositions(cords, budget):
            for k in qf._compositions(cords + 1, budget - sum(m)):
                if sum(k) + sum(m) + cords == order:
                    names.append((k, m))
    return names


def exchange(k, m, gap, cord):
    """Exchange the main-string gap k[gap] with the cord m[cord].

    An isomorphism of the tree exactly when the gap and the cord are the only two
    branches at a foot, which happens at the two ends: (0, 0) and
    (len(k) - 1, len(m) - 1).
    """
    k, m = list(k), list(m)
    k[gap], m[cord] = m[cord], k[gap]
    return tuple(k), tuple(m)


def endExchanges(k, m):
    """The quipu's two end exchanges, as (k, m) pairs, including trivial ones."""
    return [exchange(k, m, 0, 0), exchange(k, m, len(k) - 1, len(m) - 1)]


# -- the end exchanges ----------------------------------------------------

@pytest.mark.parametrize("order", range(3, 10))
def test_the_end_exchanges_do_not_change_the_quipu(order):
    """k_0 <-> m_0 and k_{r+1} <-> m_r give the same tree, hence the same class."""
    checked = 0
    for k, m in quipuNames(order):
        original = qf.graphFromQuipuParameters(k, m)
        for k2, m2 in endExchanges(k, m):
            exchanged = qf.graphFromQuipuParameters(k2, m2)
            assert nx.is_isomorphic(original, exchanged), (k, m, k2, m2)
            assert qf.quipuParameters(original) == qf.quipuParameters(exchanged), (k, m, k2, m2)
            checked += 1
    assert checked > 0


def test_the_same_exchange_at_an_interior_gap_is_not_an_isomorphism():
    """Which is why it must not be applied there.

    At an interior foot the main string continues in both directions, so the gap
    and the cord are not the two branches of a degree-3 vertex -- there is a
    third.  P^(1,2,1)_(1,1,2,1) exchanged at k_1 gives a different tree of the
    same order.
    """
    k, m = (1, 1, 2, 1), (1, 2, 1)
    original = qf.graphFromQuipuParameters(k, m)
    for gap, cord in ((1, 1), (2, 2)):
        k2, m2 = exchange(k, m, gap, cord)
        interior = qf.graphFromQuipuParameters(k2, m2)
        assert interior.number_of_nodes() == original.number_of_nodes()
        assert not nx.is_isomorphic(original, interior), (k2, m2)
        assert qf.quipuParameters(original) != qf.quipuParameters(interior), (k2, m2)


def test_reading_the_main_string_backwards_does_not_change_the_quipu():
    for order in range(3, 9):
        for k, m in quipuNames(order):
            reversed_ = qf.graphFromQuipuParameters(k[::-1], m[::-1])
            assert qf.quipuParameters(qf.graphFromQuipuParameters(k, m)) == \
                qf.quipuParameters(reversed_), (k, m)


# -- canonicalisation is exactly isomorphism ------------------------------

@pytest.mark.parametrize("order", range(3, 10))
def test_canonical_parameters_are_exactly_tree_isomorphism(order):
    """Two names canonicalise together if and only if they are the same tree.

    Both directions matter and they fail differently.  If two names for one tree
    got different parameters, one class would be reported as two -- the bug the
    end exchange would cause if it were missed.  If two names for different trees
    got the same parameters, two classes would be reported as one.

    networkx decides isomorphism here, so this does not test the canonical form
    against itself.
    """
    byCanonical = collections.defaultdict(list)
    graphs = {}
    for name in quipuNames(order):
        graph = qf.graphFromQuipuParameters(*name)
        graphs[name] = graph
        byCanonical[qf.quipuParameters(graph)].append(name)

    for canonical, group in byCanonical.items():
        first = graphs[group[0]]
        for other in group[1:]:
            assert nx.is_isomorphic(first, graphs[other]), (canonical, group[0], other)

    representatives = {c: graphs[g[0]] for c, g in byCanonical.items()}
    for a, b in itertools.combinations(sorted(representatives), 2):
        assert not nx.is_isomorphic(representatives[a], representatives[b]), (a, b)

    assert sorted(byCanonical) == qf.allQuipusOfOrder(order)


@pytest.mark.slow
@pytest.mark.parametrize("order", [10, 11])
def test_canonical_parameters_are_exactly_tree_isomorphism_further_out(order):
    test_canonical_parameters_are_exactly_tree_isomorphism(order)


# -- the LNA side agrees with the quipu side ------------------------------

def longRelationLNAs(length):
    """Every LNA of that length with almost separate relations and no length-2 one.

    The setting of `cor:EquivNakayamaAlgebras`: relations of length 2 are dropped
    by its operation 2, so the orbits it describes live among these.
    """
    return [a for a in nk.LinearNakayamaAlgebra.allOfLength(length)
            if 2 not in a.relLengths and a.hasAlmostSeparateRelations()]


@pytest.mark.parametrize("length", range(4, 10))
def test_the_papers_operations_generate_exactly_the_quipu_fibre(length):
    """The two routes to the class must agree, exactly.

    `classPreservingOrbit` applies the paper's three operations to the relations;
    `quipu` names a tree and canonicalises it.  If the fibres of the naming were
    coarser than the orbits, the classification would be merging classes the
    paper keeps apart; if they were finer, it would be splitting a class the
    paper's own operations connect -- the failure mode that would make the n = 9
    count wrong.
    """
    fibres = collections.defaultdict(set)
    for a in longRelationLNAs(length):
        quipu = a.quipu()
        assert quipu is not None, a
        fibres[quipu].add(a)

    assert fibres, length
    for quipu, members in fibres.items():
        for a in members:
            assert a.classPreservingOrbit() == frozenset(members), (quipu, a)


@pytest.mark.slow
def test_the_papers_operations_generate_exactly_the_quipu_fibre_at_ten():
    test_the_papers_operations_generate_exactly_the_quipu_fibre(10)


@pytest.mark.parametrize("length", range(4, 10))
def test_a_class_holds_at_most_eight_of_these_algebras(length):
    """The paper's own bound on the orbit, which is a check on both sides."""
    for a in longRelationLNAs(length):
        assert len(a.classPreservingOrbit()) <= 8, a


@pytest.mark.parametrize("length", range(4, 10))
def test_every_quipu_of_an_order_is_named_by_one_of_these_algebras(length):
    """The naming is onto: no quipu class is missing its representative LNA."""
    named = {a.quipu() for a in longRelationLNAs(length)}
    assert named == set(qf.allQuipusOfOrder(length))


@pytest.mark.parametrize("length", range(4, 10))
def test_the_first_exchange_fixes_the_relations_end_and_the_last_its_start(length):
    """The direction of each exchange, which is easy to get backwards.

    At the first foot the relation's end vertex stays put and start and length
    trade places; at the last foot the start stays put and the length becomes
    what was left of the line.  Stating it the other way round would give a
    different, wrong, pairing of algebras.
    """
    checked = 0
    for a in longRelationLNAs(length):
        relations = a.relations()
        if not relations:
            continue

        first = a.swapFirstRelation()
        if first is not None and first != a:
            (start, arrows), (newStart, newArrows) = relations[0], first.relations()[0]
            assert start + arrows == newStart + newArrows, (a, first)
            assert (newStart, newArrows) == (arrows - 2, start + 2), (a, first)
            assert first.quipu() == a.quipu(), (a, first)
            checked += 1

        last = a.swapLastRelation()
        if last is not None and last != a:
            (start, arrows), (newStart, newArrows) = relations[-1], last.relations()[-1]
            assert newStart == start, (a, last)
            assert newArrows == length - start - arrows + 3, (a, last)
            assert last.quipu() == a.quipu(), (a, last)
            checked += 1
    if length >= 5:
        # At length 4 the only algebra here is A_{4,(1)}^{(3)}, which both
        # exchanges fix, so there is nothing to have got backwards.
        assert checked > 0


def test_a_length_two_relation_does_not_move_the_class():
    """Operation 2, on the algebra rather than on the quipu."""
    plain = nk.LinearNakayamaAlgebra(6, "3000")
    decorated = nk.LinearNakayamaAlgebra(6, "3002")
    assert decorated.withoutShortRelations() == plain
    assert decorated.quipu() == plain.quipu()
    assert decorated.classPreservingOrbit() == plain.classPreservingOrbit()


# -- the pair the workbook merged -----------------------------------------

DISPUTED = ("3060000", "3004000")


def test_the_cospectral_order_nine_pair_is_two_classes():
    """A_{9,(1,3)}^{(3,6)} and A_{9,(1,4)}^{(3,4)} are not derived equivalent.

    The hand-made workbook merged them, giving 19 classes at n = 9 where there
    are 20 (R-006).  Its reason was the Coxeter polynomial, which these two do
    share -- they are the smallest cospectral pair of quipus, F-010.

    Every symmetry that could merge them is checked here and none does: the trees
    have different diameters, so no relabelling relates them; neither end
    exchange moves either quipu onto the other; and the orbits of the paper's
    three operations are disjoint.
    """
    first, second = (nk.LinearNakayamaAlgebra(9, name) for name in DISPUTED)
    assert first.quipuName() == "P^(1,4)_(1,0,1)"
    assert second.quipuName() == "P^(1,2)_(1,1,2)"

    graphs = [qf.graphFromQuipuParameters(*a.quipu()) for a in (first, second)]
    assert [g.number_of_nodes() for g in graphs] == [9, 9]
    assert [nx.diameter(g) for g in graphs] == [6, 5]
    assert not nx.is_isomorphic(*graphs)

    # Cospectral, which is exactly why the Coxeter polynomial cannot see this.
    assert qf.adjacencySpectrumPolynomial(graphs[0]) == qf.adjacencySpectrumPolynomial(graphs[1])
    assert first.coxeterPolynomial() == second.coxeterPolynomial()

    assert first.classPreservingOrbit().isdisjoint(second.classPreservingOrbit())
    assert {a.className() for a in first.classPreservingOrbit()} == {
        "3060000", "3030000", "6000030", "0003030"}
    assert {a.className() for a in second.classPreservingOrbit()} == {"3004000", "0400030"}


def test_neither_end_exchange_relates_the_disputed_quipus():
    """Spelled out on the parameters, since this is the claim that was doubted.

    P^(1,4)_(1,0,1) exchanges to P^(1,1)_(1,0,4) at its last foot and is fixed at
    its first; P^(1,2)_(1,1,2) is fixed at both.  Neither reaches the other.
    """
    disputed = [nk.LinearNakayamaAlgebra(9, name).quipu() for name in DISPUTED]
    reachable = []
    for k, m in disputed:
        images = {qf.canonicalQuipuParameters(k2, m2) for k2, m2 in endExchanges(k, m)}
        reachable.append(images | {(k, m)})
    assert reachable[0] == {((1, 0, 1), (1, 4))}
    assert reachable[1] == {((1, 1, 2), (1, 2))}
    assert reachable[0].isdisjoint(reachable[1])


def test_the_same_split_is_already_in_the_published_table_at_n_equals_eight():
    """One vertex down, the same two shapes are two rows of the paper's own table.

    P_(1,0,3)^(1,1) and P_(1,1,2)^(1,1) are listed separately in the n <= 8
    classification of arXiv:2305.06642.  Adding a vertex to the last cord of the
    first and to the last cord of the second gives the n = 9 pair, so splitting
    that pair is not a new claim -- it is the published one, at the first order
    where the Coxeter polynomial can no longer see it.
    """
    from paper_classification import PAPER_CLASSES

    eight = PAPER_CLASSES[8]
    assert "P_(1,0,3)^(1,1)" in eight and "P_(1,1,2)^(1,1)" in eight
    a = nk.LinearNakayamaAlgebra(8, "303000")     # A_{8,(1,3)}^{(3,3)}
    b = nk.LinearNakayamaAlgebra(8, "300300")     # A_{8,(1,4)}^{(3,3)}
    assert a.quipu() == qf.canonicalQuipuParameters((1, 0, 3), (1, 1))
    assert b.quipu() == qf.canonicalQuipuParameters((1, 1, 2), (1, 1))
    assert a.quipu() != b.quipu()
    assert a.classPreservingOrbit().isdisjoint(b.classPreservingOrbit())
