"""The integer route to the Coxeter polynomial, against the symbolic one.

`invariants.coxeterPoly` builds the Cartan matrix, inverts it over the rationals
and hands sympy a symbolic characteristic polynomial.  `coxeterCoefficients`
computes the same polynomial as `det(lambda C^T + C)`, evaluated at n + 1 integer
points and interpolated, which needs no inversion and no symbols and is what the
searches over millions of algebras use.  The two must agree everywhere, or the
searches are answering a different question from the classification.
"""

import networkx as nx
import pytest
import sympy

from quivermutation import coxeterTables as ct
from quivermutation import invariants as inv
from quivermutation import nakayama as nk
from quivermutation import quipuRelations as qr
from quivermutation import treeSearch as ts


def asExpression(coefficients):
    return sympy.expand(inv.polynomialFromCoefficients(coefficients))


@pytest.mark.parametrize("length", [4, 5, 6, 7])
def test_the_integer_route_agrees_with_the_symbolic_one_on_every_lna(length):
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        symbolic = sympy.expand(inv.coxeterPoly(algebra).as_expr())
        assert asExpression(inv.coxeterKey(algebra)) == symbolic


@pytest.mark.parametrize("length", [4, 5, 6, 7, 8])
def test_the_lna_cartan_shortcut_agrees_with_counting_paths(length):
    """`coxeterTables` builds an LNA's Cartan matrix from the relation lengths."""
    for algebra in nk.LinearNakayamaAlgebra.allOfLength(length):
        assert ct.lnaCartanMatrix(length, tuple(algebra.relLengths)) == \
            inv.integerCartanMatrix(algebra)


@pytest.mark.parametrize("order", [5, 6, 7])
def test_a_tree_gets_the_same_polynomial_however_it_is_oriented(order):
    """The orientations of a tree are derived equivalent, by BGP reflection."""
    for graph in ts.treesOfOrder(order):
        fromReachability = ts.treeCoxeterKey(graph)
        algebra = ts.orientedTreeQuiver(graph)
        assert asExpression(fromReachability) == sympy.expand(
            inv.coxeterPoly(algebra).as_expr())
        # A different root gives a different orientation and the same answer.
        for root in list(graph.nodes)[1:4]:
            other = nx.relabel_nodes(graph, _swap(min(graph.nodes), root))
            assert ts.treeCoxeterKey(other) == fromReachability


def _swap(first, second):
    return {first: second, second: first}


def test_the_determinant_identity_needs_a_unimodular_cartan_matrix():
    with pytest.raises(ValueError):
        inv.coxeterCoefficients([[2, 0], [0, 1]])


@pytest.mark.parametrize("order", [6, 7])
def test_the_cartan_bitmask_agrees_with_counting_paths(order):
    """`quipuRelations` reads a Cartan matrix off a bitmask instead of counting.

    The bitmask is maintained incrementally as relations are added, which is what
    makes the enumeration affordable, so it has to give the matrix that the path
    counting gives on the algebra it stands for.
    """
    checked = 0
    for parameters, edges, orientation, _automorphisms in qr.orientedQuipus(order):
        succ = qr.successors(order, edges, orientation)
        paths = qr.directedPaths(order, succ)
        baseMask, candidates, killMasks, comparable = qr.relationData(order, paths, 2)
        for mask, chosen in qr.relationSets(baseMask, killMasks, comparable):
            algebra = qr.algebraFromMatch({
                'edges': tuple(edges),
                'orientation': tuple(orientation),
                'relations': tuple(candidates[index] for index in chosen),
            })
            assert qr.cartanFromMask(mask, order) == inv.integerCartanMatrix(algebra)
            checked += 1
            if checked > 400:
                return
    assert checked


@pytest.mark.parametrize("order", [6, 7])
def test_the_float_sieve_agrees_with_the_exact_polynomial(order):
    """The batched float determinant is only allowed to be a sieve, not an answer.

    It selects on `det(2 C^T + C)`, which is the polynomial at 2, so it must round
    to exactly that integer -- a sieve that missed would lose matches silently.
    """
    masks = []
    for parameters, edges, orientation, _automorphisms in qr.orientedQuipus(order):
        succ = qr.successors(order, edges, orientation)
        paths = qr.directedPaths(order, succ)
        baseMask, candidates, killMasks, comparable = qr.relationData(order, paths, 2)
        masks.extend(mask for mask, _chosen in qr.relationSets(baseMask, killMasks, comparable))
        if len(masks) > 2000:
            break
    values = qr._fingerprints(masks, order)
    for mask, value in zip(masks, values):
        key = inv.coxeterCoefficients(qr.cartanFromMask(mask, order))
        assert value == sum(c * qr.FINGERPRINT_POINT**d for d, c in enumerate(key))


@pytest.mark.parametrize("order", [8, 9, 10, 11])
def test_the_batched_modular_values_are_the_polynomial(order):
    """The batched check must agree with the one-at-a-time exact polynomial.

    `exactValues` holds the polynomial as residues modulo three primes at
    `order + 1` points, computed over a whole batch at once.  An earlier version
    did the same thing with fraction-free elimination in int64 and silently
    overflowed at order 11 -- every match was lost, with nothing to say so -- so
    this compares it against `invariants.coxeterCoefficients`, which carries
    Python integers and cannot.
    """
    masks = []
    for _parameters, edges, orientation, _automorphisms in qr.orientedQuipus(order):
        succ = qr.successors(order, edges, orientation)
        paths = qr.directedPaths(order, succ)
        baseMask, _candidates, killMasks, comparable = qr.relationData(order, paths, 3)
        masks.extend(mask for mask, _chosen in qr.relationSets(baseMask, killMasks, comparable))
        if len(masks) > 300:
            break
    masks = masks[:300]
    values = qr.exactValues(masks, order)
    for mask, row in zip(masks, values):
        key = inv.coxeterCoefficients(qr.cartanFromMask(mask, order))
        expected = [sum(c * point**degree for degree, c in enumerate(key)) % prime
                    for point in range(order + 1) for prime in qr.MODULI]
        assert list(row) == expected


def test_the_moduli_are_enough_for_the_orders_they_are_used_at():
    """Equal residues mean equal values only while the values are small enough."""
    product = qr._moduliProduct()
    assert 2 * qr.modulusBound(11) < product
    assert 2 * qr.modulusBound(12) > product      # so order 12 must not use them
