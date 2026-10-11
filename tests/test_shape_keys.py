"""Keys for a quiver with relations that do not see its labels.

The atlas (spec docs/superpowers/specs/2026-09-24-shape-atlas-design.md) counts
shapes across walks out of *different* LNAs, where the same quiver turns up
under different vertex labels.  Everything here is about the one property that
makes that count mean anything: a key changes when the shape changes and not
when the labels do.
"""

import json
import random
from fractions import Fraction

import networkx as nx
import pytest

from quivermutation import fingerprint
from quivermutation import mutation
from quivermutation import nakayama as nk
from quivermutation import procedure
from quivermutation import shapeKeys as sk


def _square():
    """F-027's example: A_10 with one relation on 3..8, right mutation at 3.

    The result has a zero relation 2-4-3 and a commutativity square between
    4-3-8 (two arrows) and 4-5-6-7-8 (four).
    """
    return mutation.quiverMutationAtVertices(nk.LinearNakayamaAlgebra(10, "00500000"), [3])


def _quiver(arrows):
    quiver = nx.MultiDiGraph()
    for tail, head, key in arrows:
        quiver.add_edge(tail, head, key = key)
    return quiver


def _permuted(pathAlg, seed):
    vertices = sorted(pathAlg.quiver.nodes)
    images = vertices[:]
    random.Random(seed).shuffle(images)
    return sk.relabel(pathAlg, dict(zip(vertices, images)))


@pytest.mark.parametrize("seed", range(5))
def test_relabelling_leaves_every_bucket_unchanged(seed):
    algebra = _square()
    moved = _permuted(algebra, seed)
    for level in sk.LEVELS:
        assert sk.bucketOf(algebra, level) == sk.bucketOf(moved, level)


def test_relabel_is_a_relabelling_and_not_a_different_algebra():
    algebra = _square()
    identity = {vertex: vertex for vertex in algebra.quiver.nodes}
    assert fingerprint.canonicalKey(sk.relabel(algebra, identity)) == fingerprint.canonicalKey(algebra)
    moved = _permuted(algebra, 1)
    assert sorted(moved.quiver.nodes) == sorted(algebra.quiver.nodes)
    assert moved.quiver.number_of_edges() == algebra.quiver.number_of_edges()


def test_zero_and_commutativity_agree_at_L1_and_differ_at_L2():
    arrows = [(1, 2, 0), (2, 4, 0), (1, 3, 0), (3, 4, 0)]
    comm = procedure.toPathAlgebra(_quiver(arrows), [
        {((1, 2, 0), (2, 4, 0)): 1, ((1, 3, 0), (3, 4, 0)): -1}])
    zero = procedure.toPathAlgebra(_quiver(arrows), [
        {((1, 2, 0), (2, 4, 0)): 1}, {((1, 3, 0), (3, 4, 0)): 1}])
    assert sk.bucketOf(comm, 1) == sk.bucketOf(zero, 1)
    assert sk.bucketOf(comm, 2) != sk.bucketOf(zero, 2)


def test_serialise_round_trips_through_json():
    algebra = _square()
    back = sk.deserialise(json.loads(json.dumps(sk.serialise(algebra))))
    assert fingerprint.canonicalKey(back) == fingerprint.canonicalKey(algebra)
    assert sk.labelId(back) == sk.labelId(algebra)


def test_deserialise_keeps_integral_coefficients_integral():
    """`fingerprint.digest` hashes `repr`, and `repr(Fraction(1)) != repr(1)`."""
    back = sk.deserialise(json.loads(json.dumps(sk.serialise(_square()))))
    for relation in procedure.relationsFrom(back):
        for coefficient in relation.values():
            assert not (isinstance(coefficient, Fraction) and coefficient.denominator == 1)


def test_features_of_a_line():
    found = sk.features(nk.LinearNakayamaAlgebra(6, "3000"))
    assert found['isLine'] and found['isTree'] and found['isQuipu']
    assert found['cycleRank'] == 0
    assert found['relations'] == 1
    assert found['square'] is None


def test_features_of_the_square():
    found = sk.features(_square())
    assert not found['isLine'] and not found['isTree'] and not found['isQuipu']
    assert found['cycleRank'] == 1
    assert found['square'] == '2x4'
    assert found['defect'] is None


def test_features_of_a_quipu_count_its_defect():
    tree = procedure.toPathAlgebra(_quiver([(1, 2, 0), (2, 3, 0), (2, 4, 0)]), [])
    found = sk.features(tree)
    assert found['isQuipu'] and not found['isLine']
    assert found['defect'] == -1          # no relations, one cord
    assert found['quipu'].startswith('P^(')


def test_describe_names_arrows_and_relations():
    text = sk.describe(_square())
    assert '4->3' in text and 'comm' in text and 'zero' in text


@pytest.mark.parametrize("seed", range(5))
def test_relabelling_leaves_every_key_unchanged(seed):
    algebra = _square()
    moved = _permuted(algebra, seed)
    index = sk.ShapeIndex()
    for level in sk.LEVELS:
        assert index.keyOf(algebra, level) == index.keyOf(moved, level)


def test_a_WL_collision_still_gets_two_keys():
    """Two triangles and a hexagon: WL cannot separate them, VF2 must."""
    triangles = procedure.toPathAlgebra(_quiver([
        (1, 2, 0), (2, 3, 0), (1, 3, 0), (4, 5, 0), (5, 6, 0), (4, 6, 0)]), [])
    hexagon = procedure.toPathAlgebra(_quiver([
        (1, 2, 0), (2, 3, 0), (3, 4, 0), (4, 5, 0), (5, 6, 0), (1, 6, 0)]), [])
    assert sk.bucketOf(triangles, 0) == sk.bucketOf(hexagon, 0)
    index = sk.ShapeIndex()
    assert index.keyOf(triangles, 0) != index.keyOf(hexagon, 0)


def test_zero_and_commutativity_differ_as_keys_at_L2_and_L3():
    arrows = [(1, 2, 0), (2, 4, 0), (1, 3, 0), (3, 4, 0)]
    comm = procedure.toPathAlgebra(_quiver(arrows), [
        {((1, 2, 0), (2, 4, 0)): 1, ((1, 3, 0), (3, 4, 0)): -1}])
    zero = procedure.toPathAlgebra(_quiver(arrows), [
        {((1, 2, 0), (2, 4, 0)): 1}, {((1, 3, 0), (3, 4, 0)): 1}])
    index = sk.ShapeIndex()
    assert index.keyOf(comm, 1) == index.keyOf(zero, 1)
    assert index.keyOf(comm, 2) != index.keyOf(zero, 2)
    assert index.keyOf(comm, 3) != index.keyOf(zero, 3)


def _kronecker(relations):
    return procedure.toPathAlgebra(_quiver([(1, 2, 0), (1, 2, 1), (2, 3, 0)]), relations)


def test_parallel_naming_and_sign_gauge_agree_at_L3():
    first = _kronecker([{((1, 2, 0), (2, 3, 0)): 1}])
    second = _kronecker([{((1, 2, 1), (2, 3, 0)): 1}])
    plus = _kronecker([{((1, 2, 0), (2, 3, 0)): 1, ((1, 2, 1), (2, 3, 0)): 1}])
    minus = _kronecker([{((1, 2, 0), (2, 3, 0)): 1, ((1, 2, 1), (2, 3, 0)): -1}])
    index = sk.ShapeIndex()
    assert index.keyOf(first, 3) == index.keyOf(second, 3)
    assert index.keyOf(plus, 3) == index.keyOf(minus, 3)
    assert index.keyOf(first, 3) != index.keyOf(plus, 3)


def test_L3_splits_what_the_L2_skeleton_joins():
    """Two bundles 1 => 2 => 3 and two zero relations of length two on each.

    The skeletons are the same (two `zero:2` relations, each 1 -> 3 through 2),
    so L2 is one shape.  The algebras are not: in {ac, bd} each arrow out of 1
    kills one arrow out of 2, in {ac, ad} one arrow out of 1 kills both -- as
    tensors, a span of two rank-one tensors against a whole a (x) W, which no
    change of basis carries onto each other.  L3 must say two.  E-054 found L2
    and L3 equal at every length it counted; this is what keeps the L3 column
    from being a copy of L2 by construction.
    """
    def bundles(relations):
        return procedure.toPathAlgebra(
            _quiver([(1, 2, 0), (1, 2, 1), (2, 3, 0), (2, 3, 1)]),
            [{path: 1} for path in relations])

    ac, ad, bd = (((1, 2, 0), (2, 3, 0)), ((1, 2, 0), (2, 3, 1)), ((1, 2, 1), (2, 3, 1)))
    split, shared = bundles([ac, bd]), bundles([ac, ad])
    index = sk.ShapeIndex()
    assert index.keyOf(split, 2) == index.keyOf(shared, 2)
    assert index.keyOf(split, 3) != index.keyOf(shared, 3)
    assert index.capHits == 0 and index.uncanonical == 0


def test_two_different_lnas_are_two_keys_at_L3_and_one_at_L0():
    first = nk.LinearNakayamaAlgebra(6, "3000")
    second = nk.LinearNakayamaAlgebra(6, "0300")
    index = sk.ShapeIndex()
    assert index.keyOf(first, 0) == index.keyOf(second, 0)
    assert index.keyOf(first, 3) != index.keyOf(second, 3)


def test_a_key_carries_its_level():
    assert sk.ShapeIndex().keyOf(_square(), 2).startswith('L2:')


def test_the_index_keeps_representatives_as_strings():
    """Task 8a: a built graph per representative ran n = 9 out of memory."""
    index = sk.ShapeIndex()
    for level in sk.LEVELS:
        index.keyOf(_square(), level)
    for level in (0, 1, 2):
        assert all(isinstance(r, str) for reps in index._representatives[level].values()
                   for r in reps)
    assert all(isinstance(text, str) for reps in index._representatives[3].values()
               for text, _canonical in reps)
