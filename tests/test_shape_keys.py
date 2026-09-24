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
