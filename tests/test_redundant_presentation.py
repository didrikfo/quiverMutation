"""Seed 1614 of round 118 (research E-310, H-030): a redundant length-3 relation is harmless; a coefficient 32002 is not -1.

Algebra: 5 vertices (1-based here), arrows 1>2 1>3 1>5 2>4 2>5 3>4 4>5, relation 1>2>4 = 1>3>4 and its
consequence 1>2>4>5 = 1>3>4>5.  The round-118 harness encoded -1 as 32002 (= 32003 - 1); the repo works over the
rationals, where x + 32002 y and y + 32002 x are independent and force x = y = 0, so the second relation was no
longer a consequence of the first and the repo "failed" on a different algebra.  With coefficient -1 the gate
agrees with Hom-verified tilting (E-310 engine, 8 of 8 vertex/side cases) and the presentation does not matter.
"""

import networkx as nx
import pytest

from quivermutation import arrowPaths as ap
from quivermutation import procedure
from helpers import quiet

ARR = [(1, 2), (1, 3), (1, 5), (2, 4), (2, 5), (3, 4), (4, 5)]
# tilting (Hom-verified, workshop/rounds/118/skeptic_conv3.py engine) per vertex 1..5, from the E-310 run
TILT_F = [True, True, True, True, False]
TILT_R = [False, True, True, True, True]


def _quiver():
    q = nx.MultiDiGraph()
    q.add_nodes_from(range(1, 6))
    for a, b in ARR:
        q.add_edge(a, b)
    return q


def _rels(c, redundant):
    a = lambda t, h: (t, h, 0)
    r1 = {(a(1, 3), a(3, 4)): 1, (a(1, 2), a(2, 4)): c}
    r2 = {(a(1, 2), a(2, 4), a(4, 5)): 1, (a(1, 3), a(3, 4), a(4, 5)): c}
    return [r1, r2] if redundant else [r1]


def _gates(rels):
    q = _quiver()
    qo, ro = procedure._opposite(q, rels)
    return ([procedure.isMutable(q, rels, v) for v in range(1, 6)],
            [procedure.isMutable(qo, ro, v) for v in range(1, 6)])


@pytest.mark.parametrize("redundant", [False, True])
def test_gate_matches_tilting_with_coefficient_minus_one(redundant):
    assert _gates(_rels(-1, redundant)) == (TILT_F, TILT_R)


def test_redundant_relation_is_in_the_ideal_of_the_first():
    q = _quiver()
    r1, r2 = _rels(-1, True)
    assert ap.isInIdeal(q, [r1], ap.combination([((1, 3, 0), (3, 4, 0), (4, 5, 0))])) is False
    assert ap.isInIdeal(q, [r1], r2)


def test_coefficient_32002_is_not_minus_one():
    """The trap: a residue mod 32003 passed as a coefficient is an integer, and the algebra changes."""
    gates = _gates(_rels(32002, True))
    assert gates != (TILT_F, TILT_R)
    assert not ap.isInIdeal(_quiver(), _rels(32002, False), _rels(32002, True)[1])


def test_child_is_cartan_congruent_with_redundant_relation():
    q, rels = _quiver(), _rels(-1, True)
    checked = 0
    for v in (1, 2, 3, 4):
        nq, nr = quiet(procedure.mutateAtVertex, q, rels, v)
        nq, nr = quiet(procedure.reduce, nq, nr)
        if ap.hasOrientedCycle(nq):         # Cartan matrix of a cyclic child is not defined here
            continue
        checked += 1
        assert all(x == 0 for row in procedure.cartanDiscrepancy(q, rels, v, nq, nr) for x in row)
    assert checked >= 1
