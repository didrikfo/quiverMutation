"""A cyclic quiver is outside the setting where the Cartan matrix, ideal membership and
minimal generators are proved correct (paths are simple in the vertices, research E-193):
they warn, and stay silent on an acyclic quiver.  Round 067, thread T10/T5."""
import warnings

import networkx as nx
import pytest

from quivermutation import arrowPaths as ap
from quivermutation import procedure


def quiver(arrows):
    q = nx.MultiDiGraph()
    for a, b in arrows:
        q.add_edge(a, b)
    return q


def test_hasOrientedCycle():
    assert not ap.hasOrientedCycle(quiver([(1, 2), (2, 3), (1, 3), (1, 3)]))
    assert ap.hasOrientedCycle(quiver([(1, 2), (2, 1)]))
    assert ap.hasOrientedCycle(quiver([(1, 1), (1, 2)]))


def test_acyclic_is_silent():
    q = quiver([(1, 2), (2, 3)])
    rel = [{((1, 2, 0), (2, 3, 0)): 1}]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ap.cartanMatrix(q, rel)
        ap.isInIdeal(q, rel, rel[0])
        procedure._minimalGenerators(q, rel)


def test_cyclic_warns_everywhere_and_keeps_the_answer():
    q = quiver([(1, 2), (2, 3), (3, 1)])
    rel = [{((1, 2, 0), (2, 3, 0)): 1}]
    with pytest.warns(ap.CyclicQuiverWarning, match="oriented cycle"):
        m = ap.cartanMatrix(q, rel)
    with pytest.warns(ap.CyclicQuiverWarning):
        assert ap.isInIdeal(q, rel, rel[0])
    with pytest.warns(ap.CyclicQuiverWarning):
        assert procedure._minimalGenerators(q, rel) == rel
    with warnings.catch_warnings():  # the warning does not change the value
        warnings.simplefilter("ignore")
        assert ap.cartanMatrix(q, rel) == m
