"""`longSquare` of workshop/rounds/026/toolsmith_longsquare.py (E-108): a long square through a doubled arrow.

The test is on the algebra at v with exactly one out-arrow b: some relation with >= 2 paths all ending ..., b and starting at one
vertex, with distinct second-to-last arrows. The round 023 version compared predecessor *vertices*, so two paths through a doubled
arrow (same vertex list) were not distinguished and the square was missed."""
import importlib.util
import pathlib

import networkx as nx

from quivermutation import procedure as pr

_p = pathlib.Path(__file__).resolve().parents[1] / "workshop" / "rounds" / "026" / "toolsmith_longsquare.py"
_spec = importlib.util.spec_from_file_location("toolsmith_longsquare", _p)
ls = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ls)


def quiver(arrows):
    q = nx.MultiDiGraph()
    for t, h in arrows:
        q.add_edge(t, h)
    return q


def test_parallel_arrow_long_square_is_found():
    # 3 -> 2 => 5 -> 4 (2 => 5 doubled); relation: the two paths 3,2,5,4 agree. v = 5 has one out-arrow, 5 -> 4.
    q = quiver([(3, 2), (2, 5), (2, 5), (5, 4)])
    a, g0, g1, e = (3, 2, 0), (2, 5, 0), (2, 5, 1), (5, 4, 0)
    alg = pr.toPathAlgebra(q, [{(a, g0, e): 1, (a, g1, e): -1}])
    assert alg.rels == [[[3, 2, 5, 4], [3, 2, 5, 4]]]
    assert ls.longSquareOld(alg, 5) is False      # the defect: same vertex list twice
    assert ls.longSquare(alg, 5) is True


def test_ordinary_square_still_found_and_non_squares_still_not():
    q = quiver([(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)])
    sq = pr.toPathAlgebra(q, [{((1, 2, 0), (2, 4, 0), (4, 5, 0)): 1, ((1, 3, 0), (3, 4, 0), (4, 5, 0)): -1}])
    assert ls.longSquareOld(sq, 4) and ls.longSquare(sq, 4)
    short = pr.toPathAlgebra(q, [{((1, 2, 0), (2, 4, 0)): 1, ((1, 3, 0), (3, 4, 0)): -1}])   # relation stops at v
    assert not ls.longSquare(short, 4)
    mono = pr.toPathAlgebra(q, [{((1, 2, 0), (2, 4, 0), (4, 5, 0)): 1}])                      # one path
    assert not ls.longSquare(mono, 4)


def test_parallel_arrows_with_same_penultimate_arrow_are_not_a_square():
    # paths differ only before the doubled arrow: second-to-last arrows coincide, so no square at v
    q = quiver([(1, 2), (1, 2), (2, 5), (5, 4)])
    e = (5, 4, 0); g = (2, 5, 0)
    alg = pr.toPathAlgebra(q, [{((1, 2, 0), g, e): 1, ((1, 2, 1), g, e): -1}])
    assert not ls.longSquare(alg, 5)
