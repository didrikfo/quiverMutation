"""The gate (`mutationIsPossibleAtVertex`) is only a necessary condition (research E-080, E-068).

The 5-vertex algebra a>b>d, a>c>d, d>e with the single relation abde = acde has a 2-dimensional
e_a A e_d whose image under p |-> p(d>e) is 1-dimensional, so the two-term complex at d is not
tilting (Aihara-Iyama 2.32(b) = Ladkani 2.3(c)).  The gate nevertheless admits d.  This pins that
behaviour; it does not add `isTilting` to the library.  A true commutative square (abd = acd)
is the control: there the map is injective.
"""

import quivermutation as qm
from quivermutation import arrowPaths
from quivermutation import procedure
from helpers import path_algebra, quiet

ARROWS = [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]


def _dims(pa):
    """dim of the space of paths i ~> j modulo I, for the pairs of interest (1, 4) and (1, 5)."""
    rels = procedure.relationsFrom(pa)
    C = arrowPaths.cartanMatrix(pa.quiver, rels, exact=True)
    return max(C[0][3], C[3][0]), max(C[0][4], C[4][0])


def test_long_commutativity_relation_gate_admits_but_map_is_not_injective():
    pa = path_algebra(ARROWS, rels=[[(1, 2, 4, 5), (1, 3, 4, 5)]])
    assert quiet(qm.mutationIsPossibleAtVertex, pa, 4)
    dim_ad, dim_ae = _dims(pa)
    assert (dim_ad, dim_ae) == (2, 1)      # d>e is the only arrow out of d: 2 -> 1 cannot be injective


def test_true_commutative_square_map_is_injective():
    pa = path_algebra(ARROWS, rels=[[(1, 2, 4), (1, 3, 4)]])
    assert quiet(qm.mutationIsPossibleAtVertex, pa, 4)
    dim_ad, dim_ae = _dims(pa)
    assert (dim_ad, dim_ae) == (1, 1)
