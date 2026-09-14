"""Derived invariants: the Cartan matrix and the Coxeter polynomial.

The Coxeter polynomial is what the classification compares, being invariant
along any legal mutation path.  It is not a complete invariant, and research
F-010 says exactly where it fails: at cospectral quipus, the first pair of which
is at order 9.
"""

import networkx as nx
import numpy as np
import sympy
from sympy.matrices import Matrix, eye

from . import paths
from . import relationAlgebra



def cartanMatrix(pathAlg, exact = True):
    """The Cartan matrix: entry (j, i) is dim e_j (kQ/I) e_i.

    With exact=True the dimensions come from relationAlgebra, which decides
    which combinations of paths are zero by linear algebra over the ideal.  With
    exact=False they come from numberOfPathsUpToRels, which counts paths up to a
    partial closure under the commutativity relations and calls a path zero when
    a zero relation sits contiguously inside it.

    The two agree on every LNA of length <= 8 and on every quiver reached by
    walking mutations of depth <= 3 out of the LNAs of length 5 to 7, so this
    changes no published number.  They do not agree in general: see the 2x2
    commutative grid in tests/test_relation_algebra.py, where a zero relation on
    one path kills all three and only the exact version notices.

    The exact version costs 2 to 4 times as much on LNAs, which is nothing at
    the rate the pipeline calls it -- once per class, not once per mutation.
    """
    if exact:
        return relationAlgebra.cartanMatrixExact(pathAlg)
    quiv = pathAlg.quiver
    vertices = quiv.nodes
    cartanMatrix = eye(len(vertices), len(vertices))
    for i in vertices:
        for j in vertices:
            if i == j:
                cartanMatrix[j-1,i-1] = paths.numberOfPathsUpToRels(pathAlg, i, j) + 1
            else:
                cartanMatrix[j-1,i-1] = paths.numberOfPathsUpToRels(pathAlg, i, j)
    return cartanMatrix


def coxeterPoly(pathAlg, exact = True):
    """The Coxeter polynomial, the derived invariant the classification uses."""
    cartanMat = cartanMatrix(pathAlg, exact)
    cartanMatInvTrans = cartanMat.inv().transpose()
    coxeterMatrix = -cartanMatInvTrans*cartanMat
    coxeterPolynomial = coxeterMatrix.charpoly()
    return coxeterPolynomial
