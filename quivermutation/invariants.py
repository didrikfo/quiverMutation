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


def generateAllCoxeterPolynomials(length):

    # generate Kupisch series
    kupisch = [[[1]]]
    for i in range(1, length):
        kupisch.append([])
        for s in kupisch[i - 1]:
            for x in range(2, s[0] + 2):
                kupisch[i].append([x] + s)
    print('Number of different possible sets of relations: ',len(kupisch[length - 1]))

    polynomials = []
    for s in kupisch[length - 1]:
        M = sympy.Matrix([[1 if (m >= n and m < n + s[n]) else 0 for m in range(length)] for n in range(length)])
        C = - M * M.inv().transpose()
        p = sympy.factor(C.charpoly(sympy.Symbol("x")).as_expr())
        if p not in polynomials:
            polynomials.append(p)
    print('Coxeter polynomials: ',polynomials)
    print('Number of different Coxeter polynomials: ' ,len(polynomials))
    return polynomials


def coxPolyOfTree(tree):
    adjMat = sympy.Matrix(nx.adjacency_matrix(tree).todense(), dtype=int)
    triuMat = np.triu(Matrix(adjMat))
    mat = sympy.Matrix(triuMat, dtype=int) + sympy.eye(len(tree))#sympy.Matrix(np.triu(Matrix(matrixAsList))) + sympy.eye(n)
    print(np.matrix(mat))
    matInvTrans = mat.inv().transpose()
    coxeterMatrix = -matInvTrans * mat
    coxeterPolynomial = coxeterMatrix.charpoly()
    #print(coxeterPolynomial.as_expr())
    return coxeterPolynomial


def cartanMatrixForCanonicalAlgebra(pathAlg):
    quiv = pathAlg.quiver
    vertices = quiv.nodes
    cartanMatrix = eye(len(vertices), len(vertices))
    for i in vertices:
        for j in vertices:
            if i == j:
                cartanMatrix[j-1,i-1] = 1
            else:
                allPaths = list(nx.all_simple_paths(pathAlg.quiver, i, j))
                cartanMatrix[j-1,i-1] = min(len(allPaths), 2)
    return cartanMatrix


def coxeterPolyForCanonicalAlgebra(pathAlg):
    cartanMat = cartanMatrixForCanonicalAlgebra(pathAlg)
    print(np.matrix(cartanMat))
    cartanMatInvTrans = cartanMat.inv().transpose()
    coxeterMatrix = -cartanMatInvTrans*cartanMat
    coxeterPolynomial = coxeterMatrix.charpoly()
    return coxeterPolynomial


def divisors(n):
    # get factors and their counts
    factors = {}
    nn = n
    i = 2
    while i*i <= nn:
        while nn % i == 0:
            factors[i] = factors.get(i, 0) + 1
            nn //= i
        i += 1
    if nn > 1:
        factors[nn] = factors.get(nn, 0) + 1
    primes = list(factors.keys())
    # generates factors from primes[k:] subset
    def generate(k):
        if k == len(primes):
            yield 1
        else:
            rest = generate(k+1)
            prime = primes[k]
            for factor in rest:
                prime_to_i = 1
                # prime_to_i iterates prime**i values, i being all possible exponents
                for _ in range(factors[prime] + 1):
                    yield factor * prime_to_i
                    prime_to_i *= prime
    # in python3, `yield from generate(0)` would also work
    for factor in generate(0):
        yield factor
