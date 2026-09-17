"""Derived invariants: the Cartan matrix and the Coxeter polynomial.

The Coxeter polynomial is what the classification compares, being invariant
along any legal mutation path.  It is not a complete invariant, and research
F-010 says exactly where it fails: at cospectral quipus, the first pair of which
is at order 9.
"""

from sympy.matrices import eye

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
    vertices = pathAlg.quiver.nodes
    matrix = eye(len(vertices), len(vertices))
    for i in vertices:
        for j in vertices:
            counted = paths.numberOfPathsUpToRels(pathAlg, i, j)
            # The identity path at a vertex is in the algebra but not in rels.
            matrix[j - 1, i - 1] = counted + 1 if i == j else counted
    return matrix


def coxeterPoly(pathAlg, exact = True):
    """The Coxeter polynomial, the derived invariant the classification uses."""
    cartanMat = cartanMatrix(pathAlg, exact)
    cartanMatInvTrans = cartanMat.inv().transpose()
    coxeterMatrix = -cartanMatInvTrans*cartanMat
    coxeterPolynomial = coxeterMatrix.charpoly()
    return coxeterPolynomial


# -- the polynomial as integer coefficients ------------------------------
#
# `coxeterPoly` is the readable route: build the Cartan matrix, invert it over
# the rationals, and hand sympy a symbolic charpoly.  It costs milliseconds,
# which is nothing when the caller is a classification asking once per class and
# far too much when the caller is a search comparing millions of algebras against
# a table of polynomials.  The functions below are that second route, and they
# rest on one identity.
#
# For a quiver with no oriented cycles, ordering the vertices along a topological
# order makes the Cartan matrix unitriangular, so det C = 1.  Then
#
#     det(lambda I - Phi) = det(lambda I + C^-T C)
#                         = det(C^-T) det(lambda C^T + C)
#                         = det(lambda C^T + C),
#
# so the Coxeter polynomial is the determinant of a matrix whose entries are
# integers and lambda, with no inversion anywhere.  Evaluating it at n + 1
# integer points and interpolating gives the coefficients in exact integer
# arithmetic, which is both faster than the symbolic route and hashable, so a
# polynomial can be a dictionary key rather than something to compare with
# `simplify`.


def integerCartanMatrix(pathAlg):
    """The Cartan matrix as a list of lists of Python ints, vertices sorted.

    Row j, column i is dim e_j (kQ/I) e_i, as in `cartanMatrix`, computed the
    inexact way -- by counting paths that no relation kills.  For a *monomial*
    ideal on a quiver with no oriented cycles the two agree, since a path is zero
    exactly when it contains a generator, and every algebra this is used on --
    trees, quipus with zero relations, LNAs -- is of that kind.
    """
    vertices = sorted(pathAlg.quiver.nodes)
    index = {vertex: position for position, vertex in enumerate(vertices)}
    size = len(vertices)
    matrix = [[0] * size for _ in range(size)]
    for source in vertices:
        for target in vertices:
            counted = paths.numberOfPathsUpToRels(pathAlg, source, target)
            if source == target:
                counted += 1
            matrix[index[target]][index[source]] = counted
    return matrix


def coxeterCoefficients(cartan):
    """The Coxeter polynomial of a unimodular Cartan matrix, low degree first.

    `cartan` is a square matrix of ints as nested lists.  The result is a tuple
    of n + 1 ints, entry d being the coefficient of lambda^d, so it is exact,
    hashable, and usable as a dictionary key.  Raises when det C is not 1, which
    is the condition the identity above needs and which every algebra of finite
    global dimension on an acyclic quiver satisfies.
    """
    size = len(cartan)
    transpose = [[cartan[row][column] for row in range(size)] for column in range(size)]
    if _integerDeterminant(transpose) != 1:
        raise ValueError("the Cartan matrix is not unimodular, so the "
                         "determinant identity does not apply")
    values = [
        _integerDeterminant([
            [point * transpose[row][column] + cartan[row][column] for column in range(size)]
            for row in range(size)
        ])
        for point in range(size + 1)
    ]
    return _interpolate(values)


def coxeterKey(pathAlg):
    """`coxeterCoefficients` of a path algebra's Cartan matrix."""
    return coxeterCoefficients(integerCartanMatrix(pathAlg))


def polynomialFromCoefficients(coefficients, symbol = None):
    """A sympy expression for a coefficient tuple, for printing and comparing."""
    import sympy
    variable = sympy.Symbol("lambda") if symbol is None else symbol
    return sum(coefficient * variable**degree
               for degree, coefficient in enumerate(coefficients))


def formatCoefficients(coefficients):
    """The polynomial as a one-line string, highest degree first."""
    pieces = []
    for degree in range(len(coefficients) - 1, -1, -1):
        coefficient = coefficients[degree]
        if coefficient == 0:
            continue
        if degree == 0:
            term = str(abs(coefficient))
        else:
            power = "x" if degree == 1 else "x^{0}".format(degree)
            term = power if abs(coefficient) == 1 else "{0}{1}".format(abs(coefficient), power)
        sign = "-" if coefficient < 0 else ("+" if pieces else "")
        pieces.append((sign + " " + term) if pieces else (sign + term))
    return " ".join(pieces) if pieces else "0"


def _integerDeterminant(matrix):
    """The determinant of an integer matrix, by fraction-free elimination.

    Bareiss: every intermediate entry stays an integer and stays bounded by a
    minor of the original, so there is no rounding to worry about and no
    rationals to carry.  The matrix is copied, not written through.
    """
    size = len(matrix)
    if size == 0:
        return 1
    rows = [row[:] for row in matrix]
    sign = 1
    previousPivot = 1
    for step in range(size - 1):
        if rows[step][step] == 0:
            for candidate in range(step + 1, size):
                if rows[candidate][step] != 0:
                    rows[step], rows[candidate] = rows[candidate], rows[step]
                    sign = -sign
                    break
            else:
                return 0
        pivot = rows[step][step]
        for row in range(step + 1, size):
            rowEntry = rows[row][step]
            target = rows[row]
            pivotRow = rows[step]
            for column in range(step + 1, size):
                target[column] = (pivot * target[column] - rowEntry * pivotRow[column]) // previousPivot
            target[step] = 0
        previousPivot = pivot
    return sign * rows[size - 1][size - 1]


def _interpolate(values):
    """The integer coefficients of the polynomial through (0, v0), (1, v1), ...

    Newton's divided differences over the integers: the points are consecutive,
    so the k-th difference is divisible by k! and the whole thing can be done
    with exact integer division at the end.  Returns low degree first.
    """
    from fractions import Fraction

    size = len(values)
    differences = [Fraction(value) for value in values]
    for step in range(1, size):
        for index in range(size - 1, step - 1, -1):
            differences[index] = (differences[index] - differences[index - 1]) / step
    # Horner backwards through the Newton basis (x - 0)(x - 1)...(x - k + 1).
    coefficients = [Fraction(0)] * size
    for step in range(size - 1, -1, -1):
        shifted = [Fraction(0)] * size
        for degree in range(size - 1):
            shifted[degree + 1] += coefficients[degree]
            shifted[degree] -= step * coefficients[degree]
        shifted[0] += differences[step]
        coefficients = shifted
    result = []
    for coefficient in coefficients:
        if coefficient.denominator != 1:
            raise ValueError("interpolation left a non-integer coefficient")
        result.append(int(coefficient))
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)
