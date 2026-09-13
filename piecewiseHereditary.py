"""Cheap certificates that a linear Nakayama algebra is *not* piecewise hereditary.

An algebra is piecewise hereditary if its derived category is equivalent to that
of a hereditary abelian category.  By Happel's classification such a category is
either the module category of a hereditary algebra or derived equivalent to a
canonical algebra, so an LNA that is not piecewise hereditary is in **no** quipu
class -- it cannot be derived equivalent to the path algebra of any tree.

That makes these criteria useful to the classification directly: they say, without
any searching, that a class is not one of the quipu classes, which both saves
looking for a quipu and separates the class from every quipu class outright.

Both criteria are from

    K. M. Jacobsen and D. Fosse (and co-authors), "Non-piecewise hereditary
    Nakayama algebras", arXiv:2310.08346, section "Families of non-piecewise
    hereditary algebras".

A relation is written (start vertex, number of arrows); its end vertex is
start + arrows.  "Length" throughout means number of arrows, so the shortest
admissible relation has length 2.
"""


def relations(relLengths):
    """(start vertex, arrows) for each relation of an LNA, left to right."""
    return [(index + 1, arrows) for index, arrows in enumerate(relLengths) if arrows]


def _end(relation):
    return relation[0] + relation[1]


def overlapInArrows(first, second):
    """How many arrows two relations share.

    A relation (s, a) covers the arrows s, s+1, ..., s+a-1, so the shared count
    is min(end) - max(start), floored at zero.
    """
    return max(0, min(_end(first), _end(second)) - max(first[0], second[0]))


def failsBigOverlapCriterion(length, relLengths):
    """Proposition A13: a pair of relations overlapping by at least six arrows.

    Returns the witnessing pair, or None.

    The conditions, with `alpha` the leftmost of the two:

    * they overlap by at least six arrows (they need not be consecutive);
    * at least two vertices lie between the two starts, and at least two between
      the two ends;
    * no relation runs from the first or second vertex before the start of `beta`
      to the first or second vertex after the end of `alpha`.

    The smallest algebra it applies to is A_13 with relations 1 -> 10 and
    4 -> 13.
    """
    allRelations = relations(relLengths)
    for index, alpha in enumerate(allRelations):
        for beta in allRelations[index + 1:]:
            if overlapInArrows(alpha, beta) < 6:
                continue
            if beta[0] - alpha[0] < 3:
                continue
            if _end(beta) - _end(alpha) < 3:
                continue
            blocked = any(
                other[0] in (beta[0] - 1, beta[0] - 2)
                and _end(other) in (_end(alpha) + 1, _end(alpha) + 2)
                for other in allRelations
            )
            if not blocked:
                return alpha, beta
    return None


def failsFlankedPairCriterion(length, relLengths):
    """Proposition A9: an overlapping pair with a long relation on each side.

    Returns the witnessing pair, or None.

    The conditions, with `alpha` the leftmost of the two and at least 9 vertices
    in the quiver:

    1. no relation starts at the vertex directly before `alpha`'s start and ends
       at or before the vertex directly after `beta`'s start;
    2. no relation ends at the vertex directly after `beta`'s end and starts at
       or after the vertex directly before `alpha`'s end;
    3. some relation of length at least 3 lies entirely before `beta`, possibly
       sharing a vertex with it;
    4. some relation of length at least 3 lies entirely after `alpha`, possibly
       sharing a vertex with it.

    The archetype is the A_9 with relations 1 -> 4, 3 -> 6, 4 -> 7 and 6 -> 9,
    which the paper shows is the only LNA of length 9 that is not piecewise
    hereditary.
    """
    if length < 9:
        return None
    allRelations = relations(relLengths)
    for index, alpha in enumerate(allRelations):
        for beta in allRelations[index + 1:]:
            if overlapInArrows(alpha, beta) < 2:
                continue
            if any(other[0] == alpha[0] - 1 and _end(other) <= beta[0] + 1
                   for other in allRelations):
                continue
            if any(_end(other) == _end(beta) + 1 and other[0] >= _end(alpha) - 1
                   for other in allRelations):
                continue
            before = any(_end(other) <= beta[0] and other[1] >= 3 for other in allRelations)
            after = any(other[0] >= _end(alpha) and other[1] >= 3 for other in allRelations)
            if before and after:
                return alpha, beta
    return None


CRITERIA = (
    ("flanked overlapping pair (Proposition A9)", failsFlankedPairCriterion),
    ("overlap of six or more arrows (Proposition A13)", failsBigOverlapCriterion),
)


def isNotPiecewiseHereditary(length, relLengths):
    """Whether some criterion certifies the algebra is not piecewise hereditary.

    A False means only that these criteria say nothing, not that the algebra is
    piecewise hereditary -- they are sufficient conditions, not a
    characterisation.
    """
    return certificate(length, relLengths) is not None


def certificate(length, relLengths):
    """(criterion name, witnessing pair of relations), or None."""
    for name, criterion in CRITERIA:
        witness = criterion(length, relLengths)
        if witness is not None:
            return name, witness
    return None


# ---------------------------------------------------------------------------
# Canonical algebras
#
# Happel's classification says a hereditary abelian category is either the
# module category of a hereditary algebra or derived equivalent to a canonical
# algebra.  So a piecewise hereditary LNA that is *not* derived equivalent to a
# quipu algebra should be of canonical type, and the Coxeter polynomial names
# which one.
# ---------------------------------------------------------------------------


def canonicalCoxeterPolynomial(weights, symbol = None):
    """The Coxeter polynomial of the canonical algebra of a weight type.

    For weights (p_1, ..., p_t) it is

        (x - 1)^2 * prod_i (1 + x + ... + x^(p_i - 1)),

    of degree 2 + sum(p_i - 1), which is the number of vertices.
    """
    import sympy

    symbol = sympy.Symbol("lambda") if symbol is None else symbol
    product = (symbol - 1) ** 2
    for weight in weights:
        product *= sum(symbol ** power for power in range(weight))
    return sympy.expand(product)


def canonicalVertexCount(weights):
    return 2 + sum(weight - 1 for weight in weights)


def weightTypesOfOrder(order, maxParts = 6):
    """Every weight type (p_1 <= ... <= p_t), each p_i >= 2, on `order` vertices."""
    results = []

    def extend(current, smallest):
        if canonicalVertexCount(current) == order and len(current) >= 1:
            results.append(tuple(current))
        for weight in range(smallest, order + 1):
            candidate = current + [weight]
            if canonicalVertexCount(candidate) > order or len(candidate) > maxParts:
                break
            extend(candidate, weight)

    extend([], 2)
    return sorted(set(results))


def canonicalWeightType(length, coxeterPolynomial, maxParts = 6):
    """The weight type whose canonical algebra has this Coxeter polynomial, or None.

    A positive identification, unlike the non-piecewise-hereditary criteria: if a
    class' Coxeter polynomial is that of the canonical algebra of weight type
    `w`, and the class is not a quipu class, then `w` names it, and two classes
    with different weight types are certainly different classes.

    Note this is still evidence rather than proof, since the Coxeter polynomial
    does not characterise a derived equivalence class -- but it is the same
    quality of evidence as the rest of the pipeline's Coxeter reasoning, and it
    is the argument Happel and Seidel used.
    """
    import sympy

    variable = sympy.Symbol("x")
    target = sympy.expand(parseCoxeterPolynomial(coxeterPolynomial, variable))
    for weights in weightTypesOfOrder(length, maxParts):
        if canonicalCoxeterPolynomial(weights, variable) == target:
            return weights
    return None


def parseCoxeterPolynomial(polynomial, variable = None):
    """A Coxeter polynomial as a sympy expression in `variable`.

    Accepts a sympy expression or the printed string the table stores.  Those
    strings name the variable `lambda`, which Python's tokenizer will not accept
    as an identifier, so it is renamed before parsing rather than fed to sympify
    as it stands.
    """
    import sympy

    variable = sympy.Symbol("x") if variable is None else variable
    if isinstance(polynomial, str):
        return sympy.sympify(polynomial.replace("lambda", "x"), locals = {"x": variable})
    symbols = list(polynomial.free_symbols)
    if not symbols:
        return polynomial
    return polynomial.subs(symbols[0], variable)


TUBULAR_TYPES = ((2, 2, 2, 2), (3, 3, 3), (2, 4, 4), (2, 3, 6))


def isTubular(weights):
    """Whether a weight type is one of the four tubular ones.

    A canonical algebra of tubular type is the derived-equivalence boundary
    between the domestic and the wild cases, and its Coxeter polynomial has
    (x - 1)^2 as a factor.
    """
    return tuple(sorted(weights)) in TUBULAR_TYPES
