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


# ---------------------------------------------------------------------------
# Domestic weight types, and why a `C(...)` name for one is never the answer
#
# A canonical algebra is derived equivalent to a hereditary algebra exactly when
# its weight type is *domestic*: (p, q), (2, 2, n), (2, 3, 3), (2, 3, 4) or
# (2, 3, 5), whose hereditary partner has the corresponding extended Dynkin
# type -- A~, D~, E~6, E~7, E~8.  Everything else is tubular or wild, and there
# the canonical algebra is not derived equivalent to any hereditary algebra.
#
# So a class the pipeline names `C(w)` for a domestic `w` is claiming both that
# the class is of canonical type and that it is the class of a hereditary
# algebra -- and for the tree types that hereditary algebra is a quipu, so the
# class is one of the quipu classes and carries a quipu name already.  A
# domestic `C(...)` is therefore always a merge the search has not found yet,
# never a class of its own.  F-018 is the n = 9 pair where exactly that happened.
# ---------------------------------------------------------------------------

DOMESTIC_TREE_TYPES = {
    (2, 3, 3): "E~6",
    (2, 3, 4): "E~7",
    (2, 3, 5): "E~8",
}


def isDomestic(weights):
    """Whether a canonical algebra of this weight type is of domestic type.

    Those are the ones derived equivalent to a hereditary algebra: two weights,
    or (2, 2, n), or one of the three exceptional triples.
    """
    sortedWeights = tuple(sorted(weights))
    if len(sortedWeights) < 2:
        return False
    if len(sortedWeights) == 2:
        return True
    if len(sortedWeights) > 3:
        return False
    return sortedWeights[:2] == (2, 2) or sortedWeights in DOMESTIC_TREE_TYPES


def affineTypeOfDomesticWeightType(weights):
    """The extended Dynkin type of the hereditary partner, or None.

    `(p, q)` gives A~_(p+q-1), whose underlying graph is a cycle rather than a
    tree, so it is named here but `affineTreeOfDomesticWeightType` returns
    nothing for it.
    """
    sortedWeights = tuple(sorted(weights))
    if not isDomestic(sortedWeights):
        return None
    if len(sortedWeights) == 2:
        return "A~{0}".format(sum(sortedWeights) - 1)
    if sortedWeights in DOMESTIC_TREE_TYPES:
        return DOMESTIC_TREE_TYPES[sortedWeights]
    return "D~{0}".format(sortedWeights[2] + 2)


def affineTreeOfDomesticWeightType(weights):
    """The extended Dynkin diagram as an undirected graph, or None.

    None for a weight type that is not domestic, and for the two-weight types,
    whose diagram A~ is a cycle.  The graph is what identifies the class: two
    hereditary algebras of tree type are derived equivalent exactly when their
    underlying graphs are isomorphic, so this is the tree the class' quipu name
    must encode.
    """
    import networkx as nx

    sortedWeights = tuple(sorted(weights))
    if not isDomestic(sortedWeights) or len(sortedWeights) == 2:
        return None
    graph = nx.Graph()
    if sortedWeights in DOMESTIC_TREE_TYPES:
        # A single branch vertex with three arms: E~6 is 2, 2, 2, E~7 is 1, 3, 3
        # and E~8 is 1, 2, 5.
        arms = {(2, 3, 3): (2, 2, 2), (2, 3, 4): (1, 3, 3),
                (2, 3, 5): (1, 2, 5)}[sortedWeights]
        graph.add_node("centre")
        for index, arm in enumerate(arms):
            previous = "centre"
            for step in range(arm):
                node = (index, step)
                graph.add_edge(previous, node)
                previous = node
        return graph
    # D~(n+2): a path of n - 1 vertices with two leaves hung on each end.  For
    # n = 2 the path is the single vertex of D~4, which carries all four.
    n = sortedWeights[2]
    spine = list(range(n - 1))
    graph.add_nodes_from(spine)
    for left, right in zip(spine, spine[1:]):
        graph.add_edge(left, right)
    for end, side in ((spine[0], "left"), (spine[-1], "right")):
        for which in (0, 1):
            graph.add_edge(end, (side, which))
    return graph


# ---------------------------------------------------------------------------
# Propagating a certificate by removing vertices
#
# Corollary "removevertex" of arXiv:2310.08346: if a Nakayama algebra is
# piecewise hereditary, then so is the algebra obtained by removing any one
# vertex -- merging the arrows through it into a composite, and extending the
# relations that start or end there.
#
# Read contrapositively that propagates a certificate *upward*: if some
# one-vertex deletion of an algebra is not piecewise hereditary, neither is the
# algebra.  This is the useful direction to compute in.  The same corollary read
# forwards ("introducevertex") adds a vertex, but an extension is not unique --
# there are many algebras of length n+1 restricting to a given one of length n --
# so going up means enumerating a branching set of possibilities, while going
# down is one deterministic algebra per vertex.  So the recursion here deletes.
#
# What is preserved, and what is not: the corollary is about piecewise
# heredity only.  Deleting a vertex does *not* preserve the derived equivalence
# class, the Coxeter polynomial, or the quipu -- it is a statement about one
# property, and the certificate it carries is exactly "this class is not a quipu
# class", nothing more.
# ---------------------------------------------------------------------------


def removeVertex(length, relLengths, vertex):
    """The algebra obtained by deleting one vertex, per corollary removevertex.

    Returns (new length, new relation lengths), or None if the result is not an
    admissible LNA.

    Arrow j runs from vertex j to vertex j+1, so deleting an interior vertex i
    merges arrows i-1 and i into one composite.  A relation spanning both of them
    therefore loses an arrow; one spanning exactly one of them keeps its arrow
    count and so reaches one vertex further, which is the "extended relation" of
    the corollary.  Deleting an end vertex instead drops the dangling arrow, and
    the relation that starts (resp. ends) there has nowhere to extend to and is
    dropped.
    """
    if not 1 <= vertex <= length or length < 3:
        return None
    newLength = length - 1
    if newLength < 2:
        return None

    def mapArrow(arrow):
        """Old arrow index -> new index, or None if the arrow disappears."""
        if vertex == 1:
            return None if arrow == 1 else arrow - 1
        if vertex == length:
            return None if arrow == length - 1 else arrow
        if arrow in (vertex - 1, vertex):
            return vertex - 1
        return arrow - 1 if arrow > vertex else arrow

    moved = []
    for start, arrows in relations(relLengths):
        images = {mapArrow(arrow) for arrow in range(start, start + arrows)}
        images.discard(None)
        if len(images) < 2:
            # Fewer than two arrows left: no admissible relation remains.
            continue
        newStart, newArrows = min(images), len(images)
        if newStart + newArrows > newLength or max(images) - newStart + 1 != newArrows:
            return None
        moved.append((newStart, newArrows))

    moved = _minimalRelations(moved)
    newRelLengths = [0] * (newLength - 2)
    for start, arrows in moved:
        position = start - 1
        if position < 0 or position >= len(newRelLengths) or newRelLengths[position]:
            return None
        newRelLengths[position] = arrows
    if not _isAdmissible(newLength, newRelLengths):
        return None
    return newLength, newRelLengths


def _minimalRelations(candidates):
    """Drop any relation whose arrows contain another's -- it is not minimal.

    Extending two relations can leave one inside the other, and a Nakayama
    algebra's relations are always taken to be a minimal generating set.
    """
    kept = []
    for start, arrows in sorted(set(candidates)):
        span = set(range(start, start + arrows))
        if any(set(range(s, s + a)) < span for s, a in candidates if (s, a) != (start, arrows)):
            continue
        kept.append((start, arrows))
    return kept


def _isAdmissible(length, relLengths):
    current = relations(relLengths)
    if any(arrows < 2 for _start, arrows in current):
        return False
    if any(start + arrows > length for start, arrows in current):
        return False
    for earlier, later in zip(current, current[1:]):
        if earlier[0] >= later[0] or _end(earlier) >= _end(later):
            return False
    return True


def notPiecewiseHereditaryByDeletion(length, relLengths, cache = None):
    """Whether a certificate reaches this algebra, directly or by deleting vertices.

    Returns the chain of deletions leading to a directly certified algebra, as a
    list of (length, relation lengths, criterion name), innermost last, or None.

    The direct criteria are the base case; above them the recursion asks whether
    any one-vertex deletion is certified, memoising on (length, relation lengths)
    since many different algebras delete to the same smaller one.
    """
    cache = {} if cache is None else cache
    key = (length, tuple(relLengths))
    if key in cache:
        return cache[key]

    cache[key] = None            # guard against revisiting while recursing
    direct = certificate(length, relLengths)
    if direct is not None:
        result = [(length, list(relLengths), direct[0])]
        cache[key] = result
        return result

    for vertex in range(1, length + 1):
        smaller = removeVertex(length, relLengths, vertex)
        if smaller is None:
            continue
        chain = notPiecewiseHereditaryByDeletion(smaller[0], smaller[1], cache)
        if chain is not None:
            result = [(length, list(relLengths), 'delete vertex {0}'.format(vertex))] + chain
            cache[key] = result
            return result
    cache[key] = None
    return None
