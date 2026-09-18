"""Relations as combinations of paths that name their arrows.

`relationAlgebra` models a path as a sequence of **vertices**, which is enough
for every quiver this repo started from and not enough for the quivers the
mutation procedure produces.  Step 1 of arXiv:2112.08129 adds a composite arrow
`alpha beta: h -> j` for every `beta: h -> i` and `alpha: i -> j`, and nothing
stops `h -> j` from being an arrow already; step 3 adds one arrow `i* -> k` per
relation `i ~~> k`, and nothing stops two relations from sharing both ends.  So
**parallel arrows are produced by the procedure**, and a vertex sequence cannot
say which of two arrows with the same endpoints a path uses.

This module is the same algebra over a representation that can.  An **arrow** is
a triple `(tail, head, key)` -- exactly a `networkx.MultiDiGraph` edge, key and
all -- and a **path** is a tuple of arrows, so

    ((1, 2, 0), (2, 6, 1))

is the path along the arrow `1 -> 2` and then the *second* arrow `2 -> 6`.  The
empty tuple is a trivial path, at whichever vertex the context supplies, and is
the identity for the composition, which is plain concatenation.  A path carries
its own source and target, so nothing here needs the quiver to read a path back.

What it buys, beyond being able to state the quivers at all:

* **Step 5 divides by an arrow**, `r / alpha`, which is what the paper says.  In
  the vertex model it divided by the *target* of `alpha`, which silently lumped
  parallel arrows together.
* **Step 7 stops giving up.**  Its candidates are `rbar` followed by a tail, so
  it has to read a path's first arrow back as the relation it came from; where
  two relations `i ~~> k` gave two arrows `i* -> k`, the vertex model could not,
  and `procedure` skipped the step for that target.
* **The Cartan matrix counts the paths there are.**  Two parallel arrows are two
  paths, and the vertex model saw one, so every derived invariant read off it was
  wrong exactly where the quiver had parallel arrows (research F-038's second
  mechanism).

The linear algebra is the same as in `relationAlgebra` and is not repeated:
`rowReduce` there is over an arbitrary ordered basis, and a path of arrows is as
good a basis element as a path of vertices.

`lift` and `project` are the two directions between this and the vertex model.
Lifting is possible exactly when the quiver has no parallel arrows -- which is
every quiver the repo *starts* a walk from -- and projecting is lossy exactly
when it has, which is why `pathAlgebra.PathAlgebra` carries `arrowRels` and the
procedure reads that in preference to `rels`.
"""

from fractions import Fraction

import networkx as nx

from . import relationAlgebra as ra


# -- arrows and paths -----------------------------------------------------

def arrowsOf(quiver):
    """Every arrow of a quiver, as (tail, head, key) triples, sorted."""
    return sorted(quiver.edges(keys = True))


def arrowsOutOf(quiver, vertex):
    """The arrows with tail `vertex`, sorted."""
    return sorted(quiver.out_edges(vertex, keys = True))


def arrowsInto(quiver, vertex):
    """The arrows with head `vertex`, sorted."""
    return sorted(quiver.in_edges(vertex, keys = True))


def hasParallelArrows(quiver):
    """Whether two distinct arrows of the quiver share both endpoints."""
    return quiver.number_of_edges() != len({(tail, head)
                                            for tail, head, _ in quiver.edges(keys = True)})


def pathSource(path):
    """The vertex a non-empty path starts at."""
    if not path:
        raise ValueError("a trivial path does not carry its vertex")
    return path[0][0]


def pathTarget(path):
    if not path:
        raise ValueError("a trivial path does not carry its vertex")
    return path[-1][1]


def pathVertices(path):
    """The vertices a non-empty path visits, in order."""
    return [path[0][0]] + [arrow[1] for arrow in path]


def isPath(quiver, path):
    """Whether the arrows of `path` are arrows of the quiver, end to end."""
    for arrow in path:
        tail, head, key = arrow
        if not quiver.has_edge(tail, head, key):
            return False
    return all(path[i][1] == path[i + 1][0] for i in range(len(path) - 1))


def allPathsBetween(quiver, sourceVertex, targetVertex):
    """Every path from source to target, as arrow tuples, shortest first.

    Simple in the *vertices*, as everything in this repo is: a path does not
    revisit a vertex, so on a quiver with a cycle this enumerates the simple
    paths and stops.  Two parallel arrows give two paths, which is the whole
    point -- `relationAlgebra.allPathsBetween` gives one.
    """
    if sourceVertex == targetVertex:
        return [()]
    found = []

    def walk(current, visited, prefix):
        for tail, head, key in sorted(quiver.out_edges(current, keys = True)):
            if head in visited:
                continue
            arrow = (tail, head, key)
            if head == targetVertex:
                found.append(prefix + (arrow,))
            else:
                walk(head, visited | {head}, prefix + (arrow,))

    walk(sourceVertex, frozenset([sourceVertex]), ())
    return sorted(found, key = lambda p: (len(p), p))



# -- combinations ---------------------------------------------------------

def combination(terms):
    """A linear combination of arrow paths, as a dict path -> nonzero coefficient.

    `terms` is a mapping, an iterable of (path, coefficient) pairs, or an
    iterable of bare paths taken with coefficient 1.  Repeated paths are summed
    and terms that cancel are dropped, as in `relationAlgebra.combination`.

    A pair is told from a bare path by its *second* element: a coefficient is a
    number and an arrow is a tuple, and a bare path of two arrows is itself a
    two-element tuple, so the length cannot decide it.
    """
    if hasattr(terms, "items"):
        terms = terms.items()
    result = {}
    for term in terms:
        if (isinstance(term, tuple) and len(term) == 2
                and not isinstance(term[1], tuple)):
            path, coefficient = term
        else:
            path, coefficient = term, 1
        path = tuple(path)
        for arrow in path:
            if not (isinstance(arrow, tuple) and len(arrow) == 3):
                raise ValueError("a path is a tuple of (tail, head, key) arrows, "
                                 "got {0!r}".format(path))
        result[path] = result.get(path, 0) + coefficient
    return {path: c for path, c in result.items() if c != 0}


def add(*combs):
    terms = []
    for comb in combs:
        terms.extend(comb.items())
    return combination(terms)


def scale(comb, factor):
    return combination((path, coefficient * factor) for path, coefficient in comb.items())


def negate(comb):
    return scale(comb, -1)


def source(comb):
    """The common source vertex of every path, or None for the zero combination."""
    return pathSource(next(iter(comb))) if comb else None


def target(comb):
    return pathTarget(next(iter(comb))) if comb else None


def isHomogeneous(comb):
    """Whether every path in the combination has the same source and target."""
    if not comb:
        return True
    paths = list(comb)
    if any(not path for path in paths):
        return False
    return all(pathSource(p) == pathSource(paths[0]) and pathTarget(p) == pathTarget(paths[0])
               for p in paths)


def leftDivide(comb, arrow):
    """`r / alpha`: the paths of r that begin with the arrow, with it removed.

    The paper's own notation, and the point of naming arrows: in the vertex
    model this could only be `r` divided by *an* arrow to a given vertex.
    """
    return combination(
        (path[1:], coefficient)
        for path, coefficient in comb.items()
        if path and path[0] == arrow
    )


def rightDivide(comb, arrow):
    """The mirror of leftDivide: drop the last arrow, which must be `arrow`."""
    return combination(
        (path[:-1], coefficient)
        for path, coefficient in comb.items()
        if path and path[-1] == arrow
    )


def preCompose(comb, path):
    """`path` then `comb`: every path of comb with `path` glued on the front."""
    prefix = tuple(path)
    return combination((prefix + p, c) for p, c in comb.items())


def postCompose(comb, path):
    suffix = tuple(path)
    return combination((p + suffix, c) for p, c in comb.items())


def substituteArrow(comb, arrow, replacement):
    """Every occurrence of `arrow` in `comb` replaced by the combination `replacement`.

    An empty replacement means the arrow is zero in the algebra, so every path
    through it drops out.  Runs to a fixed point, since a replacement may itself
    contain the arrow's endpoints but never the arrow, which has been removed
    from the quiver by the time this is called.
    """
    while True:
        for path, coefficient in comb.items():
            if arrow not in path:
                continue
            position = path.index(arrow)
            others = {p: c for p, c in comb.items() if p != path}
            spliced = {}
            for piece, pieceCoefficient in replacement.items():
                grafted = path[:position] + piece + path[position + 1:]
                spliced[grafted] = spliced.get(grafted, 0) + coefficient * pieceCoefficient
            comb = add(others, spliced)
            break
        else:
            return comb


# -- ideals ---------------------------------------------------------------

def idealSpanningSet(quiver, relations, sourceVertex, targetVertex):
    """The combinations `v r u` spanning the ideal's part between two vertices."""
    spanning = []
    for relation in relations:
        if not relation:
            continue
        relationSource, relationTarget = source(relation), target(relation)
        for before in allPathsBetween(quiver, sourceVertex, relationSource):
            extendedLeft = preCompose(relation, before)
            for after in allPathsBetween(quiver, relationTarget, targetVertex):
                spanning.append(postCompose(extendedLeft, after))
    return [element for element in spanning if element]


def idealBasis(quiver, relations, sourceVertex, targetVertex):
    """A row-reduced basis of the ideal's part between two vertices.

    The row reduction is `relationAlgebra.rowReduce`, which is over an ordered
    basis of whatever the paths are; only the basis elements differ.
    """
    spanning = idealSpanningSet(quiver, relations, sourceVertex, targetVertex)
    return ra._rowReduce([dict(element) for element in spanning])


def reduceAgainstPivots(comb, pivots):
    """`comb` reduced modulo a row-reduced set of pivots, as residual coefficients."""
    row = {k: Fraction(v) for k, v in comb.items()}
    while row:
        head = min(row)
        if head not in pivots:
            break
        factor = row[head]
        pivotRow = pivots[head]
        row = {k: row.get(k, Fraction(0)) - factor * pivotRow.get(k, Fraction(0))
               for k in set(row) | set(pivotRow)}
        row = {k: v for k, v in row.items() if v != 0}
    return row


def isInIdeal(quiver, relations, comb):
    """Whether a combination of arrow paths is zero in the algebra."""
    if not comb:
        return True
    if not isHomogeneous(comb):
        raise ValueError("a combination must have one source and one target: {0}".format(comb))
    pivots = idealBasis(quiver, relations, source(comb), target(comb))
    return not reduceAgainstPivots(comb, pivots)


def homDimension(quiver, relations, sourceVertex, targetVertex):
    """dim of the span of the paths source -> target in the algebra.

    The number of arrow paths between the two vertices, minus the rank of the
    ideal between them.  For source == target the trivial path contributes 1 and
    is never in an admissible ideal, so the result is at least 1 there.
    """
    paths = allPathsBetween(quiver, sourceVertex, targetVertex)
    if not paths:
        return 0
    rank = len(idealBasis(quiver, relations, sourceVertex, targetVertex))
    return len(paths) - rank


# -- the cheap count, for the invariant the search compares --------------
#
# `homDimension` is exact and costs a row reduction per pair of vertices, which
# is too much for something a search calls at every node.  The cheap count is
# what `paths.numberOfPathsUpToRels` does, over arrow paths: classes of paths
# under substituting one side of a two-path relation for the other, and a class
# is zero when any of its members runs through a one-path relation.  For a
# monomial ideal that is exact, and every algebra the classification starts from
# has one.

def _occurrences(path, piece):
    """Where `piece` sits inside `path` as a contiguous run of arrows."""
    if not piece:
        return []
    return [i for i in range(len(path) - len(piece) + 1)
            if path[i:i + len(piece)] == piece]


def _containsZero(path, zeroPaths):
    return any(_occurrences(path, zero) for zero in zeroPaths)


def homDimensionByClosure(quiver, relations, sourceVertex, targetVertex):
    """The cheap count of the paths source -> target that no relation kills.

    Paths are joined when one becomes the other by substituting one side of a
    two-path relation for the other at some position, the join is closed
    transitively, and a class is dropped when any member runs through a one-path
    relation.  The closure is complete here, where `numberOfPathsUpToRels`
    approximates it by applying each relation at most once per pass.
    """
    if sourceVertex == targetVertex:
        # The trivial path, which no admissible relation touches, and nothing
        # else: closed paths are not counted, exactly as neither
        # `relationAlgebra.homDimension` nor `paths.numberOfPathsUpToRels` counts
        # them -- `nx.all_simple_paths(q, v, v)` yields only the trivial path.
        # On an acyclic quiver, which is what the classification works over,
        # there is nothing to count.  On a cyclic one the diagonal is then 1, the
        # Cartan matrix is not unimodular, and `coxeterCoefficients` raises --
        # which is what `search._coxeterKeyOrNone` relies on to end a branch at a
        # cycle rather than compare a meaningless key.
        return 1
    paths = allPathsBetween(quiver, sourceVertex, targetVertex)
    if not paths:
        return 0
    return _classesNotKilled(paths, relations)


def _classesNotKilled(paths, relations):
    """How many classes of `paths` under the commutativity relations survive."""
    zeroPaths = [next(iter(relation)) for relation in relations if len(relation) == 1]
    pairs = [tuple(relation) for relation in relations if len(relation) == 2]

    parent = {path: path for path in paths}

    def find(path):
        while parent[path] != path:
            parent[path] = parent[parent[path]]
            path = parent[path]
        return path

    def union(one, other):
        one, other = find(one), find(other)
        if one != other:
            parent[max(one, other)] = min(one, other)

    known = set(paths)
    for path in paths:
        for first, second in pairs:
            for source_, replacement in ((first, second), (second, first)):
                for position in _occurrences(path, source_):
                    grafted = path[:position] + replacement + path[position + len(source_):]
                    if grafted in known:
                        union(path, grafted)

    zeroClasses = {find(path) for path in paths if _containsZero(path, zeroPaths)}
    classes = {find(path) for path in paths}
    return len(classes - zeroClasses)


def cartanMatrix(quiver, relations, exact = False):
    """The Cartan matrix as nested lists of ints, vertices sorted.

    Entry (j, i) is dim e_j (kQ/I) e_i.  `exact` picks `homDimension` over
    `homDimensionByClosure`; the two agree on a monomial ideal and the cheap one
    is what the search can afford.
    """
    dimension = homDimension if exact else homDimensionByClosure
    vertices = sorted(quiver.nodes)
    index = {vertex: position for position, vertex in enumerate(vertices)}
    size = len(vertices)
    matrix = [[0] * size for _ in range(size)]
    for sourceVertex in vertices:
        for targetVertex in vertices:
            matrix[index[targetVertex]][index[sourceVertex]] = dimension(
                quiver, relations, sourceVertex, targetVertex)
    return matrix


# -- between this and the vertex model ------------------------------------

def liftPath(quiver, vertexPath):
    """A vertex sequence read as an arrow path.

    Raises when a step is not an arrow of the quiver, and when it is more than
    one -- which is the whole content of the restriction this module lifts: a
    vertex sequence names a path only while the quiver has no parallel arrows
    along it.
    """
    vertexPath = list(vertexPath)
    arrows = []
    for tail, head in zip(vertexPath, vertexPath[1:]):
        keys = sorted(quiver[tail][head]) if quiver.has_edge(tail, head) else []
        if not keys:
            raise ValueError("{0} -> {1} is not an arrow of the quiver".format(tail, head))
        if len(keys) > 1:
            raise ValueError("{0} -> {1} is {2} parallel arrows, so the vertex "
                             "sequence {3} does not name a path".format(
                                 tail, head, len(keys), vertexPath))
        arrows.append((tail, head, keys[0]))
    return tuple(arrows)


def projectPath(path):
    """An arrow path as the vertex sequence it runs through.

    Lossy exactly where the quiver has parallel arrows, which is why this is
    only ever the storage format.
    """
    return pathVertices(path)


def liftCombination(quiver, vertexComb):
    """A combination of vertex paths as one of arrow paths."""
    return combination((liftPath(quiver, path), coefficient)
                       for path, coefficient in vertexComb.items())


def projectCombination(comb):
    """A combination of arrow paths as one of vertex paths, dropping the keys.

    Two parallel paths project to the same vertex sequence and their
    coefficients are then added, which can cancel a relation that is not zero.
    The caller has to know that; `describesRels` is how `procedure` checks that a
    projection it is about to trust is still faithful.
    """
    return ra.combination((tuple(projectPath(path)), coefficient)
                          for path, coefficient in comb.items())


def lift(quiver, rels, reading = None):
    """The relations of the vertex model as arrow combinations.

    `reading` turns one relation in set-of-paths form into a combination of
    vertex paths; the default is `relationAlgebra.fromPathSet`, the reading the
    repo means.
    """
    reading = ra.fromPathSet if reading is None else reading
    return [liftCombination(quiver, reading(rel)) for rel in rels]


def projectToPathSets(relations):
    """Arrow relations as the `rels` of the vertex model: sorted vertex paths."""
    return [sorted(projectPath(path) for path in relation) for relation in relations]


def describesRels(relations, rels):
    """Whether these arrow relations project onto exactly `rels`.

    The check `procedure.relationsFrom` makes before trusting a cached set of
    arrow relations against the `rels` the rest of the repo keys tables by.
    """
    if len(relations) != len(rels):
        return False
    return all(sorted(projectPath(path) for path in relation)
               == sorted(list(p) for p in rel)
               for relation, rel in zip(relations, rels))


def usesOnlyArrowsOf(quiver, relations):
    """Whether every arrow named by these relations is an arrow of the quiver."""
    return all(isPath(quiver, path) for relation in relations for path in relation)


# -- malformed relations --------------------------------------------------

def isIllegalRelation(quiver, relation):
    """Whether a relation is malformed: not a path, not homogeneous, repeated.

    The arrow-model reading of `paths.isIllegalRelation`, and it says *no* in two
    places where the vertex model said yes.  Two parallel paths are two paths, so
    a relation between them is not a repeat; and a path through a parallel pair
    is a path, which `nx.is_path` on a vertex sequence cannot decide.
    """
    paths = list(relation)
    if not paths:
        return False
    if any(not path for path in paths):
        return True
    if not isHomogeneous(relation):
        return True
    for path in paths:
        if not isPath(quiver, path):
            return True
        vertices = pathVertices(path)
        if len(set(vertices)) != len(vertices):
            return True
    return False


def describeQuiver(quiver, relations):
    """A short readable dump, for the times a parallel arrow has to be looked at."""
    lines = ["vertices: {0}".format(sorted(quiver.nodes)),
             "arrows:   {0}".format(arrowsOf(quiver))]
    for relation in relations:
        lines.append("relation: {0}".format(
            " + ".join("{0}*{1}".format(coefficient, projectPath(path))
                       for path, coefficient in sorted(relation.items()))))
    return "\n".join(lines)
