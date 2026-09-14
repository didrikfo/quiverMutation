"""The mutation procedure on linear combinations of paths.

`mutation.quiverMutationAtVertex` implements steps 1-7 of arXiv:2112.08129 on
the set-of-paths model, where a relation is a list of paths and there are no
coefficients.  That model cannot say what the procedure's steps produce: step 4
produces a *sum*, step 5 a *difference*, and step 7 a combination whose
coefficients are the solution of a linear condition and are not recoverable from
the set of paths at all.  Reading a relation back out of the set-of-paths model
means guessing the signs, and guessing them wrong has cost real work before --
research R-003.

So this is the procedure again, with relations as `relationAlgebra`
combinations.  Nothing here guesses a sign: the coefficients come out of the
steps, and step 7 comes out of a kernel computation over the ideal.  An LNA's
own relations are single paths, so there is nothing to guess at the start
either, and the whole walk out of an LNA is exact.

The steps, on a quiver `Q` with relations `I` and a vertex `i`, writing `A` for
the targets of the arrows out of `i` and `B` for the sources of the arrows in:

1. each `beta: h -> i` and `alpha: i -> j` compose to an arrow `h -> j`;
2. each `alpha: i -> j` flips to `alpha*: j -> i*`;
3. each relation `r: i ~~> k` becomes an arrow `rbar: i* -> k`;
4. each `beta: h -> i` gives the relation `sum_alpha alpha* alpha beta = 0`,
   which in the new quiver is the sum of the paths `(h, j, i)` over `j` in `A`;
5. `rbar alpha* = r / alpha`, so `(j, i, k) - r/alpha = 0`, where `r / alpha` is
   the part of `r` whose paths begin with `alpha`, with `alpha` removed;
6. a relation ending at `i` gives one ending at each `j` in `A`, its last arrow
   `beta` replaced by the composite from step 1;
7. the relations out of `i*`, which is the kernel computation below.

A relation not meeting `i` at either end survives with `i` deleted from any path
that ran through it, since that composition is one arrow now.

**Step 7, and why it is a kernel.** A combination of paths out of `i*` is a
relation exactly when what it says about `Q` is already true there.  A path out
of `i*` is `rbar` followed by a path `s: k_r ~~> v`, and the procedure sends it
to `(r / alpha) s` for each `alpha` out of `i`.  So a combination
`sum_P eps_P P` of paths `i* ~~> v` is a relation iff, for every `alpha`,
`sum_P eps_P (r_P / alpha) s_P` lies in `I` -- which is a linear condition on
the `eps_P` over the rationals, and its solution space is a kernel.  The paper
states the case where every `s_P` is trivial; taking the kernel over all paths
out of `i*` gives that case and the longer ones together.

Cyclic quivers are out of scope here, as they are in `mutation`: step 3's cyclic
case is not implemented (research F-002), and the LNA search stops at the first
cycle.
"""

from fractions import Fraction

import networkx as nx

from . import pathAlgebra
from . import relationAlgebra as ra


# -- reading in and out of the set-of-paths model --------------------------

def relationsFrom(pathAlg):
    """The relations of a path algebra as combinations.

    Uses the coefficients the algebra is carrying if it has them and they still
    describe its relations; otherwise falls back to
    `relationAlgebra.fromPathSet`, which *guesses*: one path is zero, two paths
    are a difference, three or more a sum.

    For an LNA there is nothing to guess -- every relation is a single path --
    so a walk that starts at one and goes through `toPathAlgebra` at each step
    keeps exact coefficients the whole way.  The fallback matters for an algebra
    built by hand, and as the safe answer when something has edited `rels`
    directly and left the cache describing a different relation set.
    """
    stored = getattr(pathAlg, "relCombinations", None)
    if stored is not None and _describes(stored, pathAlg.rels):
        return [dict(relation) for relation in stored]
    return [ra.fromPathSet(rel) for rel in pathAlg.rels]


def _describes(combinations, rels):
    """Whether the cached combinations are the same relation sets as `rels`."""
    if len(combinations) != len(rels):
        return False
    return all(ra.toPathSet(c) == sorted(list(p) for p in rel)
               for c, rel in zip(combinations, rels))


def toPathAlgebra(quiver, relations):
    """A `PathAlgebra` in the set-of-paths model, for the rest of the repo.

    The coefficients are dropped, so this is lossy in exactly the way the model
    is.  Relations are sorted, and so are the paths within each, so that two
    runs that find the same algebra produce the same object.
    """
    algebra = pathAlgebra.PathAlgebra()
    algebra.add_vertices_from(list(quiver.nodes))
    algebra.add_arrows_from([[a, b] for a, b in quiver.edges()])
    ordered = sorted(((ra.toPathSet(rel), rel) for rel in relations if rel),
                     key=lambda pair: pair[0])
    algebra.add_rels_from([paths for paths, _ in ordered])
    algebra.relCombinations = [combination for _, combination in ordered]
    return algebra


# -- admissibility ---------------------------------------------------------

def isMutable(quiver, relations, vertex):
    """Whether the procedure may be applied at `vertex`.

    The paper's two conditions: no loop at the vertex, and
    `Hom(P_i*[1], Lambda) = 0`, which combinatorially is *for each nonzero path
    ending at `i` there is at least one arrow out of `i` whose composite with it
    is nonzero*.

    This is the exact reading of that condition -- "nonzero" is decided by
    `relationAlgebra.isInIdeal`, over the ideal, rather than by looking for a
    zero relation sitting inside the path.  `mutation.mutationIsPossibleAtVertex`
    is the inexact one, and is also deliberately stricter: it wants *every*
    arrow out of the vertex to keep the path nonzero, where the paper and this
    want one.  The two agree whenever the vertex has a single arrow out, which
    is the only case an LNA search meets.
    """
    if quiver.has_edge(vertex, vertex):
        return False
    outTargets = sorted(set(quiver.successors(vertex)))
    if not outTargets:
        return False
    if any(quiver.number_of_edges(a, b) > 1 for a, b in quiver.edges()):
        raise ValueError("parallel arrows cannot be named by a vertex sequence")

    for source in quiver.nodes:
        if source == vertex:
            continue
        for path in ra.allPathsBetween(quiver, source, vertex):
            if ra.isInIdeal(quiver, relations, ra.combination([path])):
                continue
            if not any(
                not ra.isInIdeal(quiver, relations,
                                 ra.combination([path + (target,)]))
                for target in outTargets
            ):
                return False
    return True


# -- the procedure ---------------------------------------------------------

def mutateAtVertex(quiver, relations, vertex):
    """Steps 1 to 7 at `vertex`, returning (new quiver, new relations).

    Does not reduce: like `mutation.quiverMutationAtVertex` it leaves the
    cleanup after the paper's step 7 to the caller, which is `reduce` below.
    Does not check admissibility either -- `isMutable` is that, and applying the
    rewrite without it gives a quiver that is not derived equivalent (R-005).
    """
    outTargets = sorted(set(quiver.successors(vertex)))
    inSources = sorted(set(quiver.predecessors(vertex)))
    outRelations = [r for r in relations if ra.source(r) == vertex]
    inRelations = [r for r in relations if ra.target(r) == vertex]
    throughRelations = [r for r in relations
                        if ra.source(r) != vertex and ra.target(r) != vertex]

    newQuiver = nx.MultiDiGraph()
    newQuiver.add_nodes_from(quiver.nodes)
    for source, target in quiver.edges():
        if source != vertex and target != vertex:
            newQuiver.add_edge(source, target)          # untouched
    for target in outTargets:
        newQuiver.add_edge(target, vertex)              # step 2: alpha*
    for source in inSources:
        for target in outTargets:
            newQuiver.add_edge(source, target)          # step 1: alpha beta
    relationArrow = {}
    for relation in outRelations:
        end = ra.target(relation)
        newQuiver.add_edge(vertex, end)                 # step 3: rbar
        relationArrow.setdefault(end, []).append(relation)

    newRelations = []

    # Step 4: one relation per arrow into the vertex.
    for source in inSources:
        newRelations.append(ra.combination(
            [(source, target, vertex) for target in outTargets]))

    # Step 5: rbar alpha* = r / alpha.
    for relation in outRelations:
        end = ra.target(relation)
        for target in outTargets:
            composite = ra.combination([((target, vertex, end), 1)])
            quotient = ra.leftDivide(relation, target)
            newRelations.append(ra.add(composite, ra.negate(quotient)))

    # Step 6: a relation into the vertex extends along each arrow out of it.
    for relation in inRelations:
        for target in outTargets:
            newRelations.append(ra.combination(
                (path[:-1] + (target,), coefficient)
                for path, coefficient in relation.items()))

    # A relation past the vertex keeps its paths, with the vertex spliced out.
    for relation in throughRelations:
        newRelations.append(ra.combination(
            (tuple(v for v in path if v != vertex), coefficient)
            for path, coefficient in relation.items()))

    # Step 7: the relations out of the mutated vertex.
    newRelations.extend(
        _relationsOutOfMutatedVertex(quiver, relations, vertex, newQuiver,
                                     outTargets, relationArrow))

    return newQuiver, normalise(newRelations)


def _relationsOutOfMutatedVertex(quiver, relations, vertex, newQuiver,
                                 outTargets, relationArrow):
    """Step 7, as a kernel, target by target.

    For a target `v`, the candidate paths out of `i*` are `rbar` followed by a
    path `k_r ~~> v` in the new quiver that does not run back through `i*`.  A
    combination of them is a relation iff every `alpha` out of `i` sends it into
    the old ideal, so the relations are the kernel of the map

        eps  |-->  ( sum_P eps_P (r_P / alpha) s_P  mod I )_alpha

    and only the part of that kernel not already forced by a relation to a
    nearer target is a new generator.
    """
    found = []
    for target in _targetsInOrder(newQuiver, vertex):
        candidates = _pathsOutOfMutatedVertex(newQuiver, vertex, target)
        if not candidates:
            continue
        shadow = {}
        for path in candidates:
            relation = _relationOfFirstArrow(relationArrow, path)
            if relation is None:
                # This path's first arrow cannot be read as one `rbar`, so the
                # step has nothing to say about paths to this target.  Other
                # targets may still be reachable through arrows that can, so
                # skip the target rather than abandoning the step.
                shadow = None
                break
            tail = path[1:]
            for arrow in outTargets:
                quotient = ra.leftDivide(relation, arrow)
                shadow[(path, arrow)] = ra.postCompose(quotient, tail) if quotient else {}
        if shadow is None:
            continue
        for element in _kernelOverIdeal(quiver, relations, candidates, outTargets, shadow):
            if not _isForcedByNearer(newQuiver, found, element, vertex):
                found.append(element)
    return found


def _targetsInOrder(newQuiver, vertex):
    """Every vertex reachable from the mutated vertex, nearest first."""
    lengths = nx.single_source_shortest_path_length(newQuiver, vertex)
    return [v for v, _ in sorted(lengths.items(), key=lambda kv: (kv[1], kv[0]))
            if v != vertex]


def _pathsOutOfMutatedVertex(newQuiver, vertex, target):
    """Paths from the mutated vertex to `target` that leave and do not return."""
    return [p for p in ra.allPathsBetween(newQuiver, vertex, target)
            if len(p) > 1 and vertex not in p[1:]]


def _relationOfFirstArrow(relationArrow, path):
    """The relation the path's first arrow came from, or None if it is not one.

    Step 3 can give two arrows `i* -> k` from two relations `i ~~> k`, which a
    vertex sequence cannot tell apart; when that happens there is no reading of
    the path as `rbar s` and step 7 has nothing to say.
    """
    candidates = relationArrow.get(path[1])
    if not candidates or len(candidates) > 1:
        return None
    return candidates[0]


def _kernelOverIdeal(quiver, relations, candidates, outTargets, shadow):
    """The eps with sum_P eps_P shadow[P, alpha] in the ideal, for every alpha.

    Solved as a linear system over the rationals: reduce each shadow against a
    basis of the ideal between its endpoints, and the residues are the rows.
    """
    rows = {}
    for arrow in outTargets:
        residues = {}
        for path in candidates:
            element = shadow.get((path, arrow)) or {}
            residues[path] = _reduceAgainstIdeal(quiver, relations, element, arrow)
        basis = sorted({key for r in residues.values() for key in r})
        for key in basis:
            rows[(arrow, key)] = {path: residues[path].get(key, Fraction(0))
                                  for path in candidates}
    if not rows:
        # Every candidate is sent to zero, so every combination of them is a
        # relation; the single paths generate all of that.
        return [ra.combination([(path, 1)]) for path in candidates]
    return [
        ra.combination((path, _asInteger(coefficient))
                       for path, coefficient in solution.items() if coefficient)
        for solution in _nullSpace(candidates, rows)
    ]


def _reduceAgainstIdeal(quiver, relations, element, arrow):
    """`element` reduced modulo the ideal, as a dict of residual coefficients."""
    if not element:
        return {}
    pivots = ra.idealBasis(quiver, relations, ra.source(element), ra.target(element))
    row = {k: Fraction(v) for k, v in element.items()}
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


def _nullSpace(columns, rows):
    """A basis of the null space of the matrix given as rows over `columns`."""
    order = list(columns)
    matrix = [[Fraction(row.get(column, 0)) for column in order]
              for row in rows.values()]
    pivotOf = {}
    pivotRows = []
    for row in matrix:
        row = list(row)
        for column, head in pivotOf.items():
            if row[column]:
                factor = row[column]
                row = [a - factor * b for a, b in zip(row, pivotRows[head])]
        lead = next((i for i, v in enumerate(row) if v), None)
        if lead is None:
            continue
        row = [v / row[lead] for v in row]
        for index, existing in enumerate(pivotRows):
            if existing[lead]:
                factor = existing[lead]
                pivotRows[index] = [a - factor * b for a, b in zip(existing, row)]
        pivotOf[lead] = len(pivotRows)
        pivotRows.append(row)

    free = [i for i in range(len(order)) if i not in pivotOf]
    basis = []
    for index in free:
        solution = [Fraction(0)] * len(order)
        solution[index] = Fraction(1)
        for column, head in pivotOf.items():
            solution[column] = -pivotRows[head][index]
        denominators = [v.denominator for v in solution if v]
        scale = 1
        for d in denominators:
            scale = scale * d // _gcd(scale, d)
        basis.append({order[i]: solution[i] * scale for i in range(len(order))})
    return basis


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def normalise(relations):
    """Integral rational coefficients back to ints.

    The substitution in `reduce` divides, so coefficients pass through
    `Fraction`; anything that comes out integral should say so, or two equal
    relations compare unequal and print differently.
    """
    cleaned = []
    for relation in relations:
        cleaned.append({
            path: int(coefficient) if Fraction(coefficient).denominator == 1
            else Fraction(coefficient)
            for path, coefficient in relation.items()
        })
    return [r for r in cleaned if r]


def _asInteger(value):
    value = Fraction(value)
    if value.denominator != 1:
        raise ValueError("a relation's coefficients should be integers, got {0}".format(value))
    return int(value)


def _isForcedByNearer(newQuiver, found, element, vertex):
    """Whether `element` is already a consequence of a relation to a nearer target."""
    if not found:
        return False
    target = ra.target(element)
    spanning = []
    for relation in found:
        end = ra.target(relation)
        if end == target:
            spanning.append(relation)
            continue
        for after in ra.allPathsBetween(newQuiver, end, target):
            if len(after) > 1 and vertex not in after[1:]:
                spanning.append(ra.postCompose(relation, after))
    if not spanning:
        return False
    pivots = ra._rowReduce([dict(s) for s in spanning])
    row = {k: Fraction(v) for k, v in element.items()}
    while row:
        head = min(row)
        if head not in pivots:
            return False
        factor = row[head]
        pivotRow = pivots[head]
        row = {k: row.get(k, Fraction(0)) - factor * pivotRow.get(k, Fraction(0))
               for k in set(row) | set(pivotRow)}
        row = {k: v for k, v in row.items() if v != 0}
    return True


# -- the cleanup after step 7 ---------------------------------------------

def reduce(quiver, relations):
    """The paper's "Note" after step 7, then minimality, both exactly.

    The Note: a relation containing a path of length one says that arrow equals
    a combination of longer paths, so the arrow and the relation both go, and
    every other relation has that arrow substituted out.  On combinations this
    is arithmetic rather than list surgery -- `c alpha + rest = 0` means
    `alpha = -rest / c`, and splicing that into a path is `preCompose` and
    `postCompose` -- which is the part the set-of-paths model could only
    approximate.

    Minimality: a relation is dropped when it lies in the ideal generated by the
    others, decided by linear algebra over that ideal.  That is the exact
    version of `reduction.removeRedundantRelations` and
    `reduction.removeExistingSubrelations`, which look for a syntactic
    substitution instead and so can only find the redundancies that happen to
    be visible as one.
    """
    quiver = nx.MultiDiGraph(quiver)
    relations = [dict(r) for r in relations if r]
    while True:
        substituted = _substituteOutAnArrow(quiver, relations)
        if substituted is None:
            break
        quiver, relations = substituted
    return quiver, normalise(_minimalGenerators(quiver, relations))


def _substituteOutAnArrow(quiver, relations):
    """One pass of the Note, or None when no relation has a length-one path."""
    for index, relation in enumerate(relations):
        singles = sorted(p for p in relation if len(p) == 2)
        if not singles:
            continue
        arrow = singles[0]
        coefficient = relation[arrow]
        rest = {p: c for p, c in relation.items() if p != arrow}
        replacement = ra.scale(rest, Fraction(-1, 1) / coefficient) if rest else {}

        newQuiver = nx.MultiDiGraph(quiver)
        if newQuiver.has_edge(*arrow):
            newQuiver.remove_edge(*arrow)
        newRelations = []
        for other in relations[:index] + relations[index + 1:]:
            rewritten = _substituteArrow(other, arrow, replacement)
            if rewritten:
                newRelations.append(rewritten)
        return newQuiver, newRelations
    return None


def _substituteArrow(relation, arrow, replacement):
    """`relation` with every occurrence of `arrow` replaced by `replacement`.

    An empty replacement means the arrow is zero in the algebra, so every path
    through it drops out.
    """
    while True:
        for path, coefficient in relation.items():
            position = _firstOccurrence(path, arrow)
            if position is None:
                continue
            others = {p: c for p, c in relation.items() if p != path}
            spliced = {}
            for piece, pieceCoefficient in replacement.items():
                grafted = path[:position] + piece + path[position + 2:]
                spliced[grafted] = spliced.get(grafted, 0) + coefficient * pieceCoefficient
            relation = ra.add(others, spliced)
            break
        else:
            return relation


def _firstOccurrence(path, arrow):
    """The index where `arrow` sits inside `path` as two consecutive vertices."""
    for position in range(len(path) - 1):
        if (path[position], path[position + 1]) == arrow:
            return position
    return None


def _minimalGenerators(quiver, relations):
    """The relations that the others do not already imply, simplest first.

    Sorted by longest path then by number of paths, so that of two relations
    that imply each other the simpler one is the one kept -- which is what makes
    the result stable enough to compare two runs.
    """
    ordered = sorted(
        (r for r in relations if r),
        key=lambda r: (max(len(p) for p in r), len(r), sorted(r)),
    )
    kept = []
    for index, relation in enumerate(ordered):
        others = kept + ordered[index + 1:]
        if not ra.isInIdeal(quiver, others, relation):
            kept.append(relation)
    return kept


def mutateAtVertices(quiver, relations, vertices):
    """Mutate at each vertex in turn, reducing after each.  Negative means left.

    A left mutation is a right mutation of the opposite algebra, read back --
    the same identity `mutation.leftQuiverMutationAtVertex` uses.
    """
    for vertex in vertices:
        if vertex >= 0:
            quiver, relations = mutateAtVertex(quiver, relations, vertex)
        else:
            quiver, relations = _opposite(*mutateAtVertex(
                *_opposite(quiver, relations), -vertex))
        quiver, relations = reduce(quiver, relations)
    return quiver, relations


def _opposite(quiver, relations):
    """Every arrow and every path reversed, coefficients kept."""
    reversed_ = nx.MultiDiGraph()
    reversed_.add_nodes_from(quiver.nodes)
    for source, target in quiver.edges():
        reversed_.add_edge(target, source)
    return reversed_, [
        ra.combination((path[::-1], coefficient) for path, coefficient in relation.items())
        for relation in relations
    ]


def mutateLeftAtVertex(quiver, relations, vertex):
    """Left mutation: the same procedure on the opposite algebra, read back."""
    return _opposite(*mutateAtVertex(*_opposite(quiver, relations), vertex))
