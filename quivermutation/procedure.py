"""The mutation procedure on linear combinations of paths that name their arrows.

Two things the set-of-paths model cannot say, and this module says both.

**Coefficients.**  Step 4 of arXiv:2112.08129 produces a *sum*, step 5 a
*difference*, and step 7 a combination whose coefficients solve a linear
condition and are not recoverable from a set of paths at all.  Guessing them
back has cost real work before -- research R-003, R-007.

**Which arrow.**  The procedure *produces parallel arrows*: step 1 adds a
composite `alpha beta: h -> j` whether or not `h -> j` is an arrow already, and
step 3 adds one arrow `i* -> k` per relation `i ~~> k`, so two relations sharing
both ends give two.  A path as a sequence of vertices cannot say which of two
arrows with the same endpoints it runs along, so until 2026-09-18 the engine ran
on vertex sequences and three things went wrong: step 5 divided a relation by the
*target* of an arrow rather than by the arrow, step 7 gave up on a target
whenever two relations led to it, and a relation carried past the mutated vertex
lost the difference between the new composite `alpha beta` and an arrow `h -> j`
that was there all along.  The relations are `arrowPaths` combinations now, over
paths that are tuples of `(tail, head, key)` arrows, and all three are exact.

The steps, on a quiver `Q` with relations `I` and a vertex `i`, writing `A` for
the arrows out of `i` and `B` for the arrows into it:

1. each `beta: h -> i` and `alpha: i -> j` compose to an arrow `alpha beta: h -> j`,
   one per *pair*;
2. each `alpha: i -> j` flips to `alpha*: j -> i*`;
3. each relation `r: i ~~> k` becomes an arrow `rbar: i* -> k`;
4. each `beta: h -> i` gives the relation `sum_alpha alpha* (alpha beta) = 0`;
5. `rbar alpha* = r / alpha`, where `r / alpha` is the part of `r` whose paths
   begin with the arrow `alpha`, with `alpha` removed;
6. a relation ending at `i` gives one ending at each `t(alpha)`, its last arrow
   `beta` replaced by the composite `alpha beta` of step 1;
7. the relations out of `i*`, which is the kernel computation below.

A relation not meeting `i` at either end survives with the two arrows it ran
through `i` along, `beta` then `alpha`, replaced by the one composite arrow
`alpha beta`.

**Step 7, and why it is a kernel.**  A combination of paths out of `i*` is a
relation exactly when what it says about `Q` is already true there.  A path out
of `i*` is `rbar` followed by a path `s: k_r ~~> v` in the new quiver, and the
procedure sends it to `(r / alpha) s` for each `alpha` out of `i`.  So a
combination `sum_P eps_P P` of paths `i* ~~> v` is a relation iff, for every
`alpha`, `sum_P eps_P (r_P / alpha) s_P` lies in `I` -- a linear condition on the
`eps_P` over the rationals, whose solution space is a kernel.  The paper states
the case where every `s_P` is trivial; taking the kernel over all paths out of
`i*` gives that case and the longer ones together.

The tail `s_P` is a path of the **new** quiver and the condition is about the
**old** one, so it is read back arrow by arrow before the membership test: a
carried arrow is itself, and a composite `alpha beta` is the two old arrows
`beta` then `alpha`.  A tail out of `i*` that never returns to `i*` contains no
flip and no `rbar` after the first, so it is always one or the other and the
reading is total.  `_inOldQuiver` is it.

Cyclic quivers are still out of scope: step 3's cyclic case is not implemented
(research F-002), and the LNA search stops at the first cycle.
"""

from fractions import Fraction

import networkx as nx

from . import arrowPaths as ap
from . import pathAlgebra
from . import relationAlgebra as ra


# -- reading in and out of the set-of-paths model --------------------------

def relationsFrom(pathAlg):
    """The relations of a path algebra as combinations of arrow paths.

    Prefers the arrow relations the algebra is carrying, `arrowRels`, which is
    what a mutation leaves behind and the only faithful record once the quiver
    has parallel arrows.  Failing that it lifts `rels` -- the vertex sequences
    the tables are keyed by -- which is possible exactly while no relation runs
    along a parallel pair, and guesses the coefficients the way
    `relationAlgebra.fromPathSet` does: one path is zero, two are a difference,
    three or more a sum.

    For an LNA there is nothing to guess and nothing parallel, so a walk that
    starts at one and goes through `toPathAlgebra` at each step is exact the
    whole way.  The lift is what makes an algebra built by hand work, and is the
    safe answer when something has edited `rels` directly.
    """
    stored = getattr(pathAlg, "arrowRels", None)
    if (stored is not None and ap.usesOnlyArrowsOf(pathAlg.quiver, stored)
            and ap.describesRels(stored, pathAlg.rels)):
        return [dict(relation) for relation in stored]
    legacy = getattr(pathAlg, "relCombinations", None)
    if legacy is not None and _describes(legacy, pathAlg.rels):
        return [ap.liftCombination(pathAlg.quiver, dict(relation)) for relation in legacy]
    return ap.lift(pathAlg.quiver, pathAlg.rels)


def _describes(combinations, rels):
    """Whether the cached vertex-model combinations are the same sets as `rels`."""
    if len(combinations) != len(rels):
        return False
    return all(ra.toPathSet(c) == sorted(list(p) for p in rel)
               for c, rel in zip(combinations, rels))


def toPathAlgebra(quiver, relations):
    """A `PathAlgebra` carrying both models: `arrowRels` and the projected `rels`.

    `rels` is the vertex-sequence form the tables are keyed by and two algebras
    are compared as, and it is **lossy exactly when the quiver has parallel
    arrows** -- two parallel paths write down as the same vertex sequence.  So
    `arrowRels` is set alongside it and is what `relationsFrom` reads; the
    projection is kept because a class name has to be a string.

    The quiver is copied with its keys intact, since those keys are the arrow
    names the relations are written in.  Relations are sorted by their
    projection, then by the arrow paths themselves, so that two runs that find
    the same algebra produce the same object.
    """
    algebra = pathAlgebra.PathAlgebra()
    algebra.quiver = nx.MultiDiGraph(quiver)
    ordered = sorted(
        ((sorted(ap.projectPath(path) for path in relation), relation)
         for relation in relations if relation),
        key = lambda pair: (pair[0], sorted(pair[1])),
    )
    algebra.add_rels_from([paths for paths, _ in ordered])
    algebra.arrowRels = [relation for _, relation in ordered]
    algebra.relCombinations = [ap.projectCombination(relation) for _, relation in ordered]
    return algebra


# -- admissibility ---------------------------------------------------------

def isMutable(quiver, relations, vertex, allowParallelArrows = True):
    """Whether the procedure may be applied at `vertex`.  The search's gate.

    The paper's theorem names two cases where mutation is impossible: no arrow
    out of the vertex at all, and *there is a nonzero path ending at the vertex
    whose composite with **every** arrow out of it is zero*.  This is that, read
    exactly: "nonzero" is decided by `arrowPaths.isInIdeal`, over the ideal,
    rather than by looking for a zero relation sitting inside the path, and
    "every arrow out" means every arrow, so two parallel arrows out of the vertex
    count twice.

    Two things it is worth being precise about.

    The criterion **rules mutation out, it does not rule it in.** The theorem's
    own hypothesis is `Hom(P_i*[1], Lambda) = 0`, and the paper says plainly
    that this is in general *not* equivalent to a condition on the quiver.  So
    passing this is necessary, not sufficient, and a rewrite performed on the
    strength of it can still fail to be a derived equivalence -- which is why
    research R-005 exists and why F-016 checks the Coxeter polynomial across
    every mutation this allows that its predecessor did not.

    **Parallel arrows are no longer refused.**  They used to be, and the reason
    given was honest about what it was: a restriction of this repo's model, where
    a path was a sequence of vertices, and not a restriction of the procedure.
    Since the relations name their arrows there is nothing left to refuse, and
    `allowParallelArrows = False` restores the old gate for measuring what the
    change does -- not for producing answers.
    """
    if quiver.has_edge(vertex, vertex):
        return False
    outArrows = ap.arrowsOutOf(quiver, vertex)
    if not outArrows:
        return False
    if not allowParallelArrows and ap.hasParallelArrows(quiver):
        return False

    for sourceVertex in quiver.nodes:
        if sourceVertex == vertex:
            continue
        for path in ap.allPathsBetween(quiver, sourceVertex, vertex):
            if ap.isInIdeal(quiver, relations, ap.combination([path])):
                continue
            if not any(
                not ap.isInIdeal(quiver, relations, ap.combination([path + (arrow,)]))
                for arrow in outArrows
            ):
                return False
    return True


# -- the procedure ---------------------------------------------------------

class _Arrows:
    """The arrows of the mutated quiver, named by where they came from.

    Steps 1 to 3 each add arrows, and steps 4 to 7 have to refer to them by
    name: step 4 to the composite of a given `beta` with a given `alpha`, step 5
    to a given `alpha*` and `rbar`, step 6 to the composite again.  With
    parallel arrows a name is no longer `(tail, head)`, so this hands out the
    `(tail, head, key)` triples and remembers which is which.

    Keys are assigned per pair of endpoints in one deterministic order -- the
    arrows carried over first, in their old order, then the flips, the
    composites and the relation arrows -- so the same mutation twice gives the
    same arrow names and two runs produce comparable algebras.
    """

    CARRIED, FLIP, COMPOSITE, RBAR = 0, 1, 2, 3

    def __init__(self):
        self._entries = []
        self._arrow = {}

    def declare(self, tag, tail, head, rank, order):
        self._entries.append((tail, head, rank, order, tag))

    def build(self, vertices):
        quiver = nx.MultiDiGraph()
        quiver.add_nodes_from(vertices)
        for tail, head, _rank, _order, tag in sorted(self._entries):
            key = quiver.add_edge(tail, head)
            self._arrow[tag] = (tail, head, key)
        return quiver

    def arrow(self, tag):
        return self._arrow[tag]

    def tagged(self, rank):
        return {tag: arrow for tag, arrow in self._arrow.items() if tag[0] == rank}


def mutateAtVertex(quiver, relations, vertex):
    """Steps 1 to 7 at `vertex`, returning (new quiver, new relations).

    Does not reduce: like the paper, it leaves the cleanup after step 7 to the
    caller, which is `reduce` below.  Does not check admissibility either --
    `isMutable` is that, and applying the rewrite without it gives a quiver that
    is not derived equivalent (R-005).
    """
    outArrows = ap.arrowsOutOf(quiver, vertex)
    inArrows = ap.arrowsInto(quiver, vertex)
    outRelations = [r for r in relations if ap.source(r) == vertex]
    inRelations = [r for r in relations if ap.target(r) == vertex]
    throughRelations = [r for r in relations
                        if ap.source(r) != vertex and ap.target(r) != vertex]

    arrows = _Arrows()
    for order, arrow in enumerate(ap.arrowsOf(quiver)):
        tail, head, _key = arrow
        if tail != vertex and head != vertex:
            arrows.declare((_Arrows.CARRIED, arrow), tail, head, _Arrows.CARRIED, order)
    for order, alpha in enumerate(outArrows):                    # step 2: alpha*
        arrows.declare((_Arrows.FLIP, alpha), alpha[1], vertex, _Arrows.FLIP, order)
    for order, beta in enumerate(inArrows):                      # step 1: alpha beta
        for inner, alpha in enumerate(outArrows):
            arrows.declare((_Arrows.COMPOSITE, beta, alpha), beta[0], alpha[1],
                           _Arrows.COMPOSITE, (order, inner))
    for order, relation in enumerate(outRelations):              # step 3: rbar
        arrows.declare((_Arrows.RBAR, order), vertex, ap.target(relation),
                       _Arrows.RBAR, order)
    newQuiver = arrows.build(quiver.nodes)

    composite = lambda beta, alpha: arrows.arrow((_Arrows.COMPOSITE, beta, alpha))
    flip = lambda alpha: arrows.arrow((_Arrows.FLIP, alpha))
    carried = lambda arrow: arrows.arrow((_Arrows.CARRIED, arrow))
    rbar = lambda index: arrows.arrow((_Arrows.RBAR, index))

    newRelations = []

    # Step 4: one relation per arrow into the vertex.
    for beta in inArrows:
        newRelations.append(ap.combination(
            [(composite(beta, alpha), flip(alpha)) for alpha in outArrows]))

    # Step 5: rbar alpha* = r / alpha.
    for index, relation in enumerate(outRelations):
        for alpha in outArrows:
            through = ap.combination([((flip(alpha), rbar(index)), 1)])
            quotient = _carry(ap.leftDivide(relation, alpha), carried)
            newRelations.append(ap.add(through, ap.negate(quotient)))

    # Step 6: a relation into the vertex extends along each arrow out of it.
    for relation in inRelations:
        for alpha in outArrows:
            newRelations.append(ap.combination(
                (_carryPath(path[:-1], carried) + (composite(path[-1], alpha),), coefficient)
                for path, coefficient in relation.items()))

    # A relation past the vertex keeps its paths, with `beta` then `alpha`
    # replaced by the one composite arrow.
    for relation in throughRelations:
        newRelations.append(ap.combination(
            (_spliceOutVertex(path, vertex, carried, composite), coefficient)
            for path, coefficient in relation.items()))

    # Step 7: the relations out of the mutated vertex.
    newRelations.extend(
        _relationsOutOfMutatedVertex(quiver, relations, vertex, newQuiver,
                                     outArrows, outRelations, arrows))

    return newQuiver, normalise(newRelations)


def _carryPath(path, carried):
    """An old path whose arrows all avoid the mutated vertex, under its new names."""
    return tuple(carried(arrow) for arrow in path)


def _carry(comb, carried):
    return ap.combination((_carryPath(path, carried), coefficient)
                          for path, coefficient in comb.items())


def _spliceOutVertex(path, vertex, carried, composite):
    """`path` with the two arrows it runs through `vertex` along made one.

    `beta` into the vertex followed by `alpha` out of it is the single composite
    arrow `alpha beta` of step 1 -- which is a *different* arrow from any that
    joined those endpoints before, and the whole reason the vertex model could
    not do this: after deleting the vertex from the sequence the two read alike.
    """
    for position in range(len(path) - 1):
        if path[position][1] == vertex and path[position + 1][0] == vertex:
            return (_carryPath(path[:position], carried)
                    + (composite(path[position], path[position + 1]),)
                    + _carryPath(path[position + 2:], carried))
    return _carryPath(path, carried)


def _relationsOutOfMutatedVertex(quiver, relations, vertex, newQuiver,
                                 outArrows, outRelations, arrows):
    """Step 7, as a kernel, target by target.

    For a target `v`, the candidate paths out of `i*` are `rbar` followed by a
    path `k_r ~~> v` in the new quiver that does not run back through `i*`.  A
    combination of them is a relation iff every `alpha` out of `i` sends it into
    the old ideal, so the relations are the kernel of the map

        eps  |-->  ( sum_P eps_P (r_P / alpha) s_P  mod I )_alpha

    and only the part of that kernel not already forced by a relation to a
    nearer target is a new generator.

    Every arrow out of `i*` is an `rbar`, since the flips point into it and
    nothing else touches it, so a candidate's first arrow always names the
    relation it came from.  In the vertex model two relations `i ~~> k` gave two
    arrows the sequence could not tell apart and the step was skipped for that
    target outright.
    """
    relationOfArrow = {arrows.arrow((_Arrows.RBAR, index)): relation
                       for index, relation in enumerate(outRelations)}
    oldArrow = _oldArrowReading(arrows)
    found = []
    for target in _targetsInOrder(newQuiver, vertex):
        candidates = _pathsOutOfMutatedVertex(newQuiver, vertex, target)
        if not candidates:
            continue
        shadow = {}
        for path in candidates:
            relation = relationOfArrow.get(path[0])
            if relation is None:
                raise ValueError("a path out of the mutated vertex does not begin "
                                 "with a relation arrow: {0}".format(path))
            tail = _inOldQuiver(path[1:], oldArrow)
            for alpha in outArrows:
                quotient = ap.leftDivide(relation, alpha)
                shadow[(path, alpha)] = ap.postCompose(quotient, tail) if quotient else {}
        for element in _kernelOverIdeal(quiver, relations, candidates, outArrows, shadow):
            if not _isForcedByNearer(newQuiver, found, element, vertex):
                found.append(element)
    return found


def _oldArrowReading(arrows):
    """Each new arrow as the old path it stands for, where it stands for one.

    A carried arrow is the old arrow itself; a composite `alpha beta` is the old
    two-arrow path `beta` then `alpha`.  A flip and an `rbar` stand for no path
    of the old quiver, and neither can occur in a tail out of `i*` that does not
    return to `i*`.
    """
    reading = {}
    for tag, arrow in arrows.tagged(_Arrows.CARRIED).items():
        reading[arrow] = (tag[1],)
    for tag, arrow in arrows.tagged(_Arrows.COMPOSITE).items():
        reading[arrow] = (tag[1], tag[2])
    return reading


def _inOldQuiver(path, oldArrow):
    """A path of the mutated quiver read back as a path of the old one."""
    read = ()
    for arrow in path:
        if arrow not in oldArrow:
            raise ValueError("{0} has no reading in the old quiver".format((arrow,)))
        read = read + oldArrow[arrow]
    return read


def _targetsInOrder(newQuiver, vertex):
    """Every vertex reachable from the mutated vertex, nearest first."""
    lengths = nx.single_source_shortest_path_length(newQuiver, vertex)
    return [v for v, _ in sorted(lengths.items(), key=lambda kv: (kv[1], kv[0]))
            if v != vertex]


def _pathsOutOfMutatedVertex(newQuiver, vertex, target):
    """Paths from the mutated vertex to `target` that leave and do not return."""
    return [p for p in ap.allPathsBetween(newQuiver, vertex, target)
            if p and vertex not in ap.pathVertices(p)[1:]]


def _kernelOverIdeal(quiver, relations, candidates, outArrows, shadow):
    """The eps with sum_P eps_P shadow[P, alpha] in the ideal, for every alpha.

    Solved as a linear system over the rationals: reduce each shadow against a
    basis of the ideal between its endpoints, and the residues are the rows.
    """
    rows = {}
    for alpha in outArrows:
        residues = {}
        for path in candidates:
            element = shadow.get((path, alpha)) or {}
            residues[path] = _reduceAgainstIdeal(quiver, relations, element)
        basis = sorted({key for r in residues.values() for key in r})
        for key in basis:
            rows[(alpha, key)] = {path: residues[path].get(key, Fraction(0))
                                  for path in candidates}
    if not rows:
        # Every candidate is sent to zero, so every combination of them is a
        # relation; the single paths generate all of that.
        return [ap.combination([(path, 1)]) for path in candidates]
    return [
        ap.combination((path, _asInteger(coefficient))
                       for path, coefficient in solution.items() if coefficient)
        for solution in _nullSpace(candidates, rows)
    ]


def _reduceAgainstIdeal(quiver, relations, element):
    """`element` reduced modulo the ideal, as a dict of residual coefficients."""
    if not element:
        return {}
    pivots = ap.idealBasis(quiver, relations, ap.source(element), ap.target(element))
    return ap.reduceAgainstPivots(element, pivots)


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
    target = ap.target(element)
    spanning = []
    for relation in found:
        end = ap.target(relation)
        if end == target:
            spanning.append(relation)
            continue
        for after in ap.allPathsBetween(newQuiver, end, target):
            if after and vertex not in ap.pathVertices(after)[1:]:
                spanning.append(ap.postCompose(relation, after))
    if not spanning:
        return False
    pivots = ra._rowReduce([dict(s) for s in spanning])
    return not ap.reduceAgainstPivots(element, pivots)


# -- the cleanup after step 7 ---------------------------------------------

def reduce(quiver, relations):
    """The paper's "Note" after step 7, then minimality, both exactly.

    The Note: a relation containing a path of length one says that arrow equals
    a combination of longer paths, so the arrow and the relation both go, and
    every other relation has that arrow substituted out.  On combinations this
    is arithmetic rather than list surgery -- `c alpha + rest = 0` means
    `alpha = -rest / c`, and splicing that into a path is concatenation -- which
    is the part the set-of-paths model could only approximate.  Naming the arrow
    also makes *which* arrow leaves the quiver unambiguous: `remove_edge` on a
    pair of endpoints drops an arbitrary one of a parallel pair.

    Minimality: a relation is dropped when it lies in the ideal generated by the
    others, decided by linear algebra over that ideal.
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
    """One pass of the Note, or None when no relation has a path of one arrow."""
    for index, relation in enumerate(relations):
        singles = sorted(p for p in relation if len(p) == 1)
        if not singles:
            continue
        arrow = singles[0][0]
        coefficient = relation[singles[0]]
        rest = {p: c for p, c in relation.items() if p != singles[0]}
        replacement = ap.scale(rest, Fraction(-1, 1) / coefficient) if rest else {}

        newQuiver = nx.MultiDiGraph(quiver)
        if newQuiver.has_edge(*arrow):
            newQuiver.remove_edge(*arrow)
        newRelations = []
        for other in relations[:index] + relations[index + 1:]:
            rewritten = ap.substituteArrow(other, arrow, replacement)
            if rewritten:
                newRelations.append(rewritten)
        return newQuiver, newRelations
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
        if not ap.isInIdeal(quiver, others, relation):
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
    """Every arrow and every path reversed, coefficients and arrow names kept.

    An arrow `(tail, head, key)` becomes `(head, tail, key)`, which is injective,
    so a parallel pair stays a parallel pair and every relation can be read back.
    """
    reversed_ = nx.MultiDiGraph()
    reversed_.add_nodes_from(quiver.nodes)
    for tail, head, key in quiver.edges(keys = True):
        reversed_.add_edge(head, tail, key = key)
    return reversed_, [
        ap.combination((tuple((head, tail, key) for tail, head, key in reversed(path)),
                        coefficient)
                       for path, coefficient in relation.items())
        for relation in relations
    ]


def mutateLeftAtVertex(quiver, relations, vertex):
    """Left mutation: the same procedure on the opposite algebra, read back."""
    return _opposite(*mutateAtVertex(*_opposite(quiver, relations), vertex))
