"""Reorienting a relation-free tree, which the procedure does for nothing.

For a quiver with no relations whose underlying graph is a tree, mutating at a
**source** reverses exactly the arrows at that vertex and leaves the algebra
relation-free: it is the BGP reflection, performed by the mutation procedure.
Mutating left at a **sink** is its inverse.  The repo has always used the
consequence -- that the underlying undirected tree of a relation-free quiver
identifies its class -- on the strength of the *derived* statement, that two tree
algebras are derived equivalent exactly when the trees are isomorphic.  What this
module adds is that the same thing holds one step further down: the orientations
of a tree are **mutation** equivalent, with an explicit sequence, so a merge made
on a shared hereditary form is a merge of mutation classes and not only of
derived ones.  Research F-036.

**The sequence.**  To turn one orientation into another, take an edge `e` where
they differ, and a side `A` of the tree minus `e` that contains no other
differing edge -- one always exists, since a smallest such side cannot contain
one.  Now flip **every vertex of `A` exactly once**.  Each edge inside `A` is
reversed twice, so it comes back to where it was; `e` has one endpoint in `A`, so
it is reversed once.  Nothing else moves.

The order is what makes every flip legal.  If `e` points *out* of `A`, take a
topological order of the quiver induced on `A` and mutate right along it: the
first vertex is a source of `A`, and externally only `e` leaves `A`, so it is a
source of the whole quiver; and when the walk reaches a later vertex every arrow
that came into it has already been reversed, so it is a source in its turn.  If
`e` points *into* `A`, the mirror: reverse topological order, left mutation at
each, every one a sink.  Repeat until no edge differs, at most `n` mutations per
differing edge.

**Right mutations alone are enough**, which matters because the search walks only
those: flipping sources and never sinks still reaches every orientation, since
flipping every vertex once along a topological order reverses the whole quiver
and twice brings it back.  Sinks only shorten the path.  What they do not do is
make it short: the worst distance between two orientations of a tree runs 4, 6,
9, 12, 16 at orders 4 to 8 with right mutations, 3, 3, 6, 6, 10 with both, and a
classification search runs at depth 6.  That is the whole practical point of
having the sequence -- it is a path no search of the pipeline would ever find.

Vertex labels do not move under any of this: the mutation procedure here does not
renumber at all, so the sequence can be computed on the arrows alone and then
handed to `quiverMutationAtVertices`.
"""

import itertools

import networkx as nx

from . import mutation
from . import pathAlgebra


def arrowsOf(pathAlg):
    """The quiver's arrows as a frozenset of (tail, head) pairs."""
    return frozenset(pathAlg.quiver.edges())


def underlyingTree(pathAlg):
    """The underlying undirected graph, or None if it is not a tree.

    Parallel arrows and loops disqualify it as well, since the reflection
    argument is about a tree's edges being reversed one at a time.
    """
    graph = nx.Graph()
    graph.add_nodes_from(pathAlg.quiver.nodes)
    graph.add_edges_from((tail, head) for tail, head in pathAlg.quiver.edges())
    if graph.number_of_edges() != pathAlg.quiver.number_of_edges():
        return None
    if graph.number_of_nodes() == 0 or not nx.is_tree(graph):
        return None
    return graph


def isReflectable(pathAlg):
    """Whether this is a relation-free quiver on a tree, where all of this holds."""
    return not pathAlg.rels and underlyingTree(pathAlg) is not None


def orientationsOfTree(graph):
    """Every orientation of an undirected tree, as frozensets of arrows.

    There are 2^(n-1) of them and F-036 says they are one mutation class, so this
    is the whole class of the hereditary algebra, written out.
    """
    edges = sorted(tuple(sorted(edge)) for edge in graph.edges)
    for bits in itertools.product((0, 1), repeat = len(edges)):
        yield frozenset((first, second) if bit else (second, first)
                        for (first, second), bit in zip(edges, bits))


def quiverWithArrows(vertices, arrows):
    """A relation-free path algebra on the given vertices and arrows."""
    algebra = pathAlgebra.PathAlgebra()
    algebra.add_vertices_from(sorted(vertices))
    for tail, head in sorted(arrows):
        algebra.add_arrow(tail, head)
    return algebra


def reflectionSequence(arrows, targetArrows):
    """A mutation sequence taking one orientation of a tree to another.

    Both arguments are sets of arrows on the same vertices with the same
    underlying tree.  The result is in the repo's signed convention -- a positive
    vertex is a right mutation there, a negative one a left mutation -- and every
    step is at a source or a sink respectively, which is what makes it a sequence
    of reflections rather than a general walk.

    Raises when the two are not orientations of one tree, since then the argument
    does not apply and a sequence would mean nothing.
    """
    current = {tuple(sorted(arrow)): arrow for arrow in arrows}
    target = {tuple(sorted(arrow)): arrow for arrow in targetArrows}
    if set(current) != set(target):
        raise ValueError("the two orientations are not of the same graph")
    graph = nx.Graph()
    graph.add_edges_from(current)
    if not nx.is_tree(graph):
        raise ValueError("reflections are only free on a tree")

    sequence = []
    while True:
        differing = [edge for edge in current if current[edge] != target[edge]]
        if not differing:
            return sequence
        edge, side = _smallestSideWithoutAnother(graph, differing)
        tail, head = current[edge]
        pointsOutOfSide = tail in side
        order = _flipOrder(current, side, rightMutation = pointsOutOfSide)
        for vertex in order:
            sequence.append(vertex if pointsOutOfSide else -vertex)
            _flipAt(current, vertex)


def _smallestSideWithoutAnother(graph, differing):
    """A differing edge and a side of it holding no other differing edge.

    Taking the smallest side over every differing edge gives one: a smaller
    differing edge inside it would have a side smaller still.
    """
    best = None
    for edge in differing:
        cut = nx.Graph(graph)
        cut.remove_edge(*edge)
        for component in nx.connected_components(cut):
            if best is None or len(component) < len(best[1]):
                best = (edge, component)
    edge, side = best
    return edge, side


def _flipOrder(current, side, rightMutation):
    """The order to flip the vertices of `side` in, so every flip is legal.

    A topological order of the quiver induced on the side for right mutation at
    sources, and its reverse for left mutation at sinks.
    """
    induced = nx.DiGraph()
    induced.add_nodes_from(sorted(side))
    for edge, arrow in current.items():
        if edge[0] in side and edge[1] in side:
            induced.add_edge(*arrow)
    order = list(nx.topological_sort(induced))
    return order if rightMutation else list(reversed(order))


def _flipAt(current, vertex):
    """Reverse every arrow at a vertex, in the edge-to-arrow map."""
    for edge, (tail, head) in list(current.items()):
        if vertex in edge:
            current[edge] = (head, tail)


def reorient(pathAlg, targetArrows, check = True):
    """Mutate `pathAlg` into the given orientation, and return it.

    With `check`, every step is put to `mutationIsPossibleAtVertex` first and the
    result is compared with what was asked for -- which is how the sequence is
    verified against the engine rather than trusted.
    """
    if not isReflectable(pathAlg):
        raise ValueError("reflections are free only on a relation-free tree")
    sequence = reflectionSequence(arrowsOf(pathAlg), frozenset(targetArrows))
    if not check:
        return mutation.quiverMutationAtVertices(pathAlg, sequence), sequence
    walked = pathAlg
    for step in sequence:
        vertex = abs(step)
        # Left mutation is the procedure on the opposite algebra, so that is
        # where its admissibility is decided.
        legal = (walked.canMutateAt(vertex) if step > 0
                 else pathAlgebra.dualPathAlgebra(walked).canMutateAt(vertex))
        if not legal:
            raise ValueError("the procedure refuses a mutation at {0}".format(step))
        walked = mutation.quiverMutationAtVertices(walked, [step])
        if walked.rels:
            raise ValueError("a reflection left relations behind at {0}".format(vertex))
    if arrowsOf(walked) != frozenset(targetArrows):
        raise ValueError("the sequence did not reach the orientation asked for")
    return walked, sequence


def linearOrientation(graph):
    """The arrows of the linearly oriented line on a path graph, or None.

    A tree that is a path has an orientation whose algebra is an **LNA** -- the
    line with no relations -- so a search that reaches any orientation of a path
    has reached `kA_n`, with a sequence to prove it.
    """
    if graph.number_of_nodes() and not nx.is_tree(graph):
        return None
    degrees = dict(graph.degree())
    if max(degrees.values(), default = 0) > 2:
        return None
    ends = [vertex for vertex, degree in degrees.items() if degree <= 1]
    if not ends:
        return None
    walk = nx.shortest_path(graph, min(ends), max(ends))
    return frozenset(zip(walk, walk[1:]))


def linesReachedThroughReflections(pathAlg, depth, depthAfter = None,
                                   orientationLimit = None):
    """The LNAs a search reaches when it may reorient a relation-free quiver.

    A plain search walks right mutations only, and a reorientation is a *known*
    sequence, so a search that hits a relation-free quiver need not spend depth
    looking for one: it can jump to any orientation of the tree and carry on.
    This does that -- search to `depth`, then from every orientation of every
    relation-free quiver reached, search again to `depthAfter` (the same depth by
    default).

    `orientationLimit` caps how many orientations are tried per quiver, which
    matters because there are 2^(n-1) of them.  The ones tried are the first in
    the enumeration, which is a fixed order and not a random sample.

    Returns a dict from relation string to the shortest path found, the same
    shape `search.linesReachedFrom` returns, so the two are comparable.
    """
    from . import search

    reached = dict(search.linesReachedFrom(pathAlg, depth))
    for quiver in relationFreeQuiversReached(pathAlg, depth):
        tree = underlyingTree(quiver)
        orientations = orientationsOfTree(tree)
        if orientationLimit is not None:
            orientations = itertools.islice(orientations, orientationLimit)
        for arrows in orientations:
            if arrows == arrowsOf(quiver):
                continue
            start = quiverWithArrows(tree.nodes, arrows)
            for relationString, path in search.linesReachedFrom(
                    start, depth if depthAfter is None else depthAfter).items():
                if relationString not in reached:
                    reached[relationString] = path
    return reached


def relationFreeQuiversReached(pathAlg, depth, alsoDual = True):
    """Every relation-free tree quiver a bounded search out of `pathAlg` reaches.

    One per orientation actually reached, deduplicated by its arrows, which is
    what `linesReachedThroughReflections` then reorients.

    What the search of the **opposite** algebra reaches is carried back through
    the opposite before it is recorded, so every quiver here is one `pathAlg`
    itself reaches: a right-mutation path out of the opposite is a left-mutation
    path out of `pathAlg`, and it lands on the opposite of what that search saw.
    Recording the quiver as the dual search found it would put the reversed
    orientation in the list, which is a different algebra and the wrong answer to
    "which orientations does this one reach".
    """
    import copy

    from . import search

    found = {}

    def visit(quiver, mutationVertices, dualised = False):
        if quiver.rels or underlyingTree(quiver) is None:
            return
        if dualised:
            quiver = pathAlgebra.dualPathAlgebra(quiver)
        found.setdefault(arrowsOf(quiver), copy.deepcopy(quiver))

    search.mutationSearchDepthFirst(copy.deepcopy(pathAlg), depth, [], 'reflections',
                                    printOutput = False, visitor = visit)
    if alsoDual:
        search.mutationSearchDepthFirst(
            pathAlgebra.dualPathAlgebra(pathAlg), depth, [], 'reflections',
            printOutput = False,
            visitor = lambda quiver, path: visit(quiver, path, dualised = True))
    return list(found.values())


# -- joining two algebras through the hereditary quivers they reach --------
#
# `classification.mergeReport` merges two classes when both reach the same tree.
# That the two are *derived* equivalent needs only the tree; that they are in one
# *mutation* class -- which is what a class in this repo is -- needs the two
# orientations joined, which is what this module supplies.  `mutationBridge`
# makes the merge constructive: it returns the sequence, and the caller can
# apply it and look at where it lands.


def mutationBridge(first, second, depth, alsoDual = True):
    """A mutation sequence from `first` to an algebra isomorphic to `second`.

    Both are searched for a relation-free quiver; if the two quivers have
    isomorphic underlying trees, the sequence is

        `first`'s path to its tree, a reflection sequence onto `second`'s
        orientation, and `second`'s path run backwards,

    which is a mutation path because every step is a mutation and because a right
    mutation at a vertex is undone by a left one there.  Returns a dict with the
    three legs, the whole `sequence`, and the `algebra` it actually produces, or
    None when the two searches reach no common tree at this depth -- which is not
    a proof that none exists.

    An empty `reflect` leg means the two searches happened to land on the same
    orientation and the merge needed nothing from this module; a non-empty one is
    a merge that did.

    The result is `second` up to relabelling, not on the nose: the second leg is
    read through an isomorphism of the two trees, so the vertices arrive with
    `first`'s numbering.
    """
    from . import quipuForms
    from . import search

    leftSide = _hereditaryReach(first, depth, alsoDual)
    rightSide = _hereditaryReach(second, depth, alsoDual)
    for tree, (quiverA, pathA) in leftSide.items():
        if tree not in rightSide:
            continue
        quiverB, pathB = rightSide[tree]
        mapping = _treeIsomorphism(underlyingTree(quiverB), underlyingTree(quiverA))
        if mapping is None:
            continue
        target = frozenset((mapping[tail], mapping[head])
                           for tail, head in arrowsOf(quiverB))
        _reached, reflectSteps = reorient(quiverA, target)
        # Undo `second`'s path: backwards, each step the other way round, and
        # each vertex read through the isomorphism of the two trees.
        backwards = [(-mapping[abs(step)] if step > 0 else mapping[abs(step)])
                     for step in reversed(pathB)]
        sequence = list(pathA) + list(reflectSteps) + backwards
        return {
            'tree': tree,
            'toTree': list(pathA),
            'reflect': list(reflectSteps),
            'back': backwards,
            'sequence': sequence,
            'algebra': mutation.quiverMutationAtVertices(first, sequence),
        }
    return None


def _hereditaryReach(pathAlg, depth, alsoDual):
    """One relation-free quiver per tree reached, with the path that got there.

    The dual is searched as well, and what it reaches is carried back through the
    dual, so the path returned is a path out of `pathAlg` itself: a right
    mutation on the opposite algebra is a left mutation here.
    """
    import copy

    from . import quipuForms
    from . import search

    found = {}

    def collect(sign):
        def visit(quiver, mutationVertices):
            if quiver.rels or underlyingTree(quiver) is None:
                return
            if sign < 0:
                quiver = pathAlgebra.dualPathAlgebra(quiver)
            tree = quipuForms.canonicalTreeForm(underlyingTree(quiver))
            path = [sign * step for step in mutationVertices]
            if tree not in found or len(path) < len(found[tree][1]):
                found[tree] = (copy.deepcopy(quiver), path)
        return visit

    search.mutationSearchDepthFirst(copy.deepcopy(pathAlg), depth, [], 'bridge',
                                    printOutput = False, visitor = collect(1))
    if alsoDual:
        search.mutationSearchDepthFirst(
            pathAlgebra.dualPathAlgebra(pathAlg), depth, [], 'bridge',
            printOutput = False, visitor = collect(-1))
    return found


def _treeIsomorphism(fromTree, toTree):
    """One isomorphism of undirected trees as a vertex map, or None."""
    matcher = nx.algorithms.isomorphism.GraphMatcher(fromTree, toTree)
    for mapping in matcher.isomorphisms_iter():
        return dict(mapping)
    return None
