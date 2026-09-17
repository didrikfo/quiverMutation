"""Walking the mutation graph.

`mutationSearchDepthFirst` is the one search primitive: descend through
admissible right mutations to a bounded depth, recording every quiver reached
that is again a line, and every quiver reached with no relations left.  The
second kind identifies the derived equivalence class completely, since for
hereditary algebras of tree type the underlying undirected tree is the whole of
the class -- which is what `hereditaryFormsReachedFrom` collects.

Reachability here is one-way: the search only walks right mutations, so A can
reach B at a depth where B reaches nothing.  Searching from the relation dual as
well is what covers the other direction -- see `memberAndItsDual`.
"""

import contextlib
import copy

import networkx as nx

from . import lines
from . import mutation
from . import nakayama
from . import pathAlgebra
from . import paths
from . import quipuForms
from . import reduction


def mutationSearchDepthFirst(pathAlg, depth, mutationVertices = None, quiverName = 'quiver', vertexRelabeling = None, printOutput = True, collected = None, collectedHereditary = None, visitor = None):
    """Walk mutations of pathAlg to the given depth, recording the lines found.

    Every quiver reached that is again a line is recorded as a triple
    (path algebra, mutation path, vertex numbering).  Pass a list as `collected`
    to receive those triples in memory, in the order the search visits them.

    Pass a list as `collectedHereditary` to also receive, for every quiver
    reached that has no relations left, a triple (canonical form of the
    underlying undirected graph, the quipu notation for it where it applies,
    mutation path).  Those are the hereditary algebras in the class, and they
    identify it completely.

    Pass `visitor` to see *every* quiver the search reaches, line or not: it is
    called as `visitor(pathAlg, mutationVertices)` at each node, before the
    node's children are walked.  The two collectors above are the two questions
    asked often enough to have been built in; a visitor is for the rest, such as
    asking which quivers of a given shape a class passes through.

    `quiverName` only labels the progress output.  It used to name a
    '<quiverName>DF.txt' transcript that the caller parsed back by string
    slicing; that round trip is gone -- see NOTES.md idea 11.
    """
    # These used to default to [] and {}, which Python evaluates once at
    # definition time.  The relabeling dict is filled in below and so leaked
    # between searches: a second search in the same process inherited the
    # first one's numbering, and crashed as soon as the quiver was longer.
    mutationVertices = [] if mutationVertices is None else mutationVertices
    vertexRelabeling = {} if vertexRelabeling is None else dict(vertexRelabeling)
    vertices = list(pathAlg.vertices())
    baseQuiver = copy.deepcopy(pathAlg.quiver)
    quiverAtThisDepth = copy.deepcopy(pathAlg.quiver)
    rels = copy.deepcopy(pathAlg.rels)
    relsAtThisDepth = copy.deepcopy(pathAlg.rels)
    if not bool(vertexRelabeling):
        for vertex in vertices:
            vertexRelabeling[vertex] = vertex
    longestPathLength = 0
    noCycles = not bool(list(nx.simple_cycles(baseQuiver)))
    if noCycles:
        longestPathLength = nx.dag_longest_path_length(baseQuiver)
    if printOutput:
        print('Quiver name: ', quiverName)
        print("Mutations: ", mutationVertices)
        print('Numbering: {0}'.format(vertexRelabeling))
        print("Longest path: ", longestPathLength)
        pathAlgebra.printPathAlgebra(pathAlg)
    if not bool(rels) and (collectedHereditary is not None or _SIGHTING_SINKS):
        # No relations left: the algebra is hereditary, and the underlying
        # undirected graph of its quiver is a complete derived invariant.
        graph = quipuForms.underlyingGraph(pathAlg)
        if collectedHereditary is not None:
            collectedHereditary.append((
                quipuForms.canonicalUndirectedForm(graph),
                quipuForms.formatQuipu(quipuForms.quipuParameters(graph)),
                mutationVertices[:],
            ))
        for sink in _SIGHTING_SINKS:
            sink.append(describeRelationFreeQuiver(pathAlg, mutationVertices, graph))
    if visitor is not None:
        visitor(pathAlg, mutationVertices)
    isLine = (longestPathLength == len(vertices) - 1) and (len(baseQuiver.edges) == len(vertices) - 1)
    if isLine and collected is not None:
        foundPathAlg = pathAlgebra.PathAlgebra()
        foundPathAlg.quiver = baseQuiver
        foundPathAlg.rels = rels
        collected.append((foundPathAlg, mutationVertices[:], dict(vertexRelabeling)))
    # debugVertexList = [1, 1, 2, 1, 2, 3, 5, 3, 4, 4, 5, 2, 2, 3, 6, 1, 4, 1, 2, 3, 1, 4]
    # for i in range(7, len(debugVertexList)):
    #      if mutationVertices == debugVertexList[:i]:
    #          input('Press enter to continue...')
    if depth > 0 and noCycles:
        depth = depth - 1
        for vertex in reversed(vertices):
            discardMutation = False
            pathAlg.quiver = copy.deepcopy(quiverAtThisDepth)
            pathAlg.rels = copy.deepcopy(relsAtThisDepth)
            mutationVerticesAtDepth = mutationVertices[:]
            # `mutationIsPossibleAtVertex` is the whole gate.  It used to be
            # followed here by a second test of the same idea, counting the
            # paths into the vertex against the paths through each arrow out of
            # it and refusing the vertex if any one arrow lost a path.  That is
            # the *strict* reading -- every arrow must keep every path -- where
            # the paper's theorem rules mutation out only when a path dies
            # against all of them, so the search was refusing mutations the
            # paper allows, twice over.  F-016.
            if mutation.mutationIsPossibleAtVertex(pathAlg, vertex):
                mutationVerticesAtDepth.append(vertexRelabeling[vertex])
                mutPathAlg = mutation.quiverMutationAtVertex(pathAlg, vertex)
                for rel in mutPathAlg.rels:
                    if paths.isIllegalRelation(mutPathAlg, rel):
                        discardMutation = True
                        break
                if discardMutation:
                    break
                mutPathAlg = reduction.reducePathAlgebra(mutPathAlg)
                mutationSearchDepthFirst(copy.deepcopy(mutPathAlg), depth, mutationVerticesAtDepth, quiverName, vertexRelabeling, printOutput, collected, collectedHereditary, visitor)
    return


def hereditaryFormsReachedFrom(pathAlg, depth):
    """The hereditary algebras reachable from pathAlg within `depth` mutations.

    Returns a dict mapping the canonical form of the underlying undirected graph
    to (quipu notation, shortest mutation path found to it).  An empty result
    means no relation-free quiver was reached at this depth, not that none
    exists.
    """
    found = []
    mutationSearchDepthFirst(pathAlg, depth, [], 'hereditary', printOutput = False,
                             collected = None, collectedHereditary = found)
    forms = {}
    for canonical, quipu, path in found:
        if canonical not in forms or len(path) < len(forms[canonical][1]):
            forms[canonical] = (quipu, path)
    return forms


def findHereditaryFormForClass(table, lineLength, className, maxDepth = 8, printOutput = False):
    """Search the members of one class for a relation-free quiver.

    Iterative deepening from each member in turn, returning as soon as any
    member reaches one.  Members are tried shortest-relation-string first, on
    the observation that an LNA with fewer and shorter relations tends to need
    fewer mutations to shed them all.

    Each member is searched from itself and from its relation dual.  The search
    only walks right mutations, so reachability is one-way; reversing every arrow
    is one of the class-preserving operations of arXiv:2305.06642, and
    rightMutate(dual(P)) = dual(leftMutate(P)), so a right-mutation path out of
    the dual is a left-mutation path out of the member and everything it reaches
    is still in the class.  The dual of a tree quiver is the same tree, so a form
    reached from the dual is the class' form unchanged.

    Returns the hereditary form as it should be written into the table -- the
    quipu notation where the graph is a quipu, otherwise the canonical tree
    encoding -- or '' if nothing was reached.
    """
    members = sorted(table.membersOfClass(className), key = lambda r: (len(r), r))
    for depth in range(2, maxDepth + 1):
        for relationString in members:
            for startPoint in memberAndItsDual(lineLength, relationString):
                forms = hereditaryFormsReachedFrom(startPoint, depth)
                if forms:
                    if printOutput:
                        print('class {0} reaches {1} at depth {2} from {3!r}'.format(
                            className, sorted(forms), depth, relationString))
                    return formatHereditaryForms(forms)
    return ''


def memberAndItsDual(lineLength, relationString):
    """An LNA and its relation dual, both as path algebras, without repeats."""
    algebra = nakayama.LinearNakayamaAlgebra.fromRelationString(lineLength, relationString)
    dual = algebra.relationDual()
    return [algebra] if dual == algebra else [algebra, dual]


def formatHereditaryForms(forms):
    """Render the result of hereditaryFormsReachedFrom for the table.

    All the hereditary algebras in one derived equivalence class have isomorphic
    underlying graphs, so this is normally one value.  More than one would mean
    either a bug in the mutation procedure or a non-tree in the mix, so they are
    all reported, joined by '|', rather than silently reduced to one.
    """
    return '|'.join(sorted(quipu or canonical for canonical, (quipu, _path) in forms.items()))


def hereditaryFormFromTheorem(lineLength, relationString):
    """The hereditary form of one LNA, straight from the quipu theorem.

    Returns '' when the LNA does not have almost separate relations, since
    theorem `thm:QuipuToAn` of arXiv:2305.06642 says nothing about those.  For
    the ones it does cover this is O(1), where reaching the same answer by
    mutation search costs a depth-4-to-9 traversal.
    """
    algebra = nakayama.LinearNakayamaAlgebra.fromRelationString(lineLength, relationString)
    return algebra.quipuName()


def linesReachedFrom(pathAlg, depth, alsoDual = True):
    """The LNAs a bounded mutation search out of `pathAlg` reaches.

    Returns a dict from relation string to the shortest mutation path found to
    it.  This is what settles a Coxeter-polynomial lead: the polynomial says two
    algebras *could* be derived equivalent, and a mutation path from one to the
    other says they are, since every step of the procedure is a tilting
    mutation.

    With `alsoDual`, the search is run from the opposite algebra as well, whose
    lines are read back through the dual.  Two algebras are derived equivalent
    exactly when their opposites are, and the opposite of an LNA is an LNA, so a
    line reached from the opposite is as good a witness as one reached directly
    -- and it is the only way to see what a *left* mutation path would reach,
    since the search walks right mutations only.
    """
    startPoints = [pathAlg]
    if alsoDual:
        startPoints.append(pathAlgebra.dualPathAlgebra(pathAlg))
    reached = {}
    for index, startPoint in enumerate(startPoints):
        collected = []
        mutationSearchDepthFirst(copy.deepcopy(startPoint), depth, [], 'lines',
                                 printOutput = False, collected = collected)
        for found in lines.mutationListLineCleanup(collected, printOutput = False):
            algebra = found[0] if index == 0 else pathAlgebra.dualPathAlgebra(found[0])
            # Sorted, because a relation set read off the dual comes out in the
            # reverse order and `relSetToString` writes it down as it stands --
            # which would make one LNA look like two.
            relationString = lines.relSetToString(sorted(_renumberedLine(algebra).rels))
            path = found[1]
            if relationString not in reached or len(path) < len(reached[relationString]):
                reached[relationString] = path
    return reached


def _renumberedLine(pathAlg):
    """A line algebra with its vertices renumbered 1 -> ... -> n along the line.

    `mutationListLineCleanup` already does this for what it collects; taking the
    opposite afterwards reverses the numbering, so it has to be done again.
    """
    order = nx.topological_sort(pathAlg.quiver)
    relabeling = {vertex: position for position, vertex in enumerate(order, start = 1)}
    renumbered = pathAlgebra.PathAlgebra()
    renumbered.add_vertices_from(sorted(relabeling.values()))
    for tail, head in pathAlg.quiver.edges():
        renumbered.add_arrow(relabeling[tail], relabeling[head])
    for rel in pathAlg.rels:
        renumbered.add_rel([[relabeling[vertex] for vertex in path] for path in rel])
    return renumbered


# -- relation-free sightings ---------------------------------------------
#
# Relations are what the procedure spends its time on, so a mutation that leaves
# a quiver with *none* is an event: the algebra is hereditary, and the underlying
# graph of the quiver settles its derived equivalence class outright.  The
# searches already use that -- `hereditaryFormsReachedFrom` is nothing else --
# but they use it locally and throw the rest away, which means nobody has ever
# looked at what those quivers are.  The expectation is that every one is a tree
# and almost every one a quipu; a relation-free quiver whose graph has a *cycle*
# would be a hereditary algebra of a kind no LNA class has produced, and would be
# worth stopping for.  Neither is known, and the cost of finding out is a few
# lines, since every search already passes through the place where it could be
# recorded.
#
# Sightings are opt-in: with no sink open, the branch above does no more work
# than it did before.

_SIGHTING_SINKS = []


@contextlib.contextmanager
def relationFreeSightings():
    """Record every relation-free quiver the searches inside the block reach.

    Yields the list the sightings accumulate in, one dict per sighting, in the
    order the searches visit them.  Sinks nest, so an inner block does not stop
    an outer one from seeing what it sees.

        with search.relationFreeSightings() as sightings:
            classification.classifyLength(9)
        print(search.summariseSightings(sightings))
    """
    sink = []
    _SIGHTING_SINKS.append(sink)
    try:
        yield sink
    finally:
        _SIGHTING_SINKS.remove(sink)


def describeRelationFreeQuiver(pathAlg, mutationVertices = None, graph = None):
    """What a relation-free quiver is, in the terms worth counting.

    `isTree` and `isQuipu` are about the *underlying undirected* graph, since
    that is what decides a hereditary algebra's derived equivalence class.
    `hasOrientedCycle` is about the quiver, and `parallelArrows` counts arrows
    the underlying simple graph merges -- either of those would mean the algebra
    is not the path algebra of a tree, whatever the undirected picture says.
    """
    graph = quipuForms.underlyingGraph(pathAlg) if graph is None else graph
    isTree = nx.is_tree(graph) if graph.number_of_nodes() else False
    return {
        'canonical': quipuForms.canonicalUndirectedForm(graph),
        'quipu': quipuForms.formatQuipu(quipuForms.quipuParameters(graph)),
        'vertices': graph.number_of_nodes(),
        'arrows': pathAlg.quiver.number_of_edges(),
        'isTree': isTree,
        'isQuipu': bool(isTree and quipuForms.isQuipuByDegrees(graph)),
        'isConnected': bool(graph.number_of_nodes()) and nx.is_connected(graph),
        'hasOrientedCycle': bool(list(nx.simple_cycles(pathAlg.quiver))),
        'parallelArrows': pathAlg.quiver.number_of_edges() - graph.number_of_edges(),
        'maxDegree': max((degree for _, degree in graph.degree()), default = 0),
        'path': list(mutationVertices or []),
    }


def summariseSightings(sightings):
    """Counts of the kinds of relation-free quiver a run saw.

    `quipus`, `otherTrees` and `notTrees` partition the sightings; `distinct` is
    how many isomorphism classes of underlying graph they came to, which is the
    number that says whether a run saw one thing many times or many things.
    `oddities` is the sightings that are not trees, kept in full, since those are
    the ones there is no reason to expect.
    """
    oddities = [sighting for sighting in sightings if not sighting['isTree']]
    return {
        'sightings': len(sightings),
        'distinct': len({sighting['canonical'] for sighting in sightings}),
        'quipus': sum(1 for sighting in sightings if sighting['isQuipu']),
        'otherTrees': sum(1 for sighting in sightings
                          if sighting['isTree'] and not sighting['isQuipu']),
        'notTrees': len(oddities),
        'oddities': oddities,
    }
