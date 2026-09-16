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

import copy

import networkx as nx

from . import mutation
from . import nakayama
from . import pathAlgebra
from . import paths
from . import quipuForms
from . import reduction


def mutationSearchDepthFirst(pathAlg, depth, mutationVertices = None, quiverName = 'quiver', vertexRelabeling = None, printOutput = True, collected = None, collectedHereditary = None):
    """Walk mutations of pathAlg to the given depth, recording the lines found.

    Every quiver reached that is again a line is recorded as a triple
    (path algebra, mutation path, vertex numbering).  Pass a list as `collected`
    to receive those triples in memory, in the order the search visits them.

    Pass a list as `collectedHereditary` to also receive, for every quiver
    reached that has no relations left, a triple (canonical form of the
    underlying undirected graph, the quipu notation for it where it applies,
    mutation path).  Those are the hereditary algebras in the class, and they
    identify it completely.

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
    if collectedHereditary is not None and not bool(rels):
        # No relations left: the algebra is hereditary, and the underlying
        # undirected graph of its quiver is a complete derived invariant.
        graph = quipuForms.underlyingGraph(pathAlg)
        collectedHereditary.append((
            quipuForms.canonicalUndirectedForm(graph),
            quipuForms.formatQuipu(quipuForms.quipuParameters(graph)),
            mutationVertices[:],
        ))
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
                mutationSearchDepthFirst(copy.deepcopy(mutPathAlg), depth, mutationVerticesAtDepth, quiverName, vertexRelabeling, printOutput, collected, collectedHereditary)
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
