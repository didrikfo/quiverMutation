"""The tilting mutation procedure itself.

`quiverMutationAtVertex` is steps 1-7 of the combinatorial rule of
arXiv:2112.08129, `mutationIsPossibleAtVertex` its admissibility condition, and
`quiverMutationAtVertices` the two composed with `reduction.reducePathAlgebra`
after each step.  Left mutation is right mutation conjugated by the opposite
algebra.

A mutation applied outside the admissibility condition still returns a quiver,
just not a derived equivalent one -- see research R-005 -- so callers that care
about the class must test admissibility first.
"""

import copy

import networkx as nx
import numpy as np

from . import invariants
from . import pathAlgebra
from . import paths
from . import plotting
from . import reduction



def quiverMutationAtVertex(pathAlg, vertex):
    oldQuiver = pathAlg.quiver
    vertices = pathAlg.vertices()
    mutPathAlg = pathAlgebra.PathAlgebra()
    mutPathAlg.add_vertices_from(vertices)
    vertexOutArrows = pathAlg.out_arrows(vertex)
    vertexOutRels = pathAlg.out_rels(vertex)
    allMinOutRels = []
    for rel in vertexOutRels:
        for minOutRel in paths.allMinimalRelsBetweenVertices(pathAlg, vertex, rel[0][-1]):
            allMinOutRels.append((minOutRel, vertexOutRels.index(rel)))
    arrowTargetVertices = list(pathAlg.quiver.successors(vertex))
    vertexSuccessors = list(nx.dfs_preorder_nodes(oldQuiver, vertex))
    vertexSuccessors.remove(vertex)
    for ar in vertexOutArrows:
        if ar[1] in vertexSuccessors:
            vertexSuccessors.remove(ar[1])
    for v in vertexSuccessors:
        predecessors = list(nx.dfs_preorder_nodes(nx.reverse(oldQuiver), v))
        targetPredecessorsOfVertex = paths.listIntersection(predecessors, arrowTargetVertices)
        if bool(targetPredecessorsOfVertex):
            outRelPathsWithRels = []
            for rel in vertexOutRels:
                outRelPathsWithRels.append([rel, []])
            targetRelAndPathHitByOutRelAndPath = []
            for w in targetPredecessorsOfVertex:
                nonMinOutRels = []
                for w_succ in list(nx.dfs_preorder_nodes(oldQuiver, w))[1:]:
                    nonMinOutRels.extend(paths.allRelsBetweenVertices(pathAlg, w, w_succ))
#                nonMinOutRelsZeroized = zeroizeRels(nonMinOutRels)
                for targetOutRel in nonMinOutRels:
                    if targetOutRel[0][-1] == v:
                        equivMinOutRelsUsedForTargetOutRel = []
                        for targetOutRelPath in targetOutRel:
                            for outRel in allMinOutRels:
                                for outRelPath in outRel[0]:
                                    rule7RelPath = [vertex] + targetOutRelPath[len(outRelPath) - 2:]
                                    if [vertex] + targetOutRelPath[:len(outRelPath) - 1] == outRelPath and not (rule7RelPath, outRel[1]) in equivMinOutRelsUsedForTargetOutRel:
                                        #rule7RelPath = [vertex] + targetOutRelPath[len(outRelPath) - 2:]
                                        targetRelAndPathHitByOutRelAndPath.append((targetOutRel, targetOutRelPath, outRel[0], outRelPath, rule7RelPath))
                                        equivMinOutRelsUsedForTargetOutRel.append((rule7RelPath, outRel[1]))
                                        break
            relevantTargetOutRels = []
            for rule7Tuple in targetRelAndPathHitByOutRelAndPath:
                if not rule7Tuple[0] in relevantTargetOutRels:
                    relevantTargetOutRels.append(rule7Tuple[0])
            possibleRule7RelCandidates = []
            for targetOutRel in relevantTargetOutRels:
                rule7TuplesForTargetOutRel = []
                for rule7Tuple in targetRelAndPathHitByOutRelAndPath:
                    if rule7Tuple[0] == targetOutRel:
                        rule7TuplesForTargetOutRel.append(rule7Tuple)
                potentialRule7RelCandidates = [[]]
                mutationPossibleAtRel = False
                for targetOutRelPath in targetOutRel:
                    mutationPossibleAtRel = True
                    newCandidates = []
                    for rule7Tuple in rule7TuplesForTargetOutRel:
                        if rule7Tuple[1] == targetOutRelPath:
                            for candidate in potentialRule7RelCandidates:
                                newCandidates.append(candidate + [rule7Tuple])
                    if not bool(newCandidates):
                        mutationPossibleAtRel = False
                        break
                    potentialRule7RelCandidates = newCandidates
                if mutationPossibleAtRel:
                    possibleRule7RelCandidates.extend(potentialRule7RelCandidates)
            doneTargetOutRels = []
            rule7RelCandidates = []
            for possibleCandidate in possibleRule7RelCandidates:
                skipCandidate = False
                relevantOutRels = []
                for rel in doneTargetOutRels:
                    if rel == possibleCandidate[0]:
                        skipCandidate = True
                if not skipCandidate:
                    for rule7Tuple in possibleCandidate:
                        if not rule7Tuple[2] in relevantOutRels:
                            relevantOutRels.append(rule7Tuple[2])
                    missingOutRelPathTuples = []
                    for outRel in relevantOutRels:
                        outRelPathsInCandidate = []
                        for outRelPath in outRel:
                            for rule7Tuple in possibleCandidate:
                                if rule7Tuple[3] == outRelPath and not rule7Tuple[3] in outRelPathsInCandidate:
                                    outRelPathsInCandidate.append(rule7Tuple[3])
                                    break
                        for relPath in outRel:
                            if not relPath in outRelPathsInCandidate:
                                missingOutRelPathTuples.append((outRel, relPath))
                    if not bool(missingOutRelPathTuples):
                        rule7RelCandidates.append(possibleCandidate)
                        doneTargetOutRels.append(possibleCandidate[0][0])
                    else:
                        mutationPossible = True
                        for missingPath in missingOutRelPathTuples:
                            outRelPathHasRel = False
                            for rule7Tuple in targetRelAndPathHitByOutRelAndPath:
                                if (rule7Tuple[2], rule7Tuple[3]) == missingPath and len(rule7Tuple[0]) == 1:
                                    outRelPathHasRel = True
                                    break
                            if not outRelPathHasRel:
                                mutationPossible = False
                        if mutationPossible:
                            rule7RelCandidates.append(possibleCandidate)
                            doneTargetOutRels.append(possibleCandidate[0][0])
            for candidate in rule7RelCandidates:
                rule7Rel = []
                for rule7Tuple in candidate:
                    if not rule7Tuple[4] in rule7Rel:
                        rule7Rel.append(rule7Tuple[4])
                legitRule7Relation = True
                removeRel = []
                for rel in mutPathAlg.rels:
                    if paths.isSubRelOf(rule7Rel, rel):
                        legitRule7Relation = False
                        break
                    elif paths.isSubRelOf(rel, rule7Rel):
                        removeRel = rel
                        break
                if bool(removeRel):
                    mutPathAlg.rels.remove(rel)
                if legitRule7Relation:
                    mutPathAlg.add_rel(rule7Rel)
    for ar in oldQuiver.edges:
        if ar[0] == vertex:
            mutPathAlg.add_arrow(ar[1], ar[0])
#            for rel in vertexOutRels:
            possibleRelsToAdd = []
            for rel in allMinOutRels:
                arrowTargetVertexNotInRel = True
                for relPath in rel[0]:
                    if relPath[1] == ar[1]:
                        arrowTargetVertexNotInRel = False
                        newRel = []
                        newRel.append([ar[1], ar[0], relPath[-1]])
                        for otherRelPath in rel[0]:
                            if otherRelPath[1] == ar[1]:
                                newRel.append(otherRelPath[1:])
                        possibleRelsToAdd.append((newRel, rel[1]))
#                        mutPathAlg.add_rel(newRel)
                if arrowTargetVertexNotInRel:
                    newRel = [[ar[1], ar[0], rel[0][0][-1]]]
                    possibleRelsToAdd.append((newRel, rel[1]))
#                    mutPathAlg.add_rel(newRel)
            possibleRelsToAdd.sort(reverse=True, key=lambda x:len(x[0]))
            relsToAdd = []
            for newRelNumber in range(len(vertexOutRels)):
                newRelIndex = [y[1] for y in possibleRelsToAdd].index(newRelNumber)
                relsToAdd.append(possibleRelsToAdd[newRelIndex][0])
            for relToAdd in relsToAdd:
                mutPathAlg.add_rel(relToAdd)
        elif ar[1] == vertex:
            newRel = []
            for arOut in vertexOutArrows:
                mutPathAlg.add_arrow(ar[0], arOut[1])
                newRel.append([arOut[1], vertex])
            for relPath in newRel:
                relPath.insert(0, ar[0])
            mutPathAlg.add_rel(newRel)
        elif (ar[0] != vertex) & (ar[1] != vertex):
            mutPathAlg.add_arrow(ar[0], ar[1])
    # for rel in allMinOutRels:
    #     mutPathAlg.add_arrow(rel[0][0], rel[0][-1])
    for rel in copy.deepcopy(pathAlg.rels):
        if rel[0][0] == vertex:
            mutPathAlg.add_arrow(rel[0][0], rel[0][-1])
        if rel[0][-1] == vertex:
            for v in arrowTargetVertices:
                newRel = []
                for relPath in rel:
                    newRelPath = relPath[:-1]
                    newRelPath.append(v)
                    newRel.append(newRelPath)
                mutPathAlg.add_rel(newRel)
        elif (rel[0][0] != vertex) and (rel[0][-1] != vertex):
            for relPath in rel:
                if vertex in relPath[:]:
                    relPath.remove(vertex)
            mutPathAlg.add_rel(rel)
    zeroizedRels = paths.zeroizeRels(mutPathAlg.rels)
    mutPathAlg.rels = zeroizedRels
    return mutPathAlg


def quiverMutationAtVertices(pathAlg, vertices : list, printMutationSteps = False):
    for i in range(len(vertices)):
        vertex = vertices[i]
        if vertex > 0:
            pathAlg = quiverMutationAtVertex(pathAlg, vertex)
        else:
            pathAlg = leftQuiverMutationAtVertex(pathAlg, -vertex)
        pathAlg = reduction.reducePathAlgebra(pathAlg)
        if printMutationSteps:
            print('Mutations: ', vertices[:i + 1])
            pathAlgebra.printPathAlgebra(pathAlg)
    return pathAlg


def leftQuiverMutationAtVertex(pathAlg, vertex):
    dualPathAlg = pathAlgebra.dualPathAlgebra(pathAlg)
    mutDualPathAlg = quiverMutationAtVertex(dualPathAlg, vertex)
    leftMutPathAlg = pathAlgebra.dualPathAlgebra(mutDualPathAlg)
    return leftMutPathAlg


def mutationIsPossibleAtVertex(pathAlg, vertex, allRels = None):
    """Whether the mutation procedure may be applied to pathAlg at vertex.

    This is the admissibility test of theorem 1 in arXiv:2112.08129, as the
    depth-first search has always applied it:

    * There must be an arrow out of vertex.  P_i* is the cocone of a right
      approximation of P_i by the other indecomposable projectives, so with no
      arrow out of i there is nothing to approximate by.
    * The quiver must have no pair of parallel arrows.  This is a restriction of
      this implementation rather than of the procedure: a relation is modelled
      as a set of vertex sequences, which cannot distinguish two arrows with
      the same source and target.
    * Hom(P_i*[1], Lambda) = 0, which holds iff every nonzero path ending in i
      composes nonzero with at least one arrow out of i.  A minimal zero
      relation whose last arrow starts in i, and whose truncation by that last
      arrow is itself nonzero, is a witness that it fails.

    Note that the last test is stricter than the paper's condition when vertex
    has more than one arrow out of it: it rejects the vertex as soon as one
    arrow out of it kills a nonzero path, where the paper only requires that
    some arrow out of it does not.  The two agree whenever vertex has a single
    arrow out of it, which is the only case the linear Nakayama search meets.
    """
    if not bool(pathAlg.out_arrows(vertex)):
        return False
    for ar in pathAlg.arrows():
        if ar[2] > 0:
            return False
    if allRels is None:
        allRels = paths.allRelsInPathAlgebra(pathAlg)
    for rel in allRels:
        if len(rel) == 1 and rel[0][-2] == vertex and (not [rel[0][:-1]] in allRels):
            return False
    return True


def showMutationSteps(pathAlg, mutationVertexList, firstDisplayedStep = 0):
    """Mutate step by step, printing and plotting each one, for looking by hand.

    `quiverMutationAtVertices` is the same walk without the commentary.  This
    also flags a step where the Coxeter polynomial moves, which it must not do
    along an admissible path -- if it does, the sequence left the admissible
    region (research R-005).

    It was called `quiverMutation`, which now names the package.
    """
    baseCoxPol = invariants.coxeterPoly(pathAlg)
    print(invariants.coxeterPoly(pathAlg))
    pathAlgebra.printPathAlgebra(pathAlg)
    if firstDisplayedStep == 0:
        plotting.plotQuiver(pathAlg)
    for i in range(len(mutationVertexList)):
        if mutationVertexList[i] >= 0:
            pathAlg = quiverMutationAtVertex(pathAlg, mutationVertexList[i])
        else:
            pathAlg = leftQuiverMutationAtVertex(pathAlg, -mutationVertexList[i])
        pathAlg = reduction.reducePathAlgebra(pathAlg)
        currentCoxPol = invariants.coxeterPoly(pathAlg)
        cartMat = invariants.cartanMatrix(pathAlg)
        print(np.matrix(cartMat))
        if currentCoxPol != baseCoxPol:
            print('COXETER POLYNOMIAL HAS CHANGED!')
        print('Mutations: ', mutationVertexList[0:i + 1])
        print(currentCoxPol)
        pathAlgebra.printPathAlgebra(pathAlg)
        if i + 1 >= firstDisplayedStep:
            plotting.plotQuiver(pathAlg)
    return pathAlg


def reverseMutationSequence(mutationVertices, vertexNumbering):
    reverseMutationVertices = []
    for n in range(len(mutationVertices) - 1, -1, -1):
        reverseMutationVertices.append(-getVertexNumberingKeyFromValue(vertexNumbering, mutationVertices[n]))
    return reverseMutationVertices


def getVertexNumberingKeyFromValue(vertexNumbering, vertex):
    keyList = list(vertexNumbering.keys())
    valueList = list(vertexNumbering.values())
    if vertex > 0:
        position = valueList.index(vertex)
        key = keyList[position]
    else:
        position = valueList.index(-vertex)
        key = -keyList[position]
    return key
