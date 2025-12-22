import copy
import networkx as nx
import pathAlgebraClass
from quiver_mutation_io import printPathAlgebra
from quiver_mutation_relations import allMinimalRelsBetweenVertices, allRelsBetweenVertices, isSubRelOf, reducePathAlgebra, zeroizeRels
from quiver_mutation_utils import listIntersection, getVertexNumberingKeyFromValue
from quiver_mutation_algebra import coxeterPoly, cartanMatrix
from quiver_mutation_io import plotQuiver

def quiverMutationAtVertex(pathAlg, vertex):
    oldQuiver = pathAlg.quiver
    vertices = pathAlg.vertices()
    mutPathAlg = pathAlgebraClass.PathAlgebra()
    mutPathAlg.add_vertices_from(vertices)
    vertexOutArrows = pathAlg.out_arrows(vertex)
    vertexOutRels = pathAlg.out_rels(vertex)
    allMinOutRels = []
    for rel in vertexOutRels:
        for minOutRel in allMinimalRelsBetweenVertices(pathAlg, vertex, rel[0][-1]):
            allMinOutRels.append((minOutRel, vertexOutRels.index(rel)))
    arrowTargetVertices = list(pathAlg.quiver.successors(vertex))
    vertexSuccessors = list(nx.dfs_preorder_nodes(oldQuiver, vertex))
    vertexSuccessors.remove(vertex)
    for ar in vertexOutArrows:
        if ar[1] in vertexSuccessors:
            vertexSuccessors.remove(ar[1])
    for v in vertexSuccessors:
        predecessors = list(nx.dfs_preorder_nodes(nx.reverse(oldQuiver), v))
        targetPredecessorsOfVertex = listIntersection(predecessors, arrowTargetVertices)
        if bool(targetPredecessorsOfVertex):
            outRelPathsWithRels = []
            for rel in vertexOutRels:
                outRelPathsWithRels.append([rel, []])
            targetRelAndPathHitByOutRelAndPath = []
            for w in targetPredecessorsOfVertex:
#                nonMinOutRels = nonMinimalOutRels(pathAlg, w)
                nonMinOutRels = []
                for w_succ in list(nx.dfs_preorder_nodes(oldQuiver, w))[1:]:
                    nonMinOutRels.extend(allRelsBetweenVertices(pathAlg, w, w_succ))
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
                    if isSubRelOf(rule7Rel, rel):
                        legitRule7Relation = False
                        break
                    elif isSubRelOf(rel, rule7Rel):
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
    zeroizedRels = zeroizeRels(mutPathAlg.rels)
    mutPathAlg.rels = zeroizedRels
    return mutPathAlg

def quiverMutationAtVertices(pathAlg, vertices : list, printMutationSteps = False):
    for i in range(len(vertices)):
        vertex = vertices[i]
        if vertex > 0:
            pathAlg = quiverMutationAtVertex(pathAlg, vertex)
        else:
            pathAlg = leftQuiverMutationAtVertex(pathAlg, -vertex)
        pathAlg = reducePathAlgebra(pathAlg)
        if printMutationSteps:
            print('Mutations: ', vertices[:i + 1])
            printPathAlgebra(pathAlg)
    return pathAlg

def dualPathAlgebra( pathAlg ):
    dualPathAlg = pathAlgebraClass.PathAlgebra()
    dualPathAlg.add_vertices_from(pathAlg.vertices())
    for arrow in pathAlg.arrows():
        dualPathAlg.add_arrow(arrow[1], arrow[0])
    for rel in pathAlg.rels:
        dualRel = []
        for relPath in rel:
            dualRel.append(list(reversed(relPath)))
        dualPathAlg.add_rel(dualRel)
    return dualPathAlg

def leftQuiverMutationAtVertex(pathAlg, vertex):
    dualPathAlg = dualPathAlgebra(pathAlg)
    mutDualPathAlg = quiverMutationAtVertex(dualPathAlg, vertex)
    leftMutPathAlg = dualPathAlgebra(mutDualPathAlg)
    return leftMutPathAlg

def reverseMutationSequence(mutationVertices, vertexNumbering):
    reverseMutationVertices = []
    for n in range(len(mutationVertices) - 1, -1, -1):
        reverseMutationVertices.append(-getVertexNumberingKeyFromValue(vertexNumbering, mutationVertices[n]))
    return reverseMutationVertices

def reverseMutationFromSequence(pathAlg, mutationVertices, vertexNumbering):
    reverseMutationVertices = reverseMutationSequence(mutationVertices, vertexNumbering)
    return quiverMutationAtVertices(pathAlg, reverseMutationVertices)

def quiverMutation(pathAlgebra, mutationVertexList, firstDisplayedStep = 0):
    #
    baseCoxPol = coxeterPoly(pathAlgebra)
    print(coxeterPoly(pathAlgebra))
    printPathAlgebra(pathAlgebra)
    if firstDisplayedStep == 0:
        plotQuiver(pathAlgebra)
    for i in range(len(mutationVertexList)):
        if mutationVertexList[i] >= 0:
            pathAlgebra = quiverMutationAtVertex(pathAlgebra, mutationVertexList[i])
        else:
            pathAlgebra = leftQuiverMutationAtVertex(pathAlgebra, -mutationVertexList[i])
        pathAlgebra = reducePathAlgebra(pathAlgebra)
        currentCoxPol = coxeterPoly(pathAlgebra)
        cartMat = cartanMatrix(pathAlgebra)
        print(np.matrix(cartMat))
        if currentCoxPol != baseCoxPol:
            print('COXETER POLYNOMIAL HAS CHANGED!')
        print('Mutations: ', mutationVertexList[0:i + 1])
        print(currentCoxPol)
        printPathAlgebra(pathAlgebra)
        if i + 1 >= firstDisplayedStep:
            plotQuiver(pathAlgebra)
    return pathAlgebra

def onePointExtension(pathAlgebra, arrowToAdd, relsToAdd = []):
    extendedPathAlg = pathAlgebra
    extendedPathAlg.add_arrows_from([arrowToAdd])
    extendedPathAlg.add_rels_from(relsToAdd)
    return extendedPathAlg
