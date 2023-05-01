import copy
import glob
import itertools
import ast
import matplotlib.pyplot as plt
import math
import networkx as nx
import os
import pathAlgebraClass
import sympy
import time
import csv
import numpy as np
from datetime import datetime
from sympy.interactive.printing import init_printing
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
from sympy.abc import x, y

now = datetime.now() # current date and time


def listIntersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

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

def nonMinimalOutRels(pathAlg, vertex):
    nonMinOutRels = []
    pathAlgCopy = copy.deepcopy(pathAlg)
    for rel in pathAlgCopy.out_rels(vertex):
        nonMinOutRels.extend(extendRel(pathAlgCopy, rel))
    for ar in pathAlgCopy.out_arrows(vertex):
        deeperOutRels = nonMinimalOutRels(pathAlgCopy, ar[1])
        for dRel in deeperOutRels:
            extendedRel = []
            for dRelPath in dRel:
                extendedRelPath = [vertex] + dRelPath
                extendedRel.append(extendedRelPath)
            nonMinOutRels.append(extendedRel)
        for rel in pathAlgCopy.out_rels(vertex):
            if len(rel) == 1 and rel[0][1] == ar[1]:
                for dRel in deeperOutRels:
                    targetInRel = False
                    for dRelPath in dRel:
                        if rel[0][-1] in dRelPath:
                            targetInRel = True
                            break
                    if targetInRel:
                        for dRelPath in dRel:
                            nonMinOutRels.append([[vertex] + dRelPath[:]])
            for relPath in rel:
                for dRel in copy.deepcopy(deeperOutRels):
                    extendedFirstRel = copy.deepcopy(rel)
                    extendedFirstRel.remove(relPath)
                    if len(rel) > 1 or len(dRel) > 1:
                        for dRelPath in dRel:
                            if dRelPath[:len(relPath)-1] == relPath[1:]:
                                extendedLastRel = dRel
                                extendedLastRel.remove(dRelPath)
                                for efRelPath in extendedFirstRel:
                                    efRelPath.extend(dRelPath[len(relPath) - 1:])
                                for elRelPath in extendedLastRel:
                                    elRelPath.insert(0, vertex)
                                nonMinOutRels.append(extendedFirstRel + extendedLastRel)
    uniqueNonMinOutRels = []
    for nonMinRel in nonMinOutRels:
        if not sorted(nonMinRel) in uniqueNonMinOutRels:
            uniqueNonMinOutRels.append(sorted(nonMinRel))
    return uniqueNonMinOutRels


def allRelsBetweenVertices(pathAlg, startVertex,endVertex):
    pathAlgCopy = copy.deepcopy(pathAlg)
    relsBetween = pathAlgCopy.rels_between(startVertex, endVertex)
    #allRelsBetween= pathAlgCopy.rels_between(startVertex, endVertex)
    for ar in pathAlgCopy.out_arrows(startVertex):
        dRelsBetween = allRelsBetweenVertices(pathAlg, ar[1], endVertex)
        for rel in dRelsBetween:
            newRel = []
            for relPath in rel:
                newRel.append([startVertex] + relPath)
            if not newRel in relsBetween:
                relsBetween.append(newRel)
    verticesBetween = []
    for path in nx.all_simple_paths(pathAlg.quiver, startVertex, endVertex):
        for i in path:
            if not i in verticesBetween:
                verticesBetween.append(i)
    for i in verticesBetween:
            if i != endVertex:
                for rel in pathAlgCopy.rels_between(startVertex, i):
                    for path in nx.all_simple_paths(pathAlg.quiver, i, endVertex):
                        newRel = []
                        for relPath in rel:
                            newRel.append(relPath + path[1:])
                        if not newRel in relsBetween:
                            relsBetween.append(newRel)
    for i in range(len(relsBetween) - 1):
        for j in range(i + 1, len(relsBetween)):
            if relsBetween[i] == relsBetween[j]:
                print('Duplicate rel: ', relsBetween[i])
    relSetsToApply = [[rel] for rel in relsBetween]
    allRelsBetween = relsBetween
    stillNewRels = True
    while stillNewRels:
        newRels = []
        stillNewRels = False
        for rel in allRelsBetween:
            for relPath in rel:
                for relSet in relSetsToApply:
                    newRelPaths = applyRelSetToPath(relPath, relSet)
                    relToAdd = sorted(rel[:rel.index(relPath)] + rel[rel.index(relPath) + 1:] + newRelPaths)
                    if not relToAdd in allRelsBetween and not any(relToAdd.count(x) > 1 for x in relToAdd) and relToAdd != []:
                        allRelsBetween.append(copy.deepcopy(relToAdd))
                        newRels.append(relToAdd)
                        stillNewRels = True
        relSetsToApply = [[copy.deepcopy(newRel)] for newRel in newRels]
    return allRelsBetween

def allMinimalRelsBetweenVertices(pathAlg, startVertex,endVertex):
    pathAlgCopy = copy.deepcopy(pathAlg)
    relsBetween = pathAlgCopy.rels_between(startVertex, endVertex)
    allRelsBetween= pathAlgCopy.rels_between(startVertex, endVertex)
    verticesBetween = []
    for path in nx.all_simple_paths(pathAlg.quiver, startVertex, endVertex):
        for i in path:
            if not i in verticesBetween:
                verticesBetween.append(i)
    intermideateRels = []
    for i in verticesBetween:
        for j in verticesBetween:
            if i != startVertex or j != endVertex:
                for rel in pathAlgCopy.rels_between(i, j):
                    if not rel in intermideateRels:
                        intermideateRels.append(rel)
    powerSetOfShorterRels = powerset(intermideateRels)
    relSetsToApply = []
    for relSet in powerSetOfShorterRels:
        if bool(relSet):
            relSetsToApply.append(relSet)
    for rel in relsBetween:
        newRelSetsToApply = []
        for relSet in relSetsToApply:
            for differentRel in relsBetween:
                if differentRel != rel:
                    newRelSetsToApply.append(relSet + [differentRel])
        relSetsToApply.extend(newRelSetsToApply)
        for relPath in rel:
            for relSet in relSetsToApply:
                newRelPaths = applyRelSetToPath(relPath, relSet)
                relToAdd = sorted(rel[:rel.index(relPath)] + rel[rel.index(relPath) + 1:] + newRelPaths)
                if not relToAdd in allRelsBetween and not any(relToAdd.count(x) > 1 for x in relToAdd) and relToAdd != []:
                    allRelsBetween.append(relToAdd)
    return allRelsBetween

def extendRel(pathAlg, rel):
    vertex = rel[0][-1]
    extendedRels = [rel]
    outArrows = pathAlg.out_arrows(vertex)
    for ar in outArrows:
        extendedRel = []
        for relPath in rel:
            extendedRelPath = relPath + [ar[1]]
            extendedRel.append(extendedRelPath)
        extendedRels.append(extendedRel)
        deeperExtendedRels = extendRel(pathAlg, extendedRel)
        extendedRels.extend(deeperExtendedRels)
    return extendedRels

def isSubRelOf(potentialSubRel, relation):
    isSubRel = True
    for relPath in potentialSubRel:
        if not relPath in relation:
            isSubRel = False
            break
    return isSubRel

def reducePathAlgebra(pathAlg):
    quiver = pathAlg.quiver
    redPathAlg = pathAlgebraClass.PathAlgebra()
    redPathAlg.add_vertices_from(quiver.nodes)
    redPathAlg.add_arrows_from(quiver.edges(keys=True))
    for rel in pathAlg.rels:
        if isIllegalRelation(pathAlg, rel):
            pathAlg.rels.remove(rel)
    redPathAlg.add_rels_from(pathAlg.rels)
    redPathAlg = removeDuplicateRelPaths(redPathAlg)
    noChange = False
    while not noChange:
        noChange = True
        rels = copy.deepcopy(redPathAlg.rels)
        redRels = copy.deepcopy(redPathAlg.rels)
        arrowStillExists = True
        removedRel = ([],[])
        for rel in rels:
            relIndex = redRels.index(rel)
            for relPath in rel:
                if len(relPath) == 2:
                    rel.remove(relPath)
                    removedRel = (relPath, rel)
                    redPathAlg.quiver.remove_edge(relPath[0], relPath[1])
                    redRels.remove(redRels[relIndex])
                    noChange = False
                    if redPathAlg.quiver.number_of_edges(relPath[0], relPath[1]) > redRels.count(rel):
                        arrowStillExists = True
                    else:
                        arrowStillExists = False
                    break
            if not noChange:
                break
        if not noChange and not arrowStillExists:
            relsToRemove = []
            newRels = []
            for rel in redRels:
                newRel = []
                unchangedPaths = []
                for relPath in rel:
                    newPaths = []
                    if bool(removedRel[0]):
                        if sublistExists(relPath, removedRel[0]):
                            for leftoverPath in removedRel[1]:
                                newPath = copy.deepcopy(relPath[:relPath.index(removedRel[0][1])] + leftoverPath[1:-1] + relPath[relPath.index(removedRel[0][1]):])
                                newPaths.append(newPath)
                            if not rel in relsToRemove:
                                relsToRemove.append(rel)
                    if bool(newPaths):
                        newRel.extend(newPaths)
                    else:
                        unchangedPaths.append(relPath[:])
                if bool(newRel):
                    newRel.extend(unchangedPaths)
                    newRels.append(newRel)
            for rel in relsToRemove:
                redRels.remove(rel)
            redRels.extend(newRels)
            while [] in redRels:
                redRels.remove([])
        redPathAlg.clear_rels()
        redPathAlg.add_rels_from(redRels)
        redPathAlg = removeDuplicateRelPaths(redPathAlg)
    redPathAlg = removeDuplicateRels(redPathAlg)
    redPathAlg = removeNonminimalZeroRels(redPathAlg, applyCommutativityRels=False)
    removeRedundantRelations(redPathAlg)
    removeExistingSubrelations(redPathAlg)
    noChange = False
    while not noChange:
        numberOfRelsBeforeRed = len(redPathAlg.rels)
        removeRedundantRelations(redPathAlg)
        if numberOfRelsBeforeRed == len(redPathAlg.rels):
            noChange = True
    redPathAlg = removeNonminimalZeroRels(redPathAlg)
    return redPathAlg


def removeNonminimalZeroRels(pathAlg, applyCommutativityRels = True):
    rels = pathAlg.rels
    minimalRels = []
    commutativityRels = []
    zeroRels = []
    for rel in rels:
        if len(rel) > 2:
            minimalRels.append(rel)
        elif len(rel) == 2:
            minimalRels.append(rel)
            commutativityRels.append(rel)
        else:
            zeroRels.append(rel)
    relSetsToApply = []
    if applyCommutativityRels:
        commutativityRelsPowerSet = powerset(commutativityRels)
        for relSet in commutativityRelsPowerSet:
            for perm in itertools.permutations(relSet):
                relSetsToApply.append(list(perm))
    for rel1 in zeroRels:
        removeRel = False
        for relSetToApply in relSetsToApply:
            if applyCommutativityRels:
                rel1Equivalent = applyRelSetToPath(rel1[0], relSetToApply)
            else:
                rel1Equivalent = rel1
            for rel2 in zeroRels:
                    if len(rel2[0]) < len(rel1Equivalent[0]):
                        for i in range(len(rel1Equivalent[0]) - len(rel2[0]) + 1):
                           if rel1Equivalent[0][i:i+len(rel2[0])] == rel2[0]:
                                removeRel = True
                                break
                    if removeRel:
                        break
            if removeRel:
                break
        if not removeRel:
            if not rel1 in minimalRels:
                minimalRels.append(rel1)
    pathAlg.rels = minimalRels
    return pathAlg

def reduceCommutativityRels(pathAlg):
    for rel in pathAlg.rels:
        if len(rel) > 1:
            relLength = len(rel[0])
            for i in range(1, len(rel)):
                relLength = max(relLength, len(rel[i]))
            sameStart = False
            sameEnd = False
            for i in range(1, relLength):
                for n in range(1,len(rel)):
                    if rel[n][i] != rel[0][i] and not sameStart:
                        sameStart = True
                        sameStartEnd = i
                    if rel[n][-i] != rel[0][-i] and not sameEnd:
                        sameEnd = True
                        sameEndStart = -i
                    if sameStart and sameEnd:
                        break
                if sameStart and sameEnd:
                    break
            if sameStart or sameEnd:
                newRel = []
                for relPath in rel:
                    newRel.append(relPath[sameStartEnd - 1:len(relPath) + sameEndStart + 2])
                pathAlg.rels[pathAlg.rels.index(rel)] = newRel
    return pathAlg

def removeDuplicateRels(pathAlg):
    uniqueRels = []
    for rel in pathAlg.rels:
        uniqueRel = []
        for relPath in rel:
            if not relPath in uniqueRel:
                uniqueRel.append(relPath)
        if bool(uniqueRel) and not uniqueRel in uniqueRels:
            uniqueRels.append(uniqueRel)
    pathAlg.rels = uniqueRels
    return pathAlg

def removeDuplicateRelPaths(pathAlg):
    uniqueRels = []
    for rel in pathAlg.rels:
        uniqueRelPaths = []
        for relPath in rel:
            if not relPath in uniqueRelPaths:
                uniqueRelPaths.append(relPath)
        uniqueRels.append(uniqueRelPaths)
    pathAlg.clear_rels()
    pathAlg.add_rels_from(uniqueRels)
    return pathAlg


def minimizeCommutingRelation(pathAlg, relation):
    if len(relation) >= 2:
        return pathAlg
    verticesBetween = []
    for path in nx.all_simple_paths(pathAlg.quiver, relation[0][0], relation[0][-1]):
        for vertex in path:
            if not vertex in verticesBetween:
                verticesBetween.append(vertex)
    rels = []
    for i in verticesBetween:
        for j in verticesBetween:
            for rel in allRelsBetweenVertices(pathAlg, i, j):
                rels.append(rel)
    replaceCommutingPartWith = []
    for rel in rels:
        if len(rel) >= 2 and len(rel[0]) < len(relation[0]):
            for m in range(len(rel)):
                for i in range(len(rel[m]) - 2):
                    commutativeRelStart = 0
                    if rel[m][i] in relation[0] and rel[m][i+1] not in relation[0]:
                        commutativeRelStart = rel[m][i]
                        if bool(commutativeRelStart):
                            commutativeRelEnd = 0
                            for j in range(len(rel[m]) - 1, i + 1, -1):
                                if rel[m][j] in relation[0] and rel[m][j-1] not in relation[0]:
                                    commutativeRelEnd = rel[m][j]
                                if bool(commutativeRelEnd):
                                    for relPath in rel[m+1:]:
                                        if relation[0][relation[0].index(commutativeRelStart):relation[0].index(commutativeRelEnd) + 1] == relPath[relPath.index(commutativeRelStart):relPath.index(commutativeRelEnd) + 1]:
                                            replaceCommutingPartWith = rel[m][i:j+1]
                                            break
                            if bool(replaceCommutingPartWith):
                                break
                if bool(replaceCommutingPartWith):
                    break
            if bool(replaceCommutingPartWith):
                break
    if bool(replaceCommutingPartWith):
        pathAlg.rels[pathAlg.rels.index(relation)][0][relation[0].index(replaceCommutingPartWith[0]):relation[0].index(replaceCommutingPartWith[-1])+1] = replaceCommutingPartWith
    return pathAlg

def removeRedundantRelations(pathAlg):
    #sort relations by increasing length of their longest path, then by increasing number of paths
    shortestRelsList = sorted(pathAlg.rels, key=lambda x: (len(max(x, key=len)), len(x)))
    necesarryRels = []
    while bool(shortestRelsList):
        necesarryRels.append(shortestRelsList.pop(0))
        relSetsToApply = powerset(necesarryRels)
        relsToRemove = []
        for rel in shortestRelsList:
            removeRel = False
            for relPath in rel:
                for relsToApply in relSetsToApply[1:]:
                    newPaths = applyRelSetToPath(relPath, relsToApply)
                    if sorted(newPaths) == sorted( rel[:rel.index(relPath)] + rel[rel.index(relPath) + 1:]):
                        relsToRemove.append(rel)
                        removeRel = True
                        break
                if removeRel:
                    break
        for rel in relsToRemove:
            shortestRelsList.remove(rel)
    pathAlg.rels = sorted(necesarryRels)
    return

def removeExistingSubrelations(pathAlg):
    reducedRels = []
    for rel in sorted(pathAlg.rels, key=len):
        subRelExists = False
        reducedRel = copy.deepcopy(rel)
        for redRel in reducedRels:
            for redRelPath in redRel:
                if redRelPath in rel:
                    subRelExists = True
                else:
                    subRelExists = False
                    break
            if subRelExists:
                for redRelPath in redRel:
                    reducedRel.remove(redRelPath)
        reducedRels.append(reducedRel)
    pathAlg.rels = reducedRels
    return

def zeroizeRels(rels):
    zeroRels = []
    nonZeroRels = []
    for rel in rels:
        if len(rel) == 1:
            zeroRels.append(rel)
        else:
            nonZeroRels.append(rel)
    zeroizedRels = zeroRels.copy()
    for rel in nonZeroRels:
        isZero = False
        for relPath in rel:
            isZeroPath = False
            for zeroRel in zeroRels:
                if sublistExists(relPath, zeroRel[0]):
                    rel.remove(relPath)
                    isZeroPath = True
                    break
        zeroizedRels.append(rel)
    zeroizedRelsReduced = [rel for rel in zeroizedRels if rel != []]
    return zeroizedRelsReduced


def sublistExists(list, sublist):
    for i in range(len(list)-len(sublist)+1):
        if sublist == list[i:i+len(sublist)]:
            return True #return position (i) if you wish
    return False

def isIllegalRelation(pathAlgebra, relation):
    isIllegal = False
    for n in range(1, len(relation)):
        if relation[n] == []:
            isIllegal = True
            break
        if relation[n][0] != relation[0][0] or relation[n][-1] != relation[0][-1]:
            isIllegal = True
            break
        if not nx.is_path(pathAlgebra.quiver, relation[n]):
            isIllegal = True
            break
        for i in range(0, len(relation[n]) - 1):
            for j in range(i + 1, len(relation[n])):
                if relation[n][i] == relation[n][j]:
                    isIllegal = True
                    break
        for m in range(n + 1, len(relation)):
            if relation[n] == relation[m]:
                isIllegal = True
                break
    if isIllegal:
        print('ILLEGAL RELATION!')
        print('The relation {0} '.format(relation))
        print('is illegal in the following path algebra:')
        printPathAlgebra(pathAlgebra)
        # input('Press enter to continue...')
    return isIllegal

def printPathAlgebra(pathAlg):
    print('Vertices: ', pathAlg.quiver.nodes)
    print('Arrows: ', pathAlg.quiver.edges)
    print('Relations: ', pathAlg.rels, '\n')
    return

def plotQuiver(pathAlg, showPlot = True, saveToFile = False, fileName = 'quiverPlot.png', folder = ''):
    try:
        os.mkdir(folder)
    except OSError:
        print("Creation of the directory %s failed" % folder)
    else:
        print("Successfully created the directory %s " % folder)
    savePath = folder + fileName
    nx.draw_networkx(pathAlg.quiver)
    if saveToFile:
        plt.savefig(savePath)
        plt.close()
    if showPlot:
        plt.show()
    return

def mutationSearchDepthFirst(pathAlg, depth, mutationVertices = [], quiverName = 'quiver', vertexRelabeling = {}, printOutput = True):
    vertices = list(pathAlg.vertices())
    baseQuiver = copy.deepcopy(pathAlg.quiver)
    quiverAtThisDepth = copy.deepcopy(pathAlg.quiver)
    rels = copy.deepcopy(pathAlg.rels)
    allRels = []
    for v in vertices:
        for w in vertices:
            allRels.extend(allRelsBetweenVertices(pathAlg, v, w))
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
        printPathAlgebra(pathAlg)
    fileName = '{0}DF.txt'.format(quiverName)
    with open(fileName, "a") as f:
        if (longestPathLength == len(vertices) - 1) and (len(baseQuiver.edges) == len(vertices) - 1):
            f.write('Mutations: {0}\n'.format(mutationVertices))
            f.write('Numbering: {0}\n'.format(vertexRelabeling))
            f.write("Longest path: {0}\n".format(longestPathLength))
            f.write('Vertices: {0}\n'.format(vertices))
            f.write('Arrows: {0}\n'.format(baseQuiver.edges))
            f.write('Relations: {0}\n'.format(rels))
            f.write('-\n')
        f.close()
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
            mutationPossible = True
            if not bool(pathAlg.out_arrows(vertex)):
                mutationPossible = False
            else:
                for ar in pathAlg.arrows():
                    if ar[2] > 0:
                        mutationPossible = False
                        break
            if mutationPossible:
#                for ar in pathAlg.out_arrows(vertex):
#                    for rel in pathAlg.in_rels(ar[1]):
#                for relStartVertex in vertices:
                for rel in allRels:
                    if len(rel) == 1 and rel[0][-2] == vertex and (not [rel[0][:-1]] in allRels):
                        mutationPossible = False
                        break
            if mutationPossible:
                vertexPredecessors = nx.dfs_preorder_nodes(nx.reverse(pathAlg.quiver), vertex)
                vertexImmideateSuccessors = list(pathAlg.quiver.successors(vertex))
                for v in vertexPredecessors:
                    numberOfPathsToVertexUpToRels = numberOfPathsUpToRels(pathAlg, v, vertex)
                    for w in vertexImmideateSuccessors:
                        if numberOfPathsToVertexUpToRels > numberOfPathsUpToRels(pathAlg, v, w):
                            mutationPossible = False
                            break
                    if not mutationPossible:
                        break
                # for v in vertexPredecessors:
                #     for nonMinRel in nonMinimalOutRels(pathAlg, v):
                #         if len(nonMinRel) == 1:
                #             if nonMinRel[0][-1] in vertexImmideateSuccessors and nonMinRel[0][-2] != vertex:
                #                 mutationPossible = False
                #                 break
                #     if not mutationPossible:
                #         break
            if mutationPossible:
                mutationVerticesAtDepth.append(vertexRelabeling[vertex])
                mutPathAlg = quiverMutationAtVertex(pathAlg, vertex)
                for rel in mutPathAlg.rels:
                    if isIllegalRelation(mutPathAlg, rel):
                        discardMutation = True
                        break
                if discardMutation:
                    break
                mutPathAlg = reducePathAlgebra(mutPathAlg)
                mutationSearchDepthFirst(copy.deepcopy(mutPathAlg), depth, mutationVerticesAtDepth, quiverName, vertexRelabeling, printOutput)
    return

def divisors(n):
    # get factors and their counts
    factors = {}
    nn = n
    i = 2
    while i*i <= nn:
        while nn % i == 0:
            factors[i] = factors.get(i, 0) + 1
            nn //= i
        i += 1
    if nn > 1:
        factors[nn] = factors.get(nn, 0) + 1
    primes = list(factors.keys())
    # generates factors from primes[k:] subset
    def generate(k):
        if k == len(primes):
            yield 1
        else:
            rest = generate(k+1)
            prime = primes[k]
            for factor in rest:
                prime_to_i = 1
                # prime_to_i iterates prime**i values, i being all possible exponents
                for _ in range(factors[prime] + 1):
                    yield factor * prime_to_i
                    prime_to_i *= prime
    # in python3, `yield from generate(0)` would also work
    for factor in generate(0):
        yield factor

def relabelLineAlgebra(pathAlg, currentRelabeling = {}):
    lineQuiver = pathAlg.quiver
    if not bool(currentRelabeling):
        for vertex in lineQuiver.nodes:
            currentRelabeling[vertex] = vertex
    isLineQuiver = False
    if len(lineQuiver.edges) == len(lineQuiver.nodes) - 1:
        if not bool(list(nx.simple_cycles(lineQuiver))):
            if nx.dag_longest_path_length(lineQuiver) == (len(lineQuiver.nodes) - 1):
                isLineQuiver = True
    if not isLineQuiver:
        return (lineQuiver, currentRelabeling)
    standardLineQuiver = nx.MultiDiGraph()
    standardLineQuiver.add_nodes_from(lineQuiver)
    for i in range(1, len(list(standardLineQuiver.nodes))):
        standardLineQuiver.add_edge(i, i+1)
    for vertex in lineQuiver.nodes:
        if not bool(lineQuiver.in_edges(vertex)):
            sourceVertex = vertex
            break
    vertexOrderChange = {1 : sourceVertex}
    for i in range(1, len(list(standardLineQuiver.nodes))):
        targetVertex = list(lineQuiver.out_edges(sourceVertex))[0][1]
        vertexOrderChange[i + 1] = targetVertex
        sourceVertex = targetVertex
    oldVerticesList = list(vertexOrderChange.values())
    newVerticesList = list(vertexOrderChange.keys())
    newLineRels = []
    for rel in pathAlg.rels:
        newRelPath = []
        for i in rel[0]:
            newRelPath.append(newVerticesList[oldVerticesList.index(i)])
        newLineRels.append([newRelPath])
    newLineRels.sort()
    newRelabeling = {}
    for i in range(1, len(list(standardLineQuiver.nodes)) + 1):
        newRelabeling[i] = currentRelabeling[vertexOrderChange[i]]
    pathAlg.quiver = standardLineQuiver
    pathAlg.rels = newLineRels
    return (pathAlg, newRelabeling)

def saveLinePathAlgMutation(pathAlg, mutationVertices = [], vertexRelabeling = {}, fileName = 'lineQuiver.txt'):
    quiver = pathAlg.quiver
    rels = pathAlg.rels
    vertices = quiver.nodes
    if len(quiver.edges) != len(vertices) - 1:
        return
    if bool(list(nx.simple_cycles(quiver))):
        return
    elif nx.dag_longest_path_length(quiver) != len(quiver.edges):
        return
    with open(fileName, "a") as f:
            f.write('Mutations: {0}\n'.format(mutationVertices))
            f.write('Numbering: {0}\n'.format(vertexRelabeling))
            f.write('Vertices: {0}\n'.format(quiver.nodes))
            f.write('Arrows: {0}\n'.format(quiver.edges))
            f.write('Relations: {0}\n'.format(rels))
            f.write('-\n')
    return

def readMutationsFromFile(fileName):
    rPathAlg = pathAlgebraClass.PathAlgebra()
    mutationList = []
    with open(fileName, 'r') as f:
        line = f.readline()
        while line != '':
            if line[0] == 'M':
                mutationVertices = []
                for v in list(line[12:-2].split(',')):
                    if v:
                        mutationVertices.append(int(v))
            elif line[0] == 'N':
                vertexRelabeling = {}
                for labelAsStr in list(line[12:-2].split(', ')):
                    label = list(labelAsStr.split(': '))
                    vertexRelabeling[int(label[0])] = int(label[1])
            elif line[0] == 'V':
                rQuiver = nx.MultiDiGraph()
                vertices = list(line[11:-2].split(','))
                for v in vertices:
                    rQuiver.add_node(int(v))
            elif line[0] == 'A':
                arrows = list(line[9:-2].split('), '))
                for arStr in arrows:
                    arStr = arStr.removeprefix('(')
                    arStr = arStr.removesuffix(')')
                    ar = tuple(map(int, arStr.split(', ')))
                    rQuiver.add_edge(ar[0], ar[1])
                rPathAlg.quiver = rQuiver
            elif line[0] == 'R':
                rRels = []
                if line[12] == ']':
                    relArrows = []
                else:
                    relArrows = list(line[12:-2].split(']], '))
                for relStr in relArrows:
                    if relStr[-1] != ']':
                        relStr = relStr + ']]'
                    relList = ast.literal_eval(relStr)
                    rRels.append(relList)
                rPathAlg.rels = copy.deepcopy(rRels)
            elif line[0] == '-':
                mutationList.append((copy.deepcopy(rPathAlg), mutationVertices, vertexRelabeling))
            line = f.readline()
    return mutationList

def mutationListLineCleanup(mutationList, relabelNodes = True, printOutput = True):
    modifiedList = []
    for mut in mutationList:
        quiv = copy.deepcopy(mut[0].quiver)
        vertices = quiv.nodes
        vertexRelabeling = mut[2]
        isLineQuiver = True
        if len(quiv.edges) != len(vertices) - 1:
            isLineQuiver = False
        if bool(list(nx.simple_cycles(quiv))):
            isLineQuiver = False
        elif nx.dag_longest_path_length(quiv) != len(quiv.edges):
            isLineQuiver = False
        if isLineQuiver:
            if relabelNodes:
                if printOutput:
                    print('Numbering: ', vertexRelabeling)
                    printPathAlgebra(mut[0])
                (pathAlg, newVertexRelabeling) = relabelLineAlgebra(mut[0], vertexRelabeling)
                if printOutput:
                    print('Renumbering: ', newVertexRelabeling)
                    printPathAlgebra(pathAlg)
            else:
                pathAlg = mut[0]
            rels = copy.deepcopy(pathAlg.rels)
            index = len(modifiedList)
            keepQuiver = True
            replaceQuiver = False
            for modMut in modifiedList:
                newRels = copy.deepcopy(rels)
                oldRels = list(modMut[0].rels)
                if newRels == oldRels:
                    if len(mut[1]) < len(modMut[1]):
                        replaceQuiver = True
                        index = modifiedList.index(modMut)
                        break
                    else:
                        keepQuiver = False
                        break
                elif newRels < oldRels:
                    index = modifiedList.index(modMut)
                    break
            if keepQuiver:
                if replaceQuiver:
                    modifiedList[index] = (pathAlg, mut[1], newVertexRelabeling)
                else:
                    modifiedList.insert(index, (pathAlg, mut[1], newVertexRelabeling))
    return modifiedList

def mutationListLineCleanupKeepDupes(mutationList, relabelNodes = True, discardLongerDupes = False):
    modifiedList = []
    for mut in mutationList:
        quiv = copy.deepcopy(mut[0].quiver)
        vertices = quiv.nodes
        vertexRelabeling = mut[2]
        isLineQuiver = True
        if len(quiv.edges) != len(vertices) - 1:
            isLineQuiver = False
        if bool(list(nx.simple_cycles(quiv))):
            isLineQuiver = False
        elif nx.dag_longest_path_length(quiv) != len(quiv.edges):
            isLineQuiver = False
        if isLineQuiver:
            if relabelNodes:
                (pathAlg, newVertexRelabeling) = relabelLineAlgebra(mut[0], vertexRelabeling)
            else:
                pathAlg = mut[0]
            rels = copy.deepcopy(pathAlg.rels)
            index = len(modifiedList)
            keepQuiver = True
            if discardLongerDupes:
                for modMut in modifiedList:
                    newRels = copy.deepcopy(rels)
                    oldRels = list(modMut[0].rels)
                    if newRels == oldRels and len(mut[1]) > len(modMut[1]):
                        keepQuiver = False
                        break
                    elif newRels < oldRels:
                        index = modifiedList.index(modMut)
                        break
            if keepQuiver:
                    modifiedList.insert(index, (pathAlg, mut[1], newVertexRelabeling))
    return modifiedList

def relationDualLineQuiver(quiverWrels):
    numberOfVertices = len(quiverWrels['quiver'])
    rels = quiverWrels['rels']
    dualRels = nx.MultiDiGraph()
    dualRels.add_nodes_from(rels.nodes)
    for rel in rels.edges:
        dualRels.add_edge(numberOfVertices - rel[1] + 1, numberOfVertices - rel[0] + 1)
    dualQuiverWrels = {'quiver' : quiverWrels['quiver'], 'rels' : dualRels}
    return dualQuiverWrels

def isRelationDualLineQuiver(quiverWrels1, quiverWrels2):
    numberOfVertices1 = len(quiverWrels1['quiver'])
    numberOfVertices2 = len(quiverWrels2['quiver'])
    rels1 = quiverWrels1['rels']
    rels2 = quiverWrels2['rels']
    isDualQuiver = False
    if numberOfVertices1 == numberOfVertices2:
        relSet1 = set(quiverWrels1['rels'].edges)
        relSet2 = set()
        for rel in quiverWrels2['rels'].edges:
            relSet2.add((numberOfVertices1 - rel[1] + 1, numberOfVertices1 - rel[0] + 1, 0))
        if relSet1 == relSet2:
            isDualQuiver = True
    return isDualQuiver

def saveLineRelationsToFile(fileName):
    mutationList = readMutationsFromFile('{0}.txt'.format(fileName))
    open('{0}Relations.txt'.format(fileName), 'w+').close()
    for mut in mutationList:
        relations = mut[0].rels
        with open('{0}Relations.txt'.format(fileName), "a") as f:
            f.write('{0}\n'.format(relations))
            f.close()
    return

def saveLineRelationsAndMutationsToFile(fileName, saveNumbering = False):
    mutationList = readMutationsFromFile('{0}.txt'.format(fileName))
    open('{0}RelationsAndMutations.txt'.format(fileName), 'w+').close()
    for mut in mutationList:
        mutations = mut[1]
        relations = mut[0].rels
        if saveNumbering:
            numbering = mut[2]
            with open('{0}RelationsAndMutations.txt'.format(fileName), "a") as f:
                f.write('Numbering: {0}\n'.format(numbering))
                f.write('Mutations: {0}\n'.format(mutations))
                f.write('Relations: {0}\n'.format(relations))
                f.write('\n')
                f.close()
        else:
            with open('{0}RelationsAndMutations.txt'.format(fileName), "a") as f:
                f.write('Mutations: {0}\n'.format(mutations))
                f.write('Relations: {0}\n'.format(relations))
                f.write('\n')
                f.close()
    return


def generateListOfRelations(listOfFileNames, combinedFileName = 'allRelations'):
    combinedMutationList = []
    for fileName in listOfFileNames:
        saveLineRelationsToFile(fileName)
        mutationList = readMutationsFromFile(fileName + '.txt')
        combinedMutationList.extend(mutationList[:])
    combinedMutationListClean = mutationListLineCleanup(combinedMutationList)
    open('{0}.txt'.format(combinedFileName), 'w+').close()
    for mut in combinedMutationListClean:
        saveLinePathAlgMutation(mut[0], mut[1], mut[2], combinedFileName + '.txt')
    saveLineRelationsToFile(combinedFileName)
    return

def readRelationsFromFile(fileName):
    relationSetList = []
    with open(fileName, 'r') as f:
        line = f.readline()
        while line != '':
            rRels = []
            if line[1] != ']':
                relArrows = list(line[1:-2].split(']], '))
                for relStr in relArrows:
                    if relStr[-1] != ']':
                        relStr = relStr + ']]'
                    relList = ast.literal_eval(relStr)
                    rRels.append(relList)
            relationSetList.append(rRels)
            line = f.readline()
    return relationSetList

def replaceSubPath(pathAlg, path, oldSubPath, newSubPath):
    newPath = path[:]
    if sublistExists(path, oldSubPath):
        newPath[path.index(oldSubPath[0]):path.index(oldSubPath[-1])] = newSubPath
        if not nx.is_path(pathAlg.quiver, newPath):
            print('New path is not a path in the quiver!\n')
            printPathAlgebra(pathAlg)
            print('Path: ', path)
            print('Old subpath: ', oldSubPath)
            print('New subpath: ', newSubPath)
            input('Press enter to continue...')
            return path
    else:
        print('Problem trying to replace subpath!\n')
        printPathAlgebra(pathAlg)
        print('Path: ', path)
        print('Old subpath: ', oldSubPath)
        print('New subpath: ', newSubPath)
        input('Press enter to continue...')
    return newPath

def powerset(iterable):
    "list(powerset([1,2,3])) --> [(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]"
    powerList = []
    s = list(iterable)
    for tup in itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)):
        powerList.append(list(tup))
    return powerList

def applyCommutativityRelSetToPath(path, relSet):
    newPath = path[:]
    for rel in relSet:
        if len(rel) == 2:
            if sublistExists(path, rel[0]):
                newPath[newPath.index(rel[0][0]):newPath.index(rel[0][-1])] = rel[1][:-1]
            elif sublistExists(path, rel[1]):
                newPath[newPath.index(rel[1][0]):newPath.index(rel[1][-1])] = rel[0][:-1]
    return newPath

def applyRelSetToPath(path, relSet):
    relSetCopy = copy.deepcopy(relSet)
    stillUnusedRels = True
    newPaths = [path[:]]
    while stillUnusedRels:
        relsToRemove = []
        usedRel = False
        stillUnusedRels = False
        for rel in relSetCopy:
            tempNewPaths = copy.deepcopy(newPaths)
            for newPath in tempNewPaths:
                for i in range(len(rel)):
                    if not (len(rel) == 1 and path == rel[0]):
                        if sublistExists(newPath, rel[i]):
                            for relPath in rel[:i] + rel[i + 1:]:
                                if not newPath[:newPath.index(rel[i][0])] + relPath[:-1] + newPath[newPath.index(rel[i][-1]):] in newPaths:
                                    newPaths.append(newPath[:newPath.index(rel[i][0])] + relPath[:-1] + newPath[newPath.index(rel[i][-1]):])
                            newPaths.remove(newPath)
                            usedRel = True
                            break
                while [] in newPaths:
                    newPaths.remove([])
            if usedRel:
                relsToRemove.append(rel)
        for rel in relsToRemove:
            relSetCopy.remove(rel)
        if bool(relSetCopy) and bool(relsToRemove):
            stillUnusedRels = True
    return newPaths

def pathHasZeroRel(path, relSet):
    hasZeroRel = False
    for rel in relSet:
        if len(rel) == 1:
            if sublistExists(path, rel[0]):
                hasZeroRel = True
                break
    return hasZeroRel

def numberOfPathsUpToRels(pathAlg, source, target):
    allPaths = nx.all_simple_paths(pathAlg.quiver, source, target)
    differentPaths = []
    numberOfDifferentPaths = 0
    rels = pathAlg.rels
    commutativityRels = []
    for rel in rels:
        if len(rel) == 2:
            commutativityRels.append(rel)
    allRelSets = []
    commutativeRelsPowerset = powerset(commutativityRels)
    for relSet in commutativeRelsPowerset:
        relSetPermutations = list(itertools.permutations(relSet))
        for relSetPermutation in relSetPermutations:
            relSetPermutationList = list(relSetPermutation)
            if not relSetPermutationList in allRelSets:
                allRelSets.append(relSetPermutationList)
    for path in allPaths:
        pathIsZero = pathHasZeroRel(path, rels)
        isSamePath = False
        for dPath in differentPaths:
            for relSet in allRelSets:
                cPath = applyRelSetToPath(path, relSet)[0]
                if cPath == dPath[0]:
                    isSamePath = True
                    if pathIsZero:
                        differentPaths[differentPaths.index(dPath)] = (dPath[0], pathIsZero)
                    break
            if isSamePath:
                break
        if not isSamePath:
            differentPaths.append((path, pathIsZero))
    for dPath in differentPaths:
        if not dPath[1]:
            numberOfDifferentPaths = numberOfDifferentPaths + 1
    return numberOfDifferentPaths


def cartanMatrix(pathAlg):
    quiv = pathAlg.quiver
    vertices = quiv.nodes
    cartanMatrix = eye(len(vertices), len(vertices))
    for i in vertices:
        for j in vertices:
            if i == j:
                cartanMatrix[j-1,i-1] = numberOfPathsUpToRels(pathAlg, i, j) + 1
            else:
                cartanMatrix[j-1,i-1] = numberOfPathsUpToRels(pathAlg, i, j)
    return cartanMatrix

def coxeterPoly(pathAlg):
    cartanMat = cartanMatrix(pathAlg)
    print(np.matrix(cartanMat))
    cartanMatInvTrans = cartanMat.inv().transpose()
    coxeterMatrix = -cartanMatInvTrans*cartanMat
    coxeterPolynomial = coxeterMatrix.charpoly()
    return coxeterPolynomial

def generateAllCoxeterPolynomials(length):

    # generate Kupisch series
    kupisch = [[[1]]]
    for i in range(1, length):
        kupisch.append([])
        for s in kupisch[i - 1]:
            for x in range(2, s[0] + 2):
                kupisch[i].append([x] + s)
    print('Number of different possible sets of relations: ',len(kupisch[length - 1]))

    polynomials = []
    for s in kupisch[length - 1]:
        M = sympy.Matrix([[1 if (m >= n and m < n + s[n]) else 0 for m in range(length)] for n in range(length)])
        C = - M * M.inv().transpose()
        p = sympy.factor(C.charpoly(sympy.Symbol("x")).as_expr())
        if p not in polynomials:
            polynomials.append(p)
    print('Coxeter polynomials: ',polynomials)
    print('Number of different Coxeter polynomials: ' ,len(polynomials))
    return polynomials

def generateAllPossibleLineRelations(lineLength):
    lineStart = 1
    lineStop = lineLength
    if lineLength <= 2:
        return [[]]
    elif lineLength == 2:
        return [[], [[[*range(lineStart, lineStop + 1)]]]]
    relSetList = generateAllPossibleLineRelations(lineLength - 1)
    allPossibleRelSets = relSetList[:]
    lastRelStart = 0
    for relSet in relSetList:
        if bool(relSet):
            lastRelStart = relSet[-1][0][0]
        for i in range(lastRelStart + 1, lineStop - 1):
            allPossibleRelSets.append(relSet + [[[*range(i,lineStop + 1)]]])
    allPossibleRelSets.sort()
    return allPossibleRelSets

def findMutationClassesForLine(lineLength, lineName, n_max = 5, manualDFdepth = 1, importAlreadyDoneSearches = False):
    start = time.time()
    DFgrowthFactor = 1 / math.ceil(lineLength / 2)
    mutationClassFiles = []
    allPossibleRelSets = generateAllPossibleLineRelations(lineLength)
    reachedRelSets = []
    allRelsReached = False
    firstUnreachedRelSet = []
    relSetNumber = '0'*(lineLength-2)
    edgeList = []
    for i in range(1, lineLength):
        edgeList.append((i, i+1))
    pathAlg = pathAlgebraClass.PathAlgebra()
    pathAlg.add_arrows_from(edgeList)
    relSetNumberAsList = ['0'] * (lineLength - 2)
    if importAlreadyDoneSearches:
        for name in glob.glob('{0}_*Relations.txt'.format(lineName)):
            mutationClassFiles.append(name)
            print('Filename: ', name)
        reachedAndMissingRels = collectMutationClasses(lineLength, saveToFile=True, printOutput=True)
        reachedRelSets = reachedAndMissingRels[0]
        for relSet in reachedRelSets:
            print('relSet: ', relSet)
        numberOfReachedRels = len(reachedRelSets)
        for n in range(len(allPossibleRelSets)):
            if n >= numberOfReachedRels or allPossibleRelSets[n] != reachedRelSets[n]:
                firstUnreachedRelSet = allPossibleRelSets[n]
                print('First unreached RelSet: ', firstUnreachedRelSet)
                allRelsReached = False
                break
            else:
                allRelsReached = True
        pathAlg.rels = copy.deepcopy(firstUnreachedRelSet)
        for rel in firstUnreachedRelSet:
            relSetNumberAsList[rel[0][0] - 1] = str(rel[0][-1] - rel[0][0])
            # rels.append([*range(rel[0][0], rel[0][-1])])
        relSetNumber = ''.join(relSetNumberAsList)
        print('RelSetNumber: ', relSetNumber)
    roundCount = 1
    while not allRelsReached:
        roundTimeStart = time.time()
        vertexRelabeling = {}
        maxRelLen = 0
        for numStr in relSetNumberAsList:
            maxRelLen = max(int(numStr), maxRelLen)
        DFdepth = max(maxRelLen - 1, math.floor(math.sqrt(lineLength)) + 2, manualDFdepth)
        quiverName = lineName + '_{0}'.format(relSetNumber)
        open('{0}DF.txt'.format(quiverName), 'w+').close()
        mutationSearchDepthFirst(pathAlg, DFdepth, [], quiverName, vertexRelabeling)
        print('\nFirst search for round done! \n')
        mutListDF = readMutationsFromFile('{0}DF.txt'.format(quiverName))
        mutList = mutationListLineCleanup(mutListDF)
        lapTimeStart = time.time()
        reachedRelations = []
        for mut in mutList:
            reachedRelations.append(mut[0].rels)
        checkedAlready = [False for i in range(len(mutList))]
        for n in range(1, n_max + 1):
            open('{0}DF.txt'.format(quiverName), 'w+').close()
            for mut in mutList:
                if checkedAlready[mutList.index(mut)]:
                    longestPathLength = nx.dag_longest_path_length(mut[0].quiver)
                    with open('{0}DF.txt'.format(quiverName), "a") as f:
                        f.write('Mutations: {0}\n'.format(mut[1]))
                        f.write('Numbering: {0}\n'.format(mut[2]))
                        f.write("Longest path: {0}\n".format(longestPathLength))
                        f.write('Vertices: {0}\n'.format(mut[0].quiver.nodes))
                        f.write('Arrows: {0}\n'.format(mut[0].quiver.edges))
                        f.write('Relations: {0}\n'.format(mut[0].rels))
                        f.write('-\n')
                        f.close()
                else:
                    print('Quiver name: ', quiverName)
                    print('Round {0}, lap {1}'.format(roundCount, n))
                    print('Quiver: ', quiverName)
                    print('Mutations: ', mut[1])
                    printPathAlgebra(mut[0])
                    DFdepthIncreaseByRound = max( -2, -math.floor(roundCount / math.sqrt(lineLength)))
                    DFdepthIncreaseByLap = 0
                    if roundCount > 4:
                        DFdepthIncreaseByLap = math.floor(n * DFgrowthFactor)
                    DFdepth = math.floor(math.sqrt(lineLength)) + 1 + DFdepthIncreaseByLap + DFdepthIncreaseByRound
                    mutationSearchDepthFirst(mut[0], DFdepth, mut[1], quiverName, mut[2])
            lapTimeEnd = time.time()
            print('Lap time for lap {0}: {1} s'.format(n, lapTimeEnd - lapTimeStart))
            lapTimeStart = time.time()
            mutListDF = readMutationsFromFile('{0}DF.txt'.format(quiverName))
            mutList = mutationListLineCleanup(mutListDF)
            checkedAlready = []
            open('{0}.txt'.format(quiverName), 'w+').close()
            for mut in mutList:
                saveLinePathAlgMutation(mut[0], mut[1], mut[2], '{0}.txt'.format(quiverName))
                checkedAlready.append(False)
                if mut[0].rels in reachedRelations:
                    checkedAlready[mutList.index(mut)] = True
                else:
                    reachedRelations.append(mut[0].rels)
        saveLineRelationsToFile('{0}'.format(quiverName))
        mutationClassFiles.append('{0}Relations.txt'.format(quiverName))
        reachedRelSetsThisRound = readRelationsFromFile('{0}Relations.txt'.format(quiverName))
        reachedRelSetsWithDupes = copy.deepcopy(reachedRelSets)
        reachedRelSetsWithDupes.extend((reachedRelSetsThisRound))
        reachedRelSetsWithDupes.sort()
        reachedRelSets = list(reachedRelSetsWithDupes for reachedRelSetsWithDupes,_ in itertools.groupby(reachedRelSetsWithDupes))
        numberOfReachedRels = len(reachedRelSets)
        for n in range(len(allPossibleRelSets)):
            print(allPossibleRelSets[n])
            if n >= numberOfReachedRels or allPossibleRelSets[n] != reachedRelSets[n]: #never satisfied if importAlreadyDoneSearches=True
                print(reachedRelSets[n])
                firstUnreachedRelSet = allPossibleRelSets[n]
                print('First unreached RelSet: ', firstUnreachedRelSet)
                allRelsReached = False
                break
            else:
                allRelsReached = True
        pathAlg.rels = copy.deepcopy(firstUnreachedRelSet)
        relSetNumberAsList = ['0']*(lineLength - 2)
        for rel in firstUnreachedRelSet:
            relSetNumberAsList[rel[0][0] - 1] = str(rel[0][-1] - rel[0][0])
            #rels.append([*range(rel[0][0], rel[0][-1])])
        relSetNumber = ''.join(relSetNumberAsList)
        print('RelSetNumber: ', relSetNumber)
        roundTimeStop = time.time()
        print('Round time for round {0}: {1}s'.format(roundCount, roundTimeStop - roundTimeStart))
        print('Looking at quiver {0}'.format(quiverName))
        roundCount = roundCount + 1
    open('{0}Relations.txt'.format(lineName), 'w+').close()
    for relSet in reachedRelSets:
        print('relSet being written to file: ', relSet)
        with open('{0}Relations.txt'.format(lineName), "a") as f:
            f.write('{0}\n'.format(relSet))
            f.close
    listOfRelSetLists = []
    mutationClassFiles.sort()
    for relFile in mutationClassFiles:
        print(relFile)
        listOfRelSetLists.append(readRelationsFromFile(relFile))
    mutationClassAdded = [False]*len(listOfRelSetLists)
    mutationClassesWithDupes = []
    for i in range(len(listOfRelSetLists)):
        currentMutationClassIndex = i
        if len(mutationClassesWithDupes) < i+1:
            if not mutationClassAdded[i]:
                mutationClassesWithDupes.append(listOfRelSetLists[i])
                currentMutationClassIndex = len(mutationClassesWithDupes) - 1
        if not mutationClassAdded[i]:
            relSetListGotAddedTo = True
            while relSetListGotAddedTo:
                relSetListGotAddedTo = False
                for j in range(i+1,len(listOfRelSetLists)):
                    if not mutationClassAdded[j]:
                        for relSet in listOfRelSetLists[j]:
                            if relSet in mutationClassesWithDupes[currentMutationClassIndex]:
                                mutationClassesWithDupes[currentMutationClassIndex].extend(listOfRelSetLists[j])
                                mutationClassAdded[j] = True
                                relSetListGotAddedTo = True
                                break
            mutationClassAdded[i] = True
    mutationClasses = []
    for relSetList in mutationClassesWithDupes:
        relSetList.sort()
        relSetListNoDupes = list(relSetList for relSetList, _ in itertools.groupby(relSetList))
        mutationClasses.append(relSetListNoDupes)
    if False:
        for relSetList1 in listOfRelSetLists[:]:
            print('relSetList1: ',relSetList1)
            for relSetList2 in listOfRelSetLists[listOfRelSetLists.index(relSetList1)+1:]:
                print('relSetList2: ', relSetList2)
                for relSet in relSetList2:
                    print('relSet: ', relSet)
                    print('relsSet in relSetn:', relSet in relSetList1)
                    if relSet in relSetList1:
                        relSetList1.extend(relSetList2)
                        print(relSetList1)
                        mutationClassAdded[listOfRelSetLists.index(relSetList2)] = True
                        break
            relSetList1.sort()
            relSetList = list(relSetList1 for relSetList1, _ in itertools.groupby(relSetList1))
            if not mutationClassAdded[listOfRelSetLists.index(relSetList1)]:
                mutationClasses.append(relSetList)
                mutationClassAdded[listOfRelSetLists.index(relSetList1)] = True
    for mutClass in mutationClasses:
        print('Mutation Class:')
        for relSet in mutClass:
            print(relSet)
        print()
    open('{0}MutationClasses.txt'.format(lineName), 'w+').close()
    for mutationClass in mutationClasses:
        with open('{0}MutationClasses.txt'.format(lineName), "a") as f:
            f.write('\n')
            f.write('--------------------------------------------------------------------------------------')
            f.write('\n')
            for relSet in mutationClass:
                f.write('{0}\n'.format(relSet))
            f.close
    with open('{0}MutationClasses.txt'.format(lineName), "a") as f:
        f.write('\n')
        f.write('--------------------------------------------------------------------------------------')
    print('Number of different sets of relations reached: ', len(reachedRelSets))
    print('Number of mutation classes: ', len(mutationClasses))
    end = time.time()
    print('Total runtime: {0}s'.format(end - start))
    return mutationClassFiles

def collectMutationClasses(lineLength, saveToFile = False, printOutput = False):
    lineName = 'A{0}'.format(lineLength)
    open('{0}MutationClasses.txt'.format(lineName), 'w+').close()
    mutationClassFiles = []
    for name in glob.glob('{0}_*Relations.txt'.format(lineName)):
        mutationClassFiles.append(name)
    mutationClassFiles.sort()
    listOfRelSetLists = []
    for relFile in mutationClassFiles:
        listOfRelSetLists.append(readRelationsFromFile('{0}'.format(relFile)))
    mutationClassAdded = [False]*len(listOfRelSetLists)
    mutationClassesWithDupes = []
    for i in range(len(listOfRelSetLists)):
        currentMutationClassIndex = i
        if len(mutationClassesWithDupes) < i+1:
            if not mutationClassAdded[i]:
                mutationClassesWithDupes.append(listOfRelSetLists[i])
                currentMutationClassIndex = len(mutationClassesWithDupes) - 1
        if not mutationClassAdded[i]:
            relSetListGotAddedTo = True
            while relSetListGotAddedTo:
                relSetListGotAddedTo = False
                for j in range(i+1,len(listOfRelSetLists)):
                    if not mutationClassAdded[j]:
                        for relSet in listOfRelSetLists[j]:
                            if relSet in mutationClassesWithDupes[currentMutationClassIndex]:
                                mutationClassesWithDupes[currentMutationClassIndex].extend(listOfRelSetLists[j])
                                mutationClassAdded[j] = True
                                relSetListGotAddedTo = True
                                break
            mutationClassAdded[i] = True
    mutationClasses = []
    for relSetList in mutationClassesWithDupes:
        relSetList.sort()
        relSetListNoDupes = list(relSetList for relSetList, _ in itertools.groupby(relSetList))
        mutationClasses.append(relSetListNoDupes)
    if False:
        for relSetList1 in listOfRelSetLists[:]:
            print('relSetList1:')
            for relSet in relSetList1:
                print(relSet)
            for relSetList2 in listOfRelSetLists[:]:
                if listOfRelSetLists.index(relSetList1) < listOfRelSetLists.index(relSetList2):
                    for relSet in relSetList2:
                        if relSet in relSetList1:
                            relSetList1.extend(relSetList2)
                            mutationClassAdded[listOfRelSetLists.index(relSetList2)] = True
                            break
                elif listOfRelSetLists.index(relSetList1) > listOfRelSetLists.index(relSetList2) and not mutationClassAdded[listOfRelSetLists.index(relSetList1)]:
                    for relSet in relSetList1:
                        if relSet in relSetList2:
                            mutationClasses[mutationClasses.index()]
                            mutationClassAdded[listOfRelSetLists.index(relSetList2)] = True
                            break
            relSetList1.sort()
            print('mutationClassAdded: ', mutationClassAdded)
            print('relSetList1: added =', mutationClassAdded[listOfRelSetLists.index(relSetList1)])
            for relSet in relSetList1:
                print(relSet)
            relSetList = list(relSetList1 for relSetList1, _ in itertools.groupby(relSetList1))
            print('relSetList:')
            for relSet in relSetList:
                print(relSet)
            if not mutationClassAdded[listOfRelSetLists.index(relSetList1)]:
                mutationClasses.append(relSetList)
                mutationClassAdded[listOfRelSetLists.index(relSetList1)] = True
    if printOutput:
        for mutClass in mutationClasses:
            for relSet in mutClass:
                print(relSet)
            print()
    numberOfReachedRelSets = 0
    allReachedRelSets = []
    if saveToFile:
        open('{0}MutationClasses.txt'.format(lineName), 'w+').close()
        for mutationClass in mutationClasses:
            numberOfReachedRelSets = numberOfReachedRelSets + len(mutationClass)
            allReachedRelSets.extend(mutationClass)
            with open('{0}MutationClasses.txt'.format(lineName), "a") as f:
                f.write('\n')
                f.write('--------------------------------------------------------------------------------------')
                f.write('\n')
                for relSet in mutationClass:
                    f.write('{0}\n'.format(relSet))
                f.close
    else:
        for mutationClass in mutationClasses:
            numberOfReachedRelSets = numberOfReachedRelSets + len(mutationClass)
            allReachedRelSets.extend(mutationClass)
    allReachedRelSets.sort()
    allPossibleRelSets = generateAllPossibleLineRelations(lineLength)
    missingRelSets = []
    for relSet in allPossibleRelSets:
        if not relSet in allReachedRelSets:
            missingRelSets.append(relSet)
    if printOutput:
        print('Missing rel sets: ')
        for relSet in missingRelSets:
           print(relSet)
        print('Number of different sets of relations reached: ', numberOfReachedRelSets)
        print('Number of mutation classes: ', len(mutationClasses))

    return (allReachedRelSets, missingRelSets)

def combineLineMutationFiles(lineLength):
    lineName = 'A{0}'.format(lineLength)
    open('{0}Lines.txt'.format(lineName), 'w+').close()
    mutationClassFiles = []
    for name in glob.glob('{0}_*[0-2].txt'.format(lineName)):
        print(name)
        mutationClassFiles.append(name)
        oneMutList = readMutationsFromFile(name)
        for mut in oneMutList:
            with open('{0}Lines.txt'.format(lineName), "a") as f:
                f.write('Mutations: {0}\n'.format(mut[1]))
                f.write('Numbering: {0}\n'.format(mut[2]))
                f.write('Vertices: {0}\n'.format(mut[0].vertices()))
                f.write('Arrows: {0}\n'.format(mut[0].arrows()))
                f.write('Relations: {0}\n'.format(mut[0].rels))
                f.write('-\n')
    reachedAndMissingRelations = collectMutationClasses(lineLength)
    mutListNotReduced = readMutationsFromFile('{0}Lines.txt'.format(lineName))
    mutList = []
    for mut in mutListNotReduced:
        keepMut = True
        for uniqueMut in mutList:
            if mut[0].rels == uniqueMut[0].rels:
                keepMut = False
                break
        if keepMut:
            mutList.append(mut)
    print('len(mutListNotReduced): ',len(mutListNotReduced))
    print('len(mutList): ', len(mutList))
    open('{0}Lines.txt'.format(lineName), 'w+').close()
    for mut in mutList:
        with open('{0}Lines.txt'.format(lineName), "a") as f:
            f.write('Mutations: {0}\n'.format(mut[1]))
            f.write('Numbering: {0}\n'.format(mut[2]))
            f.write('Vertices: {0}\n'.format(mut[0].vertices()))
            f.write('Arrows: {0}\n'.format(mut[0].arrows()))
            f.write('Relations: {0}\n'.format(mut[0].rels))
            f.write('-\n')
    return mutList

def makeStandardLineQuiver(lineLength, relationSet):
    pathAlg = pathAlgebraClass.PathAlgebra()
    pathAlg.add_vertices_from(range(1, lineLength + 1))
    pathAlg.add_arrows_from([[i, i+1] for i in range(1, lineLength)])
    pathAlg.add_rels_from(relationSet)
    return pathAlg

def readMutationClassesFromFile(lineLength, fileName):
    quiverList = []
    mutationClasses = []
    with open(fileName, 'r') as f:
        line = f.readline()
        while line != '':
            if line[0] == '-':
                if bool(quiverList):
                    mutationClasses.append(quiverList)
                    quiverList = []
            elif line[0] == '[':
                rRels = []
                if line[1] == ']':
                    relArrows = []
                else:
                    relArrows = list(line[1:-2].split(']], '))
                for relStr in relArrows:
                    if relStr[-1] != ']':
                        relStr = relStr + ']]'
                    relList = ast.literal_eval(relStr)
                    rRels.append(relList)
                quiverList.append(makeStandardLineQuiver(lineLength, rRels))
            line = f.readline()
    return mutationClasses

def combineMutationClasses(lineLength):
    mutClasses = readMutationClassesFromFile(lineLength, 'A{0}MutationClasses.txt'.format(lineLength))
    print(len(mutClasses))
    coxPolsForClasses = []
    for clas in mutClasses:
        coxPol = coxeterPoly(clas[0])
        coxPolsForClasses.append(coxPol)
        for quiv in clas[1:]:
            if coxeterPoly(quiv) != coxPol:
                print('DANGER! Coxeter polynomial does not match!')
                print('Class coxeter polynomial: ', coxPol)
                print('Quiver coxeter polynomial: ', coxeterPoly(quiv))
                printPathAlgebra(quiv)
                input('Press enter to continue...')
    potentiallySameClasses = []
    potentiallySameClassFound = []
    for i in range(len(coxPolsForClasses) - 1):
        potentiallySameClass = []
        if not i in potentiallySameClassFound:
            potentiallySameClass.append(i)
            for j in range(i + 1, len(coxPolsForClasses)):
                if not j in potentiallySameClassFound and coxPolsForClasses[j] == coxPolsForClasses[i]:
                    potentiallySameClass.append(j)
                    potentiallySameClassFound.append(j)
            potentiallySameClasses.append(potentiallySameClass)
    combinedMutClasses = []
    for classIndices in potentiallySameClasses:
        quiversInThisCombinedClass = []
        if len(classIndices) == 1:
            combinedMutClasses.append(mutClasses[classIndices[0]])
        else:
            quiversInThisCombinedClass.extend(mutClasses[classIndices[0]])
            for classIndex in classIndices[1:]:
                quiversInThisClass = []
                open('A{0}TempFileDF.txt'.format(lineLength), 'w+').close()
                for quiv in mutClasses[classIndex]:
                    mutationSearchDepthFirst(quiv, 6, [], 'A{0}TempFile'.format(lineLength), [])
                mutListDF = readMutationsFromFile('A{0}TempFileDF.txt'.format(lineLength))
                mutList = mutationListLineCleanup(mutListDF)
                isSameClass = False
                for mut in mutList:
                    quiversInThisClass.append(mut[0])
                    if mut[0] in quiversInThisCombinedClass:
                        isSameClass = True
                if isSameClass:
                    for quiv in quiversInThisClass:
                        if not quiv in quiversInThisCombinedClass:
                            quiversInThisCombinedClass.append(quiv)

def generateAllLineQuiversWithRelations(lineLength):
    allPossibleLineRelations = generateAllPossibleLineRelations(lineLength)
    allLineQuiversWithRelations = []
    for i in range(len(allPossibleLineRelations)):
        lineQuiver = makeStandardLineQuiver(lineLength, allPossibleLineRelations[i])
        allLineQuiversWithRelations.append(lineQuiver)
    return allLineQuiversWithRelations

def lineQuiverExample(lineLength, relationList, vertexRelabeling = {}):
    pathAlg = pathAlgebraClass.PathAlgebra()
    if len(relationList) != lineLength - 2:
        print('len(relationList) = ', len(relationList))
        print(lineLength - 2)
        print('Error! Invalid relation set.')
    elif (max(relationList) > lineLength - 1):
        print('Error! Invalid relation set.')
    else:
        vertices = list(range(1, lineLength + 1))
        if bool(vertexRelabeling):
            relabeledVertices = [vertexRelabeling[v] for v in vertices]
            vertices = relabeledVertices
        arrows = []
        for i in range(len(vertices)-1):
            arrows.append([vertices[i], vertices[i + 1]])
        rels = []
        for r in range(len(relationList)):
            if relationList[r] > 0:
                relStart = r + 1
                rel = [list(range(relStart, relStart + relationList[r] + 1))]
                if bool(vertexRelabeling):
                    relabeledRelPath = []
                    for i in range(len(rel[0])):
                        relabeledVertex = vertexRelabeling[rel[0][i]]
                        relabeledRelPath.append(relabeledVertex)
                    rel = [relabeledRelPath]
                rels.append(rel)
        pathAlg.add_vertices_from(vertices)
        pathAlg.add_arrows_from(arrows)
        pathAlg.add_rels_from(rels)
    return pathAlg

def createMutationClassCSV(lineLength):
    csvData = []
    allLineRels = generateAllPossibleLineRelations(lineLength)
    relSetsAsStrings = []
    for relSet in allLineRels:
        relSetsAsStrings.append(relSetToString(relSet))
    #relSetsAsStrings = nestedListToString(allLineRels)
    for relSetStr in relSetsAsStrings:
        print(relSetStr)
        csvData.append([relSetStr, '', '', '', ''])
    print(csvData)
    csvHeader = ['Relations', 'Mutation class', 'Mutation path from class representative', 'Coxeter polynomial', 'Numbering']
    with open('A_{0}_mutation_classes.csv'.format(lineLength), 'w') as file:
        writer = csv.writer(file)
        writer.writerow(csvHeader)
        writer.writerows(csvData)

    return

def importMutationClassCSV(filename):
    csvData = []
    with open(filename, newline='') as csvfile:
        csvReader = csv.reader(csvfile, delimiter=',')
        csvData = list(csvReader)
    #    for row in csvReader:
    #        csvData.append(', '.join(row))
    #for row in csvData:
    return csvData


def relSetToString(relSet):
    stringList = []
    for i in range(len(relSet)):
        stringInts = [str(int) for int in relSet[i][0]]
        stringOfInts = ";".join(stringInts)
        stringList.append(stringOfInts)
    joinedString = "|".join(stringList)
    return joinedString


def saveLineRelationsAndMutationsToCSV(fileName, mutationList, csvData, mutationClassName, printOutput = False):
    baseMutationVertexString = ''
    inputMutationClassName = mutationClassName
    minMutationLength = np.infty
    oldNumberingList = range(1, len(mutationClassName) + 3)
    baseVertexNumbering = {}
    if not bool(mutationList):
        return csvData
    for n in range(1, len(mutationList[0][2]) + 1):
        baseVertexNumbering[n] = n
    for mut in mutationList:
        vertexNumbering = mut[2]
        #vertexNumberingList = [vertexNumbering[n] for n in range(1, len(vertexNumbering) + 1)]
        mutationVertices = mut[1]
        coxPoly = coxeterPoly(mut[0])
        relSet = mut[0].rels
        relSetString = relSetToString(relSet)
        for i in range(len(csvData)):
            oldCSVrow = csvData[i]
            mutationLength = len(oldCSVrow[2]) + len(mut[1])
            if relSetString == oldCSVrow[0] and bool(oldCSVrow[1]) and mutationLength < minMutationLength:
                mutationClassName = oldCSVrow[1]
                oldNumberingString = oldCSVrow[4]
                if bool(oldNumberingString):
                    oldNumberingList = [int(s) for s in oldNumberingString.split(';')]
                else:
                    oldNumberingList = range(1, len(mut[0].vertices()) + 1)
                minusMutationVertexString = ''
                renumberedReverseMutationVertexString = ''
                if bool(mut[1]):
                    reverseMutationVertices = reverseMutationSequence(mutationVertices, vertexNumbering)
                    renumberedReverseMutationVerticesAsStringList = []
                    for n in range(len(reverseMutationVertices)):
                        if reverseMutationVertices[n] > 0:
                            renumberedVertex = oldNumberingList[reverseMutationVertices[n] - 1]
                        else:
                            renumberedVertex = -oldNumberingList[-reverseMutationVertices[n] - 1]
                        renumberedReverseMutationVerticesAsStringList.append(str(renumberedVertex))
                        #minusMutationVerticesAsStringList.append(str(-vertexNumbering[mutationVertices[-(n + 1)]]))
                    #minusMutationVerticesAsStringList = [str(-int) for int in reverseMutVertices]
                    renumberedReverseMutationVertexString = ';'.join(renumberedReverseMutationVerticesAsStringList)
                if bool(oldCSVrow[2]) and bool(renumberedReverseMutationVertexString):
                    baseMutationVertexString = ';'.join([oldCSVrow[2], renumberedReverseMutationVertexString])
                elif bool(renumberedReverseMutationVertexString):
                    baseMutationVertexString = renumberedReverseMutationVertexString
                else:
                    baseMutationVertexString = oldCSVrow[2]
                baseVertexNumbering = {}
                for n in range(1, len(vertexNumbering) + 1):
                    baseVertexNumbering[n] = oldNumberingList[getVertexNumberingKeyFromValue(vertexNumbering, n) - 1]
                minMutationLength = mutationLength
    for mut in mutationList:
        localVertexNumbering = mut[2]
        localVertexNumberingList = [localVertexNumbering[n] for n in range(1, len(localVertexNumbering) + 1)]
        vertexNumberingList = [baseVertexNumbering[localVertexNumbering[n]] for n in range(1, len(baseVertexNumbering) + 1)]
        vertexNumberingString = ';'.join([str(n) for n in vertexNumberingList])
        relSet = mut[0].rels
        relSetString = relSetToString(relSet)
        newRowAdded = False
        for i in range(len(csvData)):
            oldCSVrow = csvData[i]
            if relSetString == oldCSVrow[0] and not bool(oldCSVrow[1]):
                if bool(mut[1]):
                    localMutationVerticesAsStringList = [str(vertexNumberingList[localVertexNumberingList.index(int)]) for int in mut[1]]
                    localMutationVertexString = ';'.join(localMutationVerticesAsStringList)
                else:
                    localMutationVertexString = ''
                if bool(baseMutationVertexString) and bool(localMutationVertexString):
                    totalMutationVertexString = ';'.join([baseMutationVertexString, localMutationVertexString])
                elif bool(localMutationVertexString):
                    totalMutationVertexString = localMutationVertexString
                else:
                    totalMutationVertexString = baseMutationVertexString
                newCSVrow = [relSetString, mutationClassName, totalMutationVertexString, coxPoly.as_expr(), vertexNumberingString]
                csvData[i] = newCSVrow
                newRowAdded = True
            if newRowAdded:
                break
    print('as part of the class ', mutationClassName)
    if mutationClassName != inputMutationClassName:
        for row in csvData:
            if row[1] == inputMutationClassName:
                row[1] = mutationClassName
    if printOutput:
        for row in csvData:
            print(row)
    with open(fileName, 'w') as file:
        writer = csv.writer(file)
        writer.writerows(csvData)
    return csvData

def combineMutationClassesInCSVfile(filename, mutationDepth, lineLength):
    oldCSVdata = importMutationClassCSV(filename)
    newCSVdata = []
    for i in range(len(oldCSVdata)):
        baseRow = oldCSVdata[i]
        baseClass = baseRow[1]
        baseCoxPoly = baseRow[3]
        for j in range(i + 1, len(oldCSVdata)):
            currentRow = oldCSVdata[j]
            currentClass = currentRow[1]
            currentCoxPoly = currentRow[3]
            if currentCoxPoly == baseCoxPoly and currentClass != baseClass:
                newCSVdata.append(baseRow) #wrong!
                break

    mutationClassCSV = importMutationClassCSV(filename)
    numberOfCSVrows = len(mutationClassCSV)
    for n in range(mutationDepth):
        for i in range(numberOfCSVrows):
            baseRow = mutationClassCSV[i]
            baseClass = baseRow[1]
            baseCoxPoly = baseRow[3]
            for j in range(i + 1, numberOfCSVrows):
                currentRow = mutationClassCSV[j]
                currentClass = currentRow[1]
                currentMutationListOfStrings = currentRow[2].split(";")
                currentMutationList = [int(vert) for vert in currentMutationListOfStrings]
                currentCoxPoly = currentRow[3]
                if currentCoxPoly == baseCoxPoly and currentClass != baseClass and len(currentMutationListOfStrings) > mutationDepth - n:
                    relSetAsList = []
                    relSetAsListOfString = currentRow[0].split("|")
                    relSetAsListOfListOfString = [str.split(";") for str in relSetAsListOfString]
                    for stringList in relSetAsListOfListOfString:
                        if bool(stringList[0]):
                            listMap = map(int, stringList)
                            relList = list(listMap)
                            relSetAsList.append(relList)
                    lineRelList = [0]*(lineLength - 2)
                    for rel in relSetAsList:
                        lineRelList[rel[0] - 1] = len(rel) - 1
                    lineNumberStringList = [str(num) for num in lineRelList]
                    lineNumberString = "".join(lineNumberStringList)
                    quiverName = 'A{0}_{1}'.format(lineLength, lineNumberString)
                    pathAlg = lineQuiverExample(lineLength, lineRelList, currentRow[2])
                    printPathAlgebra(pathAlg)
                    open('{0}DF.txt'.format(quiverName), 'w').close()
                    mutationSearchDepthFirst(pathAlg, mutationDepth, [], quiverName)
                    mutList = readMutationsFromFile('{0}DF.txt'.format(quiverName))
                    cleanMutList = mutationListLineCleanup(mutList)
                    open('{0}.txt'.format(quiverName), 'w').close()
                    for mut in cleanMutList:
                        print('Mutations: {0}'.format(mut[1]))
                        print('Relations: {0}'.format(mut[0].rels))
                        saveLinePathAlgMutation(mut[0], mut[1], mut[2], '{0}.txt'.format(quiverName))
                    print('csv data: ', mutationClassCSV)
                    mutationClassCSV = saveLineRelationsAndMutationsToCSV('A_{0}_mutation_classes.csv'.format(lineLength), cleanMutList, mutationClassCSV, lineNumberString)
                    if mutationClassCSV[j][1] == baseClass:
                        baseQuiverRelations = [int(rel) in baseClass.split('')]
                        baseQuiver = lineQuiverExample(len(baseQuiverRelations) + 2, baseQuiverRelations)
                        quiverWithRightNumbering = quiverMutationAtVertices(baseQuiver, )
                        mutDictKeys = list()
                        reverseMutationVertices = []
                        for k in range(numberOfCSVrows):
                            additionalRow = mutationClassCSV[k]
                            additionalClass = additionalRow[1]
                            if additionalClass == currentClass:
                                mutationClassCSV[k][1] = baseClass
                                mutationClassCSV[k][2] = mutationClassCSV[j][2] + ';'
    return newCSVdata

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

def reverseMutationSequence(mutationVertices, vertexNumbering):
    reverseMutationVertices = []
    for n in range(len(mutationVertices) - 1, -1, -1):
        reverseMutationVertices.append(-getVertexNumberingKeyFromValue(vertexNumbering, mutationVertices[n]))
    return reverseMutationVertices

def reverseMutationFromSequence(pathAlg, mutationVertices, vertexNumbering):
    reverseMutationVertices = reverseMutationSequence(mutationVertices, vertexNumbering)
    return quiverMutationAtVertices(pathAlg, reverseMutationVertices)

def coxPolyOfTree(tree):
    adjMat = sympy.Matrix(nx.adjacency_matrix(tree).todense(), dtype=int)
    triuMat = np.triu(Matrix(adjMat))
    mat = sympy.Matrix(triuMat, dtype=int) + sympy.eye(len(tree))#sympy.Matrix(np.triu(Matrix(matrixAsList))) + sympy.eye(n)
    print(np.matrix(mat))
    matInvTrans = mat.inv().transpose()
    coxeterMatrix = -matInvTrans * mat
    coxeterPolynomial = coxeterMatrix.charpoly()
    #print(coxeterPolynomial.as_expr())
    return coxeterPolynomial

def generateAllKupischSeries(length):

    # generate Kupisch series
    kupisch = [[[1]]]
    for i in range(1, length):
        kupisch.append([])
        for s in kupisch[i - 1]:
            for x in range(2, s[0] + 2):
                kupisch[i].append([x] + s)
    print('Number of different possible kupisch series: ',len(kupisch[length - 1]))
    return kupisch

def expandClassWith2Rels(baseLineRelList):
    lineLength = len(baseLineRelList) + 2
    expanding2RelSetLists = []
    possible2RelPositions = []
    for i in range(len(baseLineRelList)):
        add2RelPos = False
        if baseLineRelList[i] == 2:
            add2RelPos = True
        elif baseLineRelList[i] == 0:
            add2RelPos = True
            if i >= 1:
                for j in range(i):
                    if baseLineRelList[i-1-j] > j + 2:
                        add2RelPos = False
                        break
        if add2RelPos:
            possible2RelPositions.append(i)
    relevantPowerSet = powerset(possible2RelPositions)
    for relPosSet in relevantPowerSet:
        lineRelList = baseLineRelList.copy()
        compRelPosSet = [x for x in possible2RelPositions if x not in relPosSet]
        for p in relPosSet:
             lineRelList[p] = 2
        for p in compRelPosSet:
            lineRelList[p] = 0
        if lineRelList != baseLineRelList:
            print(lineRelList)
            expanding2RelSetLists.append(lineRelList)
    expanding2RelSetMutList = []
    for lineRelList in expanding2RelSetLists:
        pathAlg = lineQuiverExample(lineLength, lineRelList)
        expanding2RelSetMutList.append((pathAlg, [], {}))
    return expanding2RelSetMutList

def expandClassFurtherWithEqualRelPairs(baseMutList):
    lineLength = len(baseMutList[0][0].vertices())
    lineRelListsToAdd = []
    expandingRelPairMutList = copy.deepcopy(baseMutList)
    for mut in baseMutList:
        pathAlg = mut[0]
        rels = pathAlg.rels
        lineRelList = [0] * (lineLength - 2)
        for rel in rels:
            lineRelList[rel[0][0] - 1] = len(rel[0]) - 1
        lineRelListsToAddFromThisMut = []
        noCandidatePairs = True
        for i, j in enumerate(lineRelList[:-1]):
            if j >= 3 and j==lineRelList[i+1]:
                noCandidatePairs = False
        if noCandidatePairs:
            break
        if not lineRelList in lineRelListsToAdd:
            addMutations = False
            for i, j in enumerate(lineRelList[:-1]):
                addMutations = False
                numberOfForwardMutations = 0
                numberOfBackwardMutations = 0
                if j >= 3 and j == lineRelList[i+1]:
                    addMutations = True
                    for k in range(0, j-2):
                        if lineRelList[i+2+k] > 0:
                            addMutations = False
                    if addMutations:
                        for k in range(i+j+1, lineLength - 2):
                            if lineRelList[k] == 0:
                                numberOfForwardMutations += 1
                            else:
                                break
                        numberOfBackwardMutations = i + 1
                        for k in range(0,i):
                            if lineRelList[k] >= 2:
                                numberOfBackwardMutations = i + 1 - (k + lineRelList[i])
                                if numberOfBackwardMutations < 0:
                                    addMutations = False
                                    break
                if addMutations:
                    lineRelListToAdd = lineRelList.copy()
                    for k in range(numberOfForwardMutations):
                        lineRelListToAdd[i+k] = 0
                        lineRelListToAdd[i+k+2] = j
                        print(lineRelListToAdd)
                        lineRelListsToAddFromThisMut.append(lineRelListToAdd.copy())
                    lineRelListToAdd = lineRelList.copy()
                    for k in range(numberOfBackwardMutations):
                        lineRelListToAdd[i+2-k] = 0
                        lineRelListToAdd[i-k] = j
                        print(lineRelListToAdd)
                        lineRelListsToAddFromThisMut.append(lineRelListToAdd.copy())
        lineRelListsToAdd.extend(lineRelListsToAddFromThisMut)
    for lineRelListToAdd in lineRelListsToAdd:
        pathAlg = lineQuiverExample(lineLength, lineRelListToAdd)
        if not (pathAlg, [], {}) in expandingRelPairMutList:
            expandingRelPairMutList.append((pathAlg, [], {}))
    return expandingRelPairMutList

def expandAllClassesWithEasyRels(lineLength, startRow = 0):
    mutationClassCSV = importMutationClassCSV('A_{0}_mutation_classes.csv'.format(lineLength))
    numberOfCSVrows = len(mutationClassCSV)
    for i in range(startRow + 1, numberOfCSVrows):
        row = mutationClassCSV[i]
        print('csv row: ', row)
        if not (bool(row[1]) or bool(row[2])):
            relSetAsList = []
            relSetAsListOfString = row[0].split("|")
            relSetAsListOfListOfString = [str.split(";") for str in relSetAsListOfString]
            for stringList in relSetAsListOfListOfString:
                if bool(stringList[0]):
                    listMap = map(int, stringList)
                    relList = list(listMap)
                    print(relList)
                    relSetAsList.append(relList)
            print(relSetAsList)
            lineRelList = [0] * (lineLength - 2)
            for rel in relSetAsList:
                lineRelList[rel[0] - 1] = len(rel) - 1
            lineNumberStringList = [str(num) for num in lineRelList]
            lineNumberString = "".join(lineNumberStringList)
            mutListWith2Rels = expandClassWith2Rels(lineRelList)
            if bool(mutListWith2Rels):
                mutListWith2RelsAndPairs = expandClassFurtherWithEqualRelPairs(mutListWith2Rels)
                mutationClassCSV = saveLineRelationsAndMutationsToCSV('A_{0}_mutation_classes.csv'.format(lineLength),
                                                                  mutListWith2RelsAndPairs, mutationClassCSV, lineNumberString)
    return

def mutationSearch(lineLength, mutationDepthStart, startRow, createNewCSVfile = False, printMutations = False):
    mutationDepth = mutationDepthStart
    if createNewCSVfile:
        createMutationClassCSV(lineLength)
    mutationClassCSV = importMutationClassCSV('A_{0}_mutation_classes.csv'.format(lineLength))
    numberOfCSVrows = len(mutationClassCSV)
    for i in range(numberOfCSVrows):
        row = mutationClassCSV[i + startRow - (numberOfCSVrows - 1)]
        print('csv row: ', row)
        if not bool(row[1]):
            relSetAsList = []
            relSetAsListOfString = row[0].split("|")
            relSetAsListOfListOfString = [str.split(";") for str in relSetAsListOfString]
            for stringList in relSetAsListOfListOfString:
                if bool(stringList[0]):
                    listMap = map(int, stringList)
                    relList = list(listMap)
                    relSetAsList.append(relList)
            print('Relations in class: ', relSetAsList)
            lineRelList = [0] * (lineLength - 2)
            for rel in relSetAsList:
                lineRelList[rel[0] - 1] = len(rel) - 1
            lineNumberStringList = [str(num) for num in lineRelList]
            lineNumberString = "".join(lineNumberStringList)
            print('Working on class {0}'.format(lineNumberString))
            quiverName = 'A{0}_{1}'.format(lineLength, lineNumberString)
            pathAlg = lineQuiverExample(lineLength, lineRelList, row[2])
            if printMutations:
                printPathAlgebra(pathAlg)
            open('{0}DF.txt'.format(quiverName), 'w').close()
            print('Mutation depth: {0}'.format(mutationDepth))
            mutationSearchDepthFirst(pathAlg, mutationDepth, [], quiverName, printOutput=printMutations)
            mutList = readMutationsFromFile('{0}DF.txt'.format(quiverName))
            cleanMutList = mutationListLineCleanup(mutList, printOutput=printMutations)
            open('{0}.txt'.format(quiverName), 'w').close()
            for mut in cleanMutList:
                print('Mutations: {0}'.format(mut[1]))
                print('Relations: {0}'.format(mut[0].rels))
                saveLinePathAlgMutation(mut[0], mut[1], mut[2], '{0}.txt'.format(quiverName))
            print('Writing class {0} to csv file'.format(lineNumberString))
            mutationClassCSV = saveLineRelationsAndMutationsToCSV('A_{0}_mutation_classes.csv'.format(lineLength),
                                                                  cleanMutList, mutationClassCSV, lineNumberString)
        mutationDepth = np.maximum(mutationDepthStart - np.floor(np.log10(i+1)), 2)
    return

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

def convertLineFromCSVnotation( lineLength, lineInCSVnotation ):
    #lineInCSVnotation'1;2;3;4;5|2;3;4;5;6|4;5;6;7;8;9;10'
    listOfRels = []
    listOfRelsAsStr = lineInCSVnotation.split('|')
    for relString in listOfRelsAsStr:
        relAsList = [ int(i) for i in relString.split(';') ]
        listOfRels.append([relAsList])
    lineQuiver = makeStandardLineQuiver(lineLength, listOfRels)
    return lineQuiver

def generateAllQuipusUpToLength( length ):
    timeStart = time.time()
    quipusOfAllLengths = [[nx.path_graph(1)]]
    pathGraphs = []
    for i in range(length):
        pathGraphs.append(nx.path_graph(i+1))
        quipusOfThisLength = []
        for j in range(int(np.floor((i+3)/2)), i+1 ):
            pathGraph = pathGraphs[j]
            for vertexSet in itertools.combinations(range(1,j), i + 1 - j):
                heightOneQuipu = pathGraph.copy()
                k = 1
                for v in vertexSet:
                    heightOneQuipu.add_edge(v, j+k)
                    k += 1
                isNewQuipu = True
                for Q in quipusOfThisLength:
                    if nx.is_isomorphic(heightOneQuipu, Q):
                        isNewQuipu = False
                        break
                if isNewQuipu:
                    quipusOfThisLength.append(heightOneQuipu)
        for quipu in quipusOfAllLengths[i]:
            for v in quipu.nodes():
                if quipu.degree(v) <= 1:
                    longerQuipu = quipu.copy()
                    longerQuipu.add_edge(v, i+1)
                    isNewQuipu = True
                    for Q in quipusOfThisLength:
                        if nx.is_isomorphic(longerQuipu, Q):
                            isNewQuipu = False
                            break
                    if isNewQuipu:
                        quipusOfThisLength.append(longerQuipu)
        quipusOfAllLengths.append([])
        for newQuipu in quipusOfThisLength:
            quipusOfAllLengths[i+1].append(newQuipu.copy())
        timeEnd = time.time()
        print('Generated ', len(quipusOfThisLength), ' quipus of length ', i+2, ' in ',  timeEnd - timeStart, 's.')
    return quipusOfAllLengths

def generateAllHeightOneQuipus( mainStringLength, mainStringLengthStart = 1 ):
    startTime = time.time()
    heightOneQuipusWithDupes = [nx.path_graph(1)]
    for i in range(mainStringLengthStart,mainStringLength):
        heightOneQuipusWithDupes.append(nx.path_graph(i+1))
        vertexPowerset = powerset(range(1,i))
        j = 0
        while j < 2**(i-1):
            vertexSet = vertexPowerset[j]
            setToRemove = {}
            for v in vertexSet:
                setToRemove.append(i-v)
            if setToRemove:
                print('before', vertexPowerset)
                vertexPowerset = vertexPowerset - setToRemove
                print('after', vertexPowerset)
                j += 1
            j += 1
        for vertexSet in vertexPowerset:
            heightOneQuipu = nx.path_graph(i+1)
            if bool(vertexSet):
                k = 1
                for v in vertexSet:
                    heightOneQuipu.add_edge(v, i+k)
                    k += 1
                heightOneQuipusWithDupes.append(heightOneQuipu)
        print('Done with ', i+1, ' in ', time.time() - startTime, 's')
    heightOneQuipus = []
    print('Found ', len(heightOneQuipusWithDupes), ' height one quipus with dupes.')
    for Q1 in heightOneQuipusWithDupes:
        isNewQuipu = True
        for Q2 in heightOneQuipus:
            if nx.is_isomorphic(Q1, Q2):
                isNewQuipu = False
                break
        if isNewQuipu:
            heightOneQuipus.append(Q1)
    return heightOneQuipus

def saveQuipusToCSV( quipusList, fileName, overwrirteFile = True ):
    if overwrirteFile:
        open(fileName, 'w+').close()
    with open(fileName, 'a') as f:
        for Q in quipusList:
            f.write('{0}\n'.format(str(Q.edges)))
        f.close()
    return

def cartanMatrixForCanonicalAlgebra(pathAlg):
    quiv = pathAlg.quiver
    vertices = quiv.nodes
    cartanMatrix = eye(len(vertices), len(vertices))
    for i in vertices:
        for j in vertices:
            if i == j:
                cartanMatrix[j-1,i-1] = 1
            else:
                allPaths = list(nx.all_simple_paths(pathAlg.quiver, i, j))
                cartanMatrix[j-1,i-1] = min(len(allPaths), 2)
    return cartanMatrix

def coxeterPolyForCanonicalAlgebra(pathAlg):
    cartanMat = cartanMatrixForCanonicalAlgebra(pathAlg)
    print(np.matrix(cartanMat))
    cartanMatInvTrans = cartanMat.inv().transpose()
    coxeterMatrix = -cartanMatInvTrans*cartanMat
    coxeterPolynomial = coxeterMatrix.charpoly()
    return coxeterPolynomial

