"""Cleaning up after a mutation.

`quiverMutationAtVertex` applies steps 1-7 of the procedure of
arXiv:2112.08129 and stops there, so its result routinely contains relations the
paper's "Note" after step 7 says to cancel, relations that are not minimal, and
duplicates.  `reducePathAlgebra` is that cleanup, run to a fixed point; the rest
of this module is its individual passes, in the order it applies them.
"""

import copy
import itertools


from . import pathAlgebra
from . import paths



def reducePathAlgebra(pathAlg):
    quiver = pathAlg.quiver
    redPathAlg = pathAlgebra.PathAlgebra()
    redPathAlg.add_vertices_from(quiver.nodes)
    redPathAlg.add_arrows_from(quiver.edges(keys=True))
    for rel in pathAlg.rels:
        if paths.isIllegalRelation(pathAlg, rel):
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
                        if paths.sublistExists(relPath, removedRel[0]):
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
        commutativityRelsPowerSet = paths.powerset(commutativityRels)
        for relSet in commutativityRelsPowerSet:
            for perm in itertools.permutations(relSet):
                relSetsToApply.append(list(perm))
    for rel1 in zeroRels:
        removeRel = False
        for relSetToApply in relSetsToApply:
            if applyCommutativityRels:
                rel1Equivalent = paths.applyRelSetToPath(rel1[0], relSetToApply)
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


def removeRedundantRelations(pathAlg):
    #sort relations by increasing length of their longest path, then by increasing number of paths
    shortestRelsList = sorted(pathAlg.rels, key=lambda x: (len(max(x, key=len)), len(x)))
    necesarryRels = []
    while bool(shortestRelsList):
        necesarryRels.append(shortestRelsList.pop(0))
        relSetsToApply = paths.powerset(necesarryRels)
        relsToRemove = []
        for rel in shortestRelsList:
            removeRel = False
            for relPath in rel:
                for relsToApply in relSetsToApply[1:]:
                    newPaths = paths.applyRelSetToPath(relPath, relsToApply)
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
