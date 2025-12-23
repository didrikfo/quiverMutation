import copy
import itertools
import networkx as nx
import path_algebra_class
from quiver_mutation_io import print_path_algebra
from quiver_mutation_utils import powerset, sublist_exists

def non_minimal_out_rels(pathAlg, vertex):
    nonMinOutRels = []
    pathAlgCopy = copy.deepcopy(pathAlg)
    for rel in pathAlgCopy.out_rels(vertex):
        nonMinOutRels.extend(extend_rel(pathAlgCopy, rel))
    for ar in pathAlgCopy.out_arrows(vertex):
        deeperOutRels = non_minimal_out_rels(pathAlgCopy, ar[1])
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

def all_rels_between_vertices(pathAlg, startVertex,endVertex):
    pathAlgCopy = copy.deepcopy(pathAlg)
    relsBetween = pathAlgCopy.rels_between(startVertex, endVertex)
    #allRelsBetween= pathAlgCopy.rels_between(startVertex, endVertex)
    for ar in pathAlgCopy.out_arrows(startVertex):
        dRelsBetween = all_rels_between_vertices(pathAlg, ar[1], endVertex)
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
                    newRelPaths = apply_rel_set_to_path(relPath, relSet)
                    relToAdd = sorted(rel[:rel.index(relPath)] + rel[rel.index(relPath) + 1:] + newRelPaths)
                    if not relToAdd in allRelsBetween and not any(relToAdd.count(x) > 1 for x in relToAdd) and relToAdd != []:
                        allRelsBetween.append(copy.deepcopy(relToAdd))
                        newRels.append(relToAdd)
                        stillNewRels = True
        relSetsToApply = [[copy.deepcopy(newRel)] for newRel in newRels]
    return allRelsBetween

def all_minimal_rels_between_vertices(pathAlg, startVertex,endVertex):
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
                newRelPaths = apply_rel_set_to_path(relPath, relSet)
                relToAdd = sorted(rel[:rel.index(relPath)] + rel[rel.index(relPath) + 1:] + newRelPaths)
                if not relToAdd in allRelsBetween and not any(relToAdd.count(x) > 1 for x in relToAdd) and relToAdd != []:
                    allRelsBetween.append(relToAdd)
    return allRelsBetween

def extend_rel(pathAlg, rel):
    vertex = rel[0][-1]
    extendedRels = [rel]
    outArrows = pathAlg.out_arrows(vertex)
    for ar in outArrows:
        extendedRel = []
        for relPath in rel:
            extendedRelPath = relPath + [ar[1]]
            extendedRel.append(extendedRelPath)
        extendedRels.append(extendedRel)
        deeperExtendedRels = extend_rel(pathAlg, extendedRel)
        extendedRels.extend(deeperExtendedRels)
    return extendedRels

def is_sub_rel_of(potentialSubRel, relation):
    isSubRel = True
    for relPath in potentialSubRel:
        if not relPath in relation:
            isSubRel = False
            break
    return isSubRel

def reduce_path_algebra(pathAlg):
    quiver = pathAlg.quiver
    redPathAlg = path_algebra_class.PathAlgebra()
    redPathAlg.add_vertices_from(quiver.nodes)
    redPathAlg.add_arrows_from(quiver.edges(keys=True))
    for rel in pathAlg.rels:
        if is_illegal_relation(pathAlg, rel):
            pathAlg.rels.remove(rel)
    redPathAlg.add_rels_from(pathAlg.rels)
    redPathAlg = remove_duplicate_rel_paths(redPathAlg)
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
                        if sublist_exists(relPath, removedRel[0]):
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
        redPathAlg = remove_duplicate_rel_paths(redPathAlg)
    redPathAlg = remove_duplicate_rels(redPathAlg)
    redPathAlg = remove_nonminimal_zero_rels(redPathAlg, applyCommutativityRels=False)
    remove_redundant_relations(redPathAlg)
    remove_existing_subrelations(redPathAlg)
    noChange = False
    while not noChange:
        numberOfRelsBeforeRed = len(redPathAlg.rels)
        remove_redundant_relations(redPathAlg)
        if numberOfRelsBeforeRed == len(redPathAlg.rels):
            noChange = True
    redPathAlg = remove_nonminimal_zero_rels(redPathAlg)
    return redPathAlg

def remove_nonminimal_zero_rels(pathAlg, applyCommutativityRels = True):
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
                rel1Equivalent = apply_rel_set_to_path(rel1[0], relSetToApply)
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

def reduce_commutativity_rels(pathAlg):
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

def remove_duplicate_rels(pathAlg):
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

def remove_duplicate_rel_paths(pathAlg):
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

def minimize_commuting_relation(pathAlg, relation):
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
            for rel in all_rels_between_vertices(pathAlg, i, j):
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

def remove_redundant_relations(pathAlg):
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
                    newPaths = apply_rel_set_to_path(relPath, relsToApply)
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

def remove_existing_subrelations(pathAlg):
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

def zeroize_rels(rels):
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
                if sublist_exists(relPath, zeroRel[0]):
                    rel.remove(relPath)
                    isZeroPath = True
                    break
        zeroizedRels.append(rel)
    zeroizedRelsReduced = [rel for rel in zeroizedRels if rel != []]
    return zeroizedRelsReduced

def is_illegal_relation(pathAlgebra, relation):
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
        print_path_algebra(pathAlgebra)
        # input('Press enter to continue...')
    return isIllegal

def replace_sub_path(pathAlg, path, oldSubPath, newSubPath):
    newPath = path[:]
    if sublist_exists(path, oldSubPath):
        newPath[path.index(oldSubPath[0]):path.index(oldSubPath[-1])] = newSubPath
        if not nx.is_path(pathAlg.quiver, newPath):
            print('New path is not a path in the quiver!\n')
            print_path_algebra(pathAlg)
            print('Path: ', path)
            print('Old subpath: ', oldSubPath)
            print('New subpath: ', newSubPath)
            input('Press enter to continue...')
            return path
    else:
        print('Problem trying to replace subpath!\n')
        print_path_algebra(pathAlg)
        print('Path: ', path)
        print('Old subpath: ', oldSubPath)
        print('New subpath: ', newSubPath)
        input('Press enter to continue...')
    return newPath

def apply_commutativity_rel_set_to_path(path, relSet):
    newPath = path[:]
    for rel in relSet:
        if len(rel) == 2:
            if sublist_exists(path, rel[0]):
                newPath[newPath.index(rel[0][0]):newPath.index(rel[0][-1])] = rel[1][:-1]
            elif sublist_exists(path, rel[1]):
                newPath[newPath.index(rel[1][0]):newPath.index(rel[1][-1])] = rel[0][:-1]
    return newPath

def apply_rel_set_to_path(path, relSet):
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
                        if sublist_exists(newPath, rel[i]):
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

def path_has_zero_rel(path, relSet):
    hasZeroRel = False
    for rel in relSet:
        if len(rel) == 1:
            if sublist_exists(path, rel[0]):
                hasZeroRel = True
                break
    return hasZeroRel

def number_of_paths_up_to_rels(pathAlg, source, target):
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
        pathIsZero = path_has_zero_rel(path, rels)
        isSamePath = False
        for dPath in differentPaths:
            for relSet in allRelSets:
                cPath = apply_rel_set_to_path(path, relSet)[0]
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
