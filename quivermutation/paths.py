"""Paths and relations inside a path algebra.

The layer everything else is built on: enumerating the paths between two
vertices, the relations between them, which relations are minimal, and what a
relation set does to a path.  A path is a list of vertices and a relation is a
list of paths (see NOTES.md, "The model of a path algebra"), so all of this is
list manipulation with the quiver consulted for the arrows.

Nothing here mutates or reduces; `mutation` and `reduction` do that.
"""

import copy
import itertools

import networkx as nx

from . import pathAlgebra



def listIntersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


def sublistExists(list, sublist):
    for i in range(len(list)-len(sublist)+1):
        if sublist == list[i:i+len(sublist)]:
            return True #return position (i) if you wish
    return False


def powerset(iterable):
    "list(powerset([1,2,3])) --> [(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]"
    powerList = []
    s = list(iterable)
    for tup in itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)):
        powerList.append(list(tup))
    return powerList


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
    # networkx >= 3.1 yields the trivial length-zero path when source == target.
    # This function counts non-trivial paths only (the trivial path is accounted
    # for separately by the +1 on the diagonal of the Cartan matrix), so drop it.
    allPaths = (p for p in nx.all_simple_paths(pathAlg.quiver, source, target) if len(p) > 1)
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


def allRelsBetweenVertices(pathAlg, startVertex, endVertex, visited = None):
    """Every relation from startVertex to endVertex, minimal or not.

    `visited` carries the vertices already on the current recursion path.  The
    recursion follows the arrows out of startVertex, and without that set a
    cycle in the quiver makes it descend forever: allRelsBetweenVertices on a
    quiver with a cycle raised RecursionError.  Since the rest of the module
    works with simple paths throughout, not revisiting a vertex is also the
    right semantics, and on an acyclic quiver it changes nothing -- no vertex
    can repeat on a path there anyway.

    This is the crash that stopped the length-12 run.  mutationSearchDepthFirst
    calls allRelsInPathAlgebra at the top of every node, and it used to do so
    before testing the quiver for cycles, so the first mutation that produced a
    cyclic quiver killed the search on the following node.  Both halves are
    fixed: the recursion is bounded here, and the search now tests for cycles
    first.
    """
    # This used to deep-copy the whole path algebra, including its networkx
    # graph, on every one of its recursive calls, which accounted for most of
    # the run time of a class search.  The quiver is only read here, so the
    # relations are all that need copying, and rels_between already returns a
    # fresh list.
    visited = frozenset() if visited is None else visited
    visited = visited | {startVertex}
    relsBetween = [copy.deepcopy(rel) for rel in pathAlg.rels_between(startVertex, endVertex)]
    for ar in pathAlg.out_arrows(startVertex):
        if ar[1] in visited:
            continue
        dRelsBetween = allRelsBetweenVertices(pathAlg, ar[1], endVertex, visited)
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
                for rel in pathAlg.rels_between(startVertex, i):
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
    # As in allRelsBetweenVertices, the quiver is only read, so copying the
    # relations is enough and copying the graph with it was pure overhead.
    relsBetween = [copy.deepcopy(rel) for rel in pathAlg.rels_between(startVertex, endVertex)]
    allRelsBetween = [copy.deepcopy(rel) for rel in pathAlg.rels_between(startVertex, endVertex)]
    verticesBetween = []
    for path in nx.all_simple_paths(pathAlg.quiver, startVertex, endVertex):
        for i in path:
            if not i in verticesBetween:
                verticesBetween.append(i)
    intermideateRels = []
    for i in verticesBetween:
        for j in verticesBetween:
            if i != startVertex or j != endVertex:
                for rel in pathAlg.rels_between(i, j):
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


def extendRel(pathAlg, rel, visited = None):
    """Every way of extending a relation forward along the arrows out of its end.

    As in allRelsBetweenVertices, `visited` bounds the recursion so a cycle in
    the quiver does not make it descend forever.  It is seeded with every vertex
    the relation already passes through, so an extension never doubles back into
    the relation itself.

    Note that each extension is returned twice, since the recursion's own result
    starts with the relation it was given.  That was harmless when
    `nonMinimalOutRels` was the only caller and deduped at the end; that caller
    is gone, so a caller of this now has to dedupe for itself.
    """
    vertex = rel[0][-1]
    if visited is None:
        # Seed with every vertex the relation already passes through, so an
        # extension cannot double back into it.
        visited = frozenset(v for relPath in rel for v in relPath)
    visited = visited | {vertex}
    extendedRels = [rel]
    outArrows = pathAlg.out_arrows(vertex)
    for ar in outArrows:
        if ar[1] in visited:
            continue
        extendedRel = []
        for relPath in rel:
            extendedRelPath = relPath + [ar[1]]
            extendedRel.append(extendedRelPath)
        extendedRels.append(extendedRel)
        deeperExtendedRels = extendRel(pathAlg, extendedRel, visited)
        extendedRels.extend(deeperExtendedRels)
    return extendedRels


def isSubRelOf(potentialSubRel, relation):
    isSubRel = True
    for relPath in potentialSubRel:
        if not relPath in relation:
            isSubRel = False
            break
    return isSubRel


def allRelsInPathAlgebra(pathAlg):
    """Every relation between every ordered pair of vertices, minimal or not."""
    allRels = []
    vertices = list(pathAlg.vertices())
    for v in vertices:
        for w in vertices:
            allRels.extend(allRelsBetweenVertices(pathAlg, v, w))
    return allRels


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


def isIllegalRelation(pathAlg, relation):
    isIllegal = False
    for n in range(1, len(relation)):
        if relation[n] == []:
            isIllegal = True
            break
        if relation[n][0] != relation[0][0] or relation[n][-1] != relation[0][-1]:
            isIllegal = True
            break
        if not nx.is_path(pathAlg.quiver, relation[n]):
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
        pathAlgebra.printPathAlgebra(pathAlg)
        # input('Press enter to continue...')
    return isIllegal
