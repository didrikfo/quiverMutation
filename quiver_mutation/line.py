import copy
import csv
import glob
import itertools
import math
import time
import networkx as nx
import numpy as np
from . import path_algebra_class
from .coxeter_ploynomial import coxeter_poly
from .core import quiver_mutation_at_vertex, quiver_mutation_at_vertices
from .csv_io import print_path_algebra, read_mutations_from_file, read_relations_from_file
from .relations import (
    all_rels_between_vertices,
    is_illegal_relation,
    non_minimal_out_rels,
    number_of_paths_up_to_rels,
    reduce_path_algebra,
)
from .utils import rel_set_to_string

def mutation_search_depth_first(pathAlg, depth, mutationVertices = [], quiverName = 'quiver', vertexRelabeling = {}, printOutput = True):
    vertices = list(pathAlg.vertices())
    baseQuiver = copy.deepcopy(pathAlg.quiver)
    quiverAtThisDepth = copy.deepcopy(pathAlg.quiver)
    rels = copy.deepcopy(pathAlg.rels)
    allRels = []
    for v in vertices:
        for w in vertices:
            allRels.extend(all_rels_between_vertices(pathAlg, v, w))
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
        print_path_algebra(pathAlg)
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
                    numberOfPathsToVertexUpToRels = number_of_paths_up_to_rels(pathAlg, v, vertex)
                    for w in vertexImmideateSuccessors:
                        if numberOfPathsToVertexUpToRels > number_of_paths_up_to_rels(pathAlg, v, w):
                            mutationPossible = False
                            break
                    if not mutationPossible:
                        break
                # for v in vertexPredecessors:
                #     for nonMinRel in non_minimal_out_rels(pathAlg, v):
                #         if len(nonMinRel) == 1:
                #             if nonMinRel[0][-1] in vertexImmideateSuccessors and nonMinRel[0][-2] != vertex:
                #                 mutationPossible = False
                #                 break
                #     if not mutationPossible:
                #         break
            if mutationPossible:
                mutationVerticesAtDepth.append(vertexRelabeling[vertex])
                mutPathAlg = quiver_mutation_at_vertex(pathAlg, vertex)
                for rel in mutPathAlg.rels:
                    if is_illegal_relation(mutPathAlg, rel):
                        discardMutation = True
                        break
                if discardMutation:
                    break
                mutPathAlg = reduce_path_algebra(mutPathAlg)
                mutation_search_depth_first(copy.deepcopy(mutPathAlg), depth, mutationVerticesAtDepth, quiverName, vertexRelabeling, printOutput)
    return

def relabel_line_algebra(pathAlg, currentRelabeling = {}):
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

def save_line_path_alg_mutation(pathAlg, mutationVertices = [], vertexRelabeling = {}, fileName = 'lineQuiver.txt'):
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

def mutation_list_line_cleanup(mutationList, relabelNodes = True, printOutput = True):
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
                    print_path_algebra(mut[0])
                (pathAlg, newVertexRelabeling) = relabel_line_algebra(mut[0], vertexRelabeling)
                if printOutput:
                    print('Renumbering: ', newVertexRelabeling)
                    print_path_algebra(pathAlg)
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

def mutation_list_line_cleanup_keep_dupes(mutationList, relabelNodes = True, discardLongerDupes = False):
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
                (pathAlg, newVertexRelabeling) = relabel_line_algebra(mut[0], vertexRelabeling)
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

def relation_dual_line_quiver(quiverWrels):
    numberOfVertices = len(quiverWrels['quiver'])
    rels = quiverWrels['rels']
    dualRels = nx.MultiDiGraph()
    dualRels.add_nodes_from(rels.nodes)
    for rel in rels.edges:
        dualRels.add_edge(numberOfVertices - rel[1] + 1, numberOfVertices - rel[0] + 1)
    dualQuiverWrels = {'quiver' : quiverWrels['quiver'], 'rels' : dualRels}
    return dualQuiverWrels

def is_relation_dual_line_quiver(quiverWrels1, quiverWrels2):
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

def save_line_relations_to_file(fileName):
    mutationList = read_mutations_from_file('{0}.txt'.format(fileName))
    open('{0}Relations.txt'.format(fileName), 'w+').close()
    for mut in mutationList:
        relations = mut[0].rels
        with open('{0}Relations.txt'.format(fileName), "a") as f:
            f.write('{0}\n'.format(relations))
            f.close()
    return

def save_line_relations_and_mutations_to_file(fileName, saveNumbering = False):
    mutationList = read_mutations_from_file('{0}.txt'.format(fileName))
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

def generate_list_of_relations(listOfFileNames, combinedFileName = 'allRelations'):
    combinedMutationList = []
    for fileName in listOfFileNames:
        save_line_relations_to_file(fileName)
        mutationList = read_mutations_from_file(fileName + '.txt')
        combinedMutationList.extend(mutationList[:])
    combinedMutationListClean = mutation_list_line_cleanup(combinedMutationList)
    open('{0}.txt'.format(combinedFileName), 'w+').close()
    for mut in combinedMutationListClean:
        save_line_path_alg_mutation(mut[0], mut[1], mut[2], combinedFileName + '.txt')
    save_line_relations_to_file(combinedFileName)
    return

def generate_all_possible_line_relations(lineLength):
    lineStart = 1
    lineStop = lineLength
    if lineLength <= 2:
        return [[]]
    elif lineLength == 2:
        return [[], [[[*range(lineStart, lineStop + 1)]]]]
    relSetList = generate_all_possible_line_relations(lineLength - 1)
    allPossibleRelSets = relSetList[:]
    lastRelStart = 0
    for relSet in relSetList:
        if bool(relSet):
            lastRelStart = relSet[-1][0][0]
        for i in range(lastRelStart + 1, lineStop - 1):
            allPossibleRelSets.append(relSet + [[[*range(i,lineStop + 1)]]])
    allPossibleRelSets.sort()
    return allPossibleRelSets

def collect_mutation_classes(lineLength, saveToFile = False, printOutput = False):
    lineName = 'A{0}'.format(lineLength)
    open('{0}MutationClasses.txt'.format(lineName), 'w+').close()
    mutationClassFiles = []
    for name in glob.glob('{0}_*Relations.txt'.format(lineName)):
        mutationClassFiles.append(name)
    mutationClassFiles.sort()
    listOfRelSetLists = []
    for relFile in mutationClassFiles:
        listOfRelSetLists.append(read_relations_from_file('{0}'.format(relFile)))
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
    allPossibleRelSets = generate_all_possible_line_relations(lineLength)
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

def combine_line_mutation_files(lineLength):
    lineName = 'A{0}'.format(lineLength)
    open('{0}Lines.txt'.format(lineName), 'w+').close()
    mutationClassFiles = []
    for name in glob.glob('{0}_*[0-2].txt'.format(lineName)):
        print(name)
        mutationClassFiles.append(name)
        oneMutList = read_mutations_from_file(name)
        for mut in oneMutList:
            with open('{0}Lines.txt'.format(lineName), "a") as f:
                f.write('Mutations: {0}\n'.format(mut[1]))
                f.write('Numbering: {0}\n'.format(mut[2]))
                f.write('Vertices: {0}\n'.format(mut[0].vertices()))
                f.write('Arrows: {0}\n'.format(mut[0].arrows()))
                f.write('Relations: {0}\n'.format(mut[0].rels))
                f.write('-\n')
    reachedAndMissingRelations = collect_mutation_classes(lineLength)
    mutListNotReduced = read_mutations_from_file('{0}Lines.txt'.format(lineName))
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

def make_standard_line_quiver(lineLength, relationSet):
    pathAlg = path_algebra_class.PathAlgebra()
    pathAlg.add_vertices_from(range(1, lineLength + 1))
    pathAlg.add_arrows_from([[i, i+1] for i in range(1, lineLength)])
    pathAlg.add_rels_from(relationSet)
    return pathAlg

def read_mutation_classes_from_file(lineLength, fileName):
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
                quiverList.append(make_standard_line_quiver(lineLength, rRels))
            line = f.readline()
    return mutationClasses

def combine_mutation_classes(lineLength):
    mutClasses = read_mutation_classes_from_file(lineLength, 'A{0}MutationClasses.txt'.format(lineLength))
    print(len(mutClasses))
    coxPolsForClasses = []
    for clas in mutClasses:
        coxPol = coxeter_poly(clas[0])
        coxPolsForClasses.append(coxPol)
        for quiv in clas[1:]:
            if coxeter_poly(quiv) != coxPol:
                print('DANGER! Coxeter polynomial does not match!')
                print('Class coxeter polynomial: ', coxPol)
                print('Quiver coxeter polynomial: ', coxeter_poly(quiv))
                print_path_algebra(quiv)
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
                    mutation_search_depth_first(quiv, 6, [], 'A{0}TempFile'.format(lineLength), [])
                mutListDF = read_mutations_from_file('A{0}TempFileDF.txt'.format(lineLength))
                mutList = mutation_list_line_cleanup(mutListDF)
                isSameClass = False
                for mut in mutList:
                    quiversInThisClass.append(mut[0])
                    if mut[0] in quiversInThisCombinedClass:
                        isSameClass = True
                if isSameClass:
                    for quiv in quiversInThisClass:
                        if not quiv in quiversInThisCombinedClass:
                            quiversInThisCombinedClass.append(quiv)

def generate_all_line_quivers_with_relations(lineLength):
    allPossibleLineRelations = generate_all_possible_line_relations(lineLength)
    allLineQuiversWithRelations = []
    for i in range(len(allPossibleLineRelations)):
        lineQuiver = make_standard_line_quiver(lineLength, allPossibleLineRelations[i])
        allLineQuiversWithRelations.append(lineQuiver)
    return allLineQuiversWithRelations

def line_quiver_example(lineLength, relationList, vertexRelabeling = {}):
    pathAlg = path_algebra_class.PathAlgebra()
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

def create_mutation_class_csv(lineLength):
    csvData = []
    allLineRels = generate_all_possible_line_relations(lineLength)
    relSetsAsStrings = []
    for relSet in allLineRels:
        relSetsAsStrings.append(rel_set_to_string(relSet))
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

def import_mutation_class_csv(filename):
    csvData = []
    with open(filename, newline='') as csvfile:
        csvReader = csv.reader(csvfile, delimiter=',')
        csvData = list(csvReader)
    #    for row in csvReader:
    #        csvData.append(', '.join(row))
    #for row in csvData:
    return csvData

def save_line_relations_and_mutations_to_csv(fileName, mutationList, csvData, mutationClassName, printOutput = False):
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
        coxPoly = coxeter_poly(mut[0])
        relSet = mut[0].rels
        relSetString = rel_set_to_string(relSet)
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
                    reverseMutationVertices = reverse_mutation_sequence(mutationVertices, vertexNumbering)
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
                    baseVertexNumbering[n] = oldNumberingList[get_vertex_numbering_key_from_value(vertexNumbering, n) - 1]
                minMutationLength = mutationLength
    for mut in mutationList:
        localVertexNumbering = mut[2]
        localVertexNumberingList = [localVertexNumbering[n] for n in range(1, len(localVertexNumbering) + 1)]
        vertexNumberingList = [baseVertexNumbering[localVertexNumbering[n]] for n in range(1, len(baseVertexNumbering) + 1)]
        vertexNumberingString = ';'.join([str(n) for n in vertexNumberingList])
        relSet = mut[0].rels
        relSetString = rel_set_to_string(relSet)
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

def combine_mutation_classes_in_csvfile(filename, mutationDepth, lineLength):
    oldCSVdata = import_mutation_class_csv(filename)
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

    mutationClassCSV = import_mutation_class_csv(filename)
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
                    pathAlg = line_quiver_example(lineLength, lineRelList, currentRow[2])
                    print_path_algebra(pathAlg)
                    open('{0}DF.txt'.format(quiverName), 'w').close()
                    mutation_search_depth_first(pathAlg, mutationDepth, [], quiverName)
                    mutList = read_mutations_from_file('{0}DF.txt'.format(quiverName))
                    cleanMutList = mutation_list_line_cleanup(mutList)
                    open('{0}.txt'.format(quiverName), 'w').close()
                    for mut in cleanMutList:
                        print('Mutations: {0}'.format(mut[1]))
                        print('Relations: {0}'.format(mut[0].rels))
                        save_line_path_alg_mutation(mut[0], mut[1], mut[2], '{0}.txt'.format(quiverName))
                    print('csv data: ', mutationClassCSV)
                    mutationClassCSV = save_line_relations_and_mutations_to_csv('A_{0}_mutation_classes.csv'.format(lineLength), cleanMutList, mutationClassCSV, lineNumberString)
                    if mutationClassCSV[j][1] == baseClass:
                        baseQuiverRelations = [int(rel) in baseClass.split('')]
                        baseQuiver = line_quiver_example(len(baseQuiverRelations) + 2, baseQuiverRelations)
                        quiverWithRightNumbering = quiver_mutation_at_vertices(baseQuiver, )
                        mutDictKeys = list()
                        reverseMutationVertices = []
                        for k in range(numberOfCSVrows):
                            additionalRow = mutationClassCSV[k]
                            additionalClass = additionalRow[1]
                            if additionalClass == currentClass:
                                mutationClassCSV[k][1] = baseClass
                                mutationClassCSV[k][2] = mutationClassCSV[j][2] + ';'
    return newCSVdata

def convert_line_from_csvnotation( lineLength, lineInCSVnotation ):
    #lineInCSVnotation'1;2;3;4;5|2;3;4;5;6|4;5;6;7;8;9;10'
    listOfRels = []
    listOfRelsAsStr = lineInCSVnotation.split('|')
    for relString in listOfRelsAsStr:
        relAsList = [ int(i) for i in relString.split(';') ]
        listOfRels.append([relAsList])
    lineQuiver = make_standard_line_quiver(lineLength, listOfRels)
    return lineQuiver
