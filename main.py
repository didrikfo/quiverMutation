from itertools import chain, combinations
import itertools
import networkx as nx
import matplotlib.pyplot as plt
import math
import time
import os
import operator
import random
import copy

import numpy as np

import pathAlgebraClass
import glob
import csv
from quiverExamples import *
from quiverMutation import *
from datetime import datetime
import sympy
from sympy.interactive.printing import init_printing
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
from sympy.abc import x, y
import nodepy
from nodepy import rooted_trees as rt

now = datetime.now() # current date and time

mode = 3

# mode 1: Mutation search
# mode 2: investigate mutation sequence
# mode 3: generate trees

startTime = time.time()


#PurePoly(lambda**12 + lambda**11 - 2*lambda**10 - 4*lambda**9 - 4*lambda**8 - 4*lambda**7 - 4*lambda**6 - 4*lambda**5 - 4*lambda**4 - 4*lambda**3 - 2*lambda**2 + lambda + 1, lambda, domain='ZZ')


if mode == 1:
    mutationSearch(lineLength = 12, mutationDepthStart = 2, startRow = 30000)

elif mode == 2: # Investigate mutation sequence
    pathAlg = pathAlgebraClass.PathAlgebra()
    pathAlg.add_paths_from([[1,2,3,4,5,6,7,8]])
    pathAlg.add_rels_from([[[1,2,3]],[[2,3,4]],[[3,4,5]],[[4,5,6]],[[5,6,7,8]]])#
    #pathAlg = convertLineFromCSVnotation(len(pathAlg.vertices()), '1;2;3;4;5;6|3;4;5;6;7;8|5;6;7;8;9;10')

    mutationVertexList = [1,1,-4,-6,-6,-6,-8,-4,-8,5,7,4,3,1,1,2] #
    l = len(mutationVertexList)
    firstDisplayStep = 0 # 0 if you want to display all mutation steps, len(mutationVertexList) if you only want the last
    #pathAlg = onePointExtension(pathAlg, [11,5], [[[11,5,6]]])
    pathAlg = quiverMutation(pathAlg, mutationVertexList, firstDisplayStep)

    #pathAlg = onePointExtension(pathAlg, [11, 13], [[[5,6,7,12,1,8,9,10,11,13]]])
    #mutationVertexList = [-13, -13, -11, -11]  #
    #quiverMutation(pathAlg, mutationVertexList, firstDisplayStep)

    if False:
        baseCoxPol = coxeterPoly(pathAlg)
        print(coxeterPoly(pathAlg))
        printPathAlgebra(pathAlg)
        if firstDisplayStep == 0:
            plotQuiver(pathAlg)
        for i in range(len(mutationVertexList)):
            if mutationVertexList[i] >= 0:
                pathAlg = quiverMutationAtVertex(pathAlg, mutationVertexList[i])
            else:
                pathAlg = leftQuiverMutationAtVertex(pathAlg, -mutationVertexList[i])
            pathAlg = reducePathAlgebra(pathAlg)
            currentCoxPol = coxeterPoly(pathAlg)
            cartMat = cartanMatrix(pathAlg)
            print(np.matrix(cartMat))
            if currentCoxPol != baseCoxPol:
                print('COXETER POLYNOMIAL HAS CHANGED!')
            print('Mutations: ', mutationVertexList[0:i+1])
            print(currentCoxPol)
            printPathAlgebra(pathAlg)
            if i+1 >= firstDisplayStep:
                plotQuiver(pathAlg)
elif mode == 3:
    #startTime = time.time()
    #number_of_quipus = count_quipus(18)
    #for i in range(2,len(number_of_quipus)):
    #    print(number_of_quipus[i], 'quipus of length', i )
    #print('Total runtime', time.time() - startTime)
    print('Code started running at', time.strftime('%H:%M:%S',  time.gmtime(time.time()+7200)))
    for length in range(2,35):
        startTime = time.time()
        #quipus = generateAllQuipusGPT(length)
        numberOfQuipusOfLength = count_quipusV1(length)
        endTime = time.time()
        print(numberOfQuipusOfLength, 'quipus of length', length,  'generated in', endTime - startTime, 's')


    if False:
        displayTrees = False

        h1quipus = generateAllHeightOneQuipus(17)

        print(len(h1quipus), ' length one quipus in ', time.time() - startTime, 's')

        saveQuipusToCSV(h1quipus, 'quipusCSVTest', overwrirteFile=True)

        i = 1
        s = 1
        while s > 0:
            s = sum(1 for Q in h1quipus if len(Q) == i)
            print(s, 'height one quipus of length ', i, ' in ', time.time() - startTime, 's')
            i += 1


        quipusLists = generateAllQuipusUpToLength(40)
        open('test_quipu_file.txt', 'w+').close()
        with open('test_quipu_file.txt', "a") as f:
            for quipuList in quipusLists:
                f.write('Edges of quipus of length {0}:\n'.format(len(quipuList[0])))
                for quipu in quipuList:
                    f.write('{0}\n'.format(quipu.edges))
                f.write('\n')
            f.close()

    #    for quipuList in quipusLists:
    #        print('Quipus of length ', len(quipuList[0].nodes))
    #        for quipu in quipuList:
    #            print(quipu.edges)
    #            nx.draw(quipu)
    #            plt.show()


        for n in range(23, 40):
            if False:
                trees = list(nx.generators.nonisomorphic_trees(n))
                quipus = []
                for tree in trees:
                        maxDeg = 0
                        for deg in tree.degree:
                            if deg[1] > maxDeg:
                                maxDeg = deg[1]
                        if maxDeg <= 3:
                            deg3Vertices = []
                            for v in tree.nodes:
                                if tree.degree(v) == 3:
                                    deg3Vertices.append(v)
                            longestDeg3Path = []
                            for v1 in deg3Vertices:
                                for v2 in deg3Vertices[deg3Vertices.index(v1):]:
                                    pathBetween = nx.shortest_path(tree, v1,v2)
                                    if len(pathBetween) > len(longestDeg3Path):
                                        longestDeg3Path = pathBetween
                            isQuipu = True
                            for v in deg3Vertices:
                                if not v in longestDeg3Path:
                                    isQuipu = False
                                    break
                            if isQuipu:
                                quipus.append(tree)
    ###
            if False:
                trees = list(nx.generators.nonisomorphic_trees(n))
                quipus = []
                for tree in trees:
                    isQuipu = False
                    maxDeg = 0
                    for deg in tree.degree:
                        if deg[1] > maxDeg:
                            maxDeg = deg[1]
                    if maxDeg <= 3:
                        isQuipu = True
                        deg3SubGraph = nx.Graph()
                        deg3Vertices = []
                        for v in tree.nodes:
                            if tree.degree(v) == 3:
                                deg3Vertices.append(v)
                        if len(deg3Vertices) > 3:
                            nx.add_path(deg3SubGraph, nx.shortest_path(tree,deg3Vertices[0],deg3Vertices[1]))
                            for v1 in deg3Vertices[2:]:
                                pathToAddAsGraph = nx.Graph()
                                nx.add_path(pathToAddAsGraph, nx.shortest_path(tree, deg3Vertices[0], v1))
                                deg3SubGraph = nx.compose(deg3SubGraph, pathToAddAsGraph)
                            correspondingPathGraph = nx.path_graph(len(deg3SubGraph))
                            if not nx.is_isomorphic(deg3SubGraph, correspondingPathGraph):
                                isQuipu = False
                    if isQuipu:
                        quipus.append(tree)
            #numberOfTrees = len(quipus)
            #print('Generated ', len(quipus), ' quipus of length ', n)
    ###
            if displayTrees:
                nRows = int(np.ceil(np.sqrt(numberOfTrees) * np.sqrt(2))) - 1
                nCols = int(np.floor(np.sqrt(numberOfTrees) * (1 / np.sqrt(2)))) + 1
                # nRows = numberOfTrees
                fig, axes = plt.subplots(nrows=nRows, ncols=nCols)
                ax = axes.flatten()
                for i in range(nRows * nCols):
                    if i < numberOfTrees:
                        nx.draw_networkx(quipus[i], ax=ax[i], node_size=10, with_labels=False)
                        coxPoly = coxPolyOfTree(trees[i])
                        # ax[i].title.set_text(coxPoly.as_expr())
                    ax[i].set_axis_off()
                plt.subplot_tool()
                plt.show()
elif mode == 4: # Investigate mutation sequence
    pathAlg = pathAlgebraClass.PathAlgebra()
    pathAlg.add_paths_from([[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]])
    pathAlg.add_rels_from([[[3,4,5,6,7,8,9,10,11,12]],[[8,9,10,11,12,13,14,15,16]],[[15,16,17]],[[16,17,18]]])#
    #pathAlg = onePointExtension(pathAlg, [18,21], [[[7,8,9,10,11,12,13,14,15,16,17,18,21]]])
    mutationVertexList = [8,8,9,9,10,10]
    l = len(mutationVertexList)
    firstDisplayStep = l-2 # 0 to display all mutation steps, len(mutationVertexList) only the last, np.infty none
    pathAlg = quiverMutation(pathAlg, mutationVertexList, firstDisplayStep)

    #pathAlg = onePointExtension(pathAlg, [17,19], [[[5,6,7,8,16,17,19]]])
    #mutationVertexList = [-19,-19,-17,-17,-18,-18,-10,-11,-12,-13,-14,-15,-16,-16,-16]
    #firstDisplayStep = len(mutationVertexList) # 0 to display all mutation steps, len(mutationVertexList) only the last, np.infty none
    #pathAlg = quiverMutation(pathAlg, mutationVertexList, firstDisplayStep)

#createMutationClassCSV(13)
#expandAllClassesWithEasyRels(12, startRow = 30500)

if False:
    n = 11
    AnCoxPols = generateAllCoxeterPolynomials(n)
    treeMatrices = nx.generators.nonisomorphic_trees(n, 'matrix')
    for matrixAsList in treeMatrices:
        tree = nx.from_numpy_matrix(np.matrix(matrixAsList))
        nx.draw(tree)
        mat = sympy.Matrix(np.triu(Matrix(matrixAsList))) + sympy.eye(n)
        print(np.matrix(mat))
        matInvTrans = mat.inv().transpose()
        coxeterMatrix = -matInvTrans * mat
        coxeterPolynomial = coxeterMatrix.charpoly()
        print(coxeterPolynomial.as_expr())
        print()
        plt.show()




if False:
    pathAlg = lineQuiverExample(9, [2, 3, 4, 4, 0, 3, 0])
    quiverName = 'A9_2344030'
    open('{0}DF.txt'.format(quiverName), 'w').close()
    mutationSearchDepthFirst(pathAlg, 6, [], quiverName)
    mutList = readMutationsFromFile('{0}DF.txt'.format(quiverName))
    cleanMutList = mutationListLineCleanup(mutList)
    open('{0}.txt'.format(quiverName), 'w').close()
    for mut in cleanMutList:
        print('Mutations: {0}'.format(mut[1]))
        print('Relations: {0}'.format(mut[0].rels))
        saveLinePathAlgMutation(mut[0], mut[1], mut[2], '{0}.txt'.format(quiverName))



if False:
    quiverName = 'A6_3300test'
    pathAlg = lineQuiverExample(6, [3,3,0,0], {1: 3, 2: 4, 3: 5, 4: 6, 5: 2, 6: 1})
    open('{0}DF.txt'.format(quiverName), 'w').close()
    mutationSearchDepthFirst(pathAlg, 8, [], quiverName)
    mutList1 = readMutationsFromFile('{0}DF.txt'.format(quiverName))
    cleanMutList1 = mutationListLineCleanup(mutList1)
    open('{0}.txt'.format(quiverName), 'w').close()
    for mut in cleanMutList1:
        print('Mutations: {0}'.format(mut[1]))
        print('Relations: {0}'.format(mut[0].rels))
        saveLinePathAlgMutation(mut[0], mut[1], mut[2], '{0}.txt'.format(quiverName))
    saveLineRelationsAndMutationsToFile(quiverName, saveNumbering=True)

if False:
    mutClasses = readMutationClassesFromFile(7, 'A7MutationClasses.txt')
    print('Number of mutation classes: ', len(mutClasses))
    for mutClass in mutClasses:
        for pathAlg in mutClass:
            coxPoly = coxeterPoly(pathAlg)
            print(coxPoly, pathAlg.rels)
        print()


if False:
    if mode == 1: # Mutation search
        lineLength = 7
        nMax = 10
        DFdepth = 4
        quiverName = 'A' + str(lineLength)
        findMutationClassesForLine(lineLength, quiverName, nMax, DFdepth, True) #(lineLength, name, n_max, importAlreadyDoneSearches)
        generateAllCoxeterPolynomials(lineLength)
    elif mode == 2: # Debug mutation search
        pathAlg = A7_40030
        print([])
        baseCoxPol = coxeterPoly(pathAlg)
        print(coxeterPoly(pathAlg))
        printPathAlgebra(pathAlg)
        mutationVertexList = [1, 4, 2, 3, 2, 1, 5, 2]
        for i in range(len(mutationVertexList)):
            pathAlg = quiverMutationAtVertex(pathAlg, mutationVertexList[i])
            pathAlg = reducePathAlgebra(pathAlg)
            currentCoxPol = coxeterPoly(pathAlg)
            cartMat = cartanMatrix(pathAlg)
            print(np.matrix(cartMat))
            if currentCoxPol != baseCoxPol:
                print('COXETER POLYNOMIAL HAS CHANGED!')
            print('Mutations: ', mutationVertexList[0:i+1])
            print(currentCoxPol)
            printPathAlgebra(pathAlg)
    elif mode == 3: # Investigate mutation class
        lineLength = len(mutClasses[0][0].vertices())
        for mutClass in mutClasses:
            for pathAlg in mutClass:
                relSetNumberAsList = ['0']*(lineLength - 2)
                for rel in pathAlg.rels:
                    relSetNumberAsList[rel[0][0] - 1] = str(rel[0][-1] - rel[0][0])
                relSetNumber = ''.join(relSetNumberAsList)
                pathAlgName = 'A{0}_{1}'.format(lineLength, relSetNumber)
                open('{0}.txt'.format(pathAlgName), 'w+').close()
                mutationSearchDepthFirst(pathAlg, 8, [], pathAlgName)
                mutListDF = readMutationsFromFile('{0}DF.txt'.format(pathAlgName))
                mutList = mutationListLineCleanup(mutListDF)
                open('{0}.txt'.format(pathAlgName), 'w+').close()
                for mut in mutList:
                    saveLinePathAlgMutation(mut[0], mut[1], mut[2], '{0}.txt'.format(pathAlgName))
                saveLineRelationsToFile(pathAlgName)




# mutationSearchDepthFirst(A5_0, 6, [], 'A5_0')
# mutListDF = readMutationsFromFile('{0}DF.txt'.format('A5_0'))
# open('{0}M.txt'.format('A5_0'), 'w+').close()
# for mut in mutListDF:
#     print(mut[0].rels)
#     saveLinePathAlgMutation(mut[0], mut[1], mut[2], 'A5_0M.txt')



#A4relations = findMutationClassesForLine(9, 'A9', 9, 6)

if False:
    start = time.time()
    vertexRelabeling = {}
    saveQuiverPlot = False
    quiverName = 'A6_303'
    myQuiverWrels = A6_303
    open('{0}DF.txt'.format(quiverName), 'w+').close()
    mutationSearchDepthFirst(myQuiverWrels, 5, [], quiverName, vertexRelabeling)
    print('\nFirst search done! \n')

    print("//////////////////////////////////////////////////////////////////////////////////////////")


    mutListDF = readMutationsFromFile('{0}DF.txt'.format(quiverName))
    mutList = mutationListLineCleanup(mutListDF)
    n_max = 7
    DFdepth = 3
    lapTimeStart = time.time()
    reachedRelations = []
    for mut in mutList:
        reachedRelations.append(mut[0]['rels'].edges)
    checkedAlready = [False for i in range(len(mutList))]
    for n in range(1, n_max + 1):
        open('{0}DF.txt'.format(quiverName), 'w+').close()
        for mut in mutList:
            if checkedAlready[mutList.index(mut)]:
                longestPathLength = nx.dag_longest_path_length(mut[0]['quiver'])
                with open('{0}DF.txt'.format(quiverName), "a") as f:
                    f.write('Mutations: {0}\n'.format(mut[1]))
                    f.write('Numbering: {0}\n'.format(mut[2]))
                    f.write("Longest path: {0}\n".format(longestPathLength))
                    f.write('Vertices: {0}\n'.format(mut[0]['quiver'].nodes))
                    f.write('Arrows: {0}\n'.format(mut[0]['quiver'].edges))
                    f.write('Relations: {0}\n'.format(mut[0]['rels'].edges))
                    f.write('-\n')
                    f.close()
            else:
                print('Mutations: ', mut[1])
                printQuiver(mut[0])
                mutationSearchDepthFirst(mut[0], DFdepth + math.ceil(n/2), mut[1], quiverName, mut[2])
        lapTimeEnd = time.time()
        print('Lap time for lap {0}: {1} s'.format(n, lapTimeEnd - lapTimeStart))
        lapTimeStart = time.time()
        mutListDF = readMutationsFromFile('{0}DF.txt'.format(quiverName))
        mutList = mutationListLineCleanup(mutListDF)
        checkedAlready = []
        open('{0}L{1}.txt'.format(quiverName, n), 'w+').close()
        for mut in mutList:
            saveLineQuiverMutation(mut[0], mut[1], mut[2], '{0}L{1}.txt'.format(quiverName, n))
            checkedAlready.append(False)
            if mut[0]['rels'].edges in reachedRelations:
                checkedAlready[mutList.index(mut)] = True
            else:
                reachedRelations.append(mut[0]['rels'].edges)

    print('Number of different sets of relations reached: ', len(reachedRelations))



    if saveQuiverPlot:
        mutIndex = random.randint(0, len(mutList))
        mutationList = mutList[mutIndex][1]

        print([])
        printQuiver(myQuiverWrels)
        fileName = quiverName
        saveToFile = True
        showPlot = False
        rootFolder = '/Users/didrikfo/OneDrive - NTNU/PhD/QuiverMutation'
        folder = rootFolder + '/' + fileName + 'm' + ''.join(map(str, mutationList)) + '/'
        plotQuiver(myQuiverWrels, showPlot, saveToFile, fileName, folder)
        fileName = fileName + 'm'
        for i in range(0, len(mutationList)):
            print(mutationList[:i+1])
            myQuiverWrels = quiverMutationAtVertex(myQuiverWrels, mutationList[i])
            myQuiverWrels = reduceQuiver(myQuiverWrels)
            printQuiver(myQuiverWrels)
            fileName = fileName + str(mutationList[i])
            plotQuiver(myQuiverWrels, showPlot, saveToFile, fileName, folder)


#generateListOfRelations(['A6_0L7', 'A6_24L7', 'A6_2mL7', 'A6_23L7', 'A6_303L7'], 'A6')


if False:
    quiverWrels = A8_223000
    print([])
    print(coxeterPoly(quiverWrels))
    printQuiver(quiverWrels)
    mutationVertexList = [7, 6, 1, 1, 1, 1, 2, 2, 2, 6, 3, 3, 5, 4, 5, 4, 4, 3, 8, 8, 3, 3, 7, 1, 2, 2, 1, 7]
    for i in range(len(mutationVertexList)):
        quiverWrels = quiverMutationAtVertex(quiverWrels, mutationVertexList[i])
        quiverWrels = reduceQuiver(quiverWrels)
        print(mutationVertexList[0:i+1])
        print(coxeterPoly(quiverWrels))
        printQuiver(quiverWrels)



# mutationVertexList = [2, 5, 5, 4, 4, 9, 2, 2, 8, 8, 4, 9, 2, 1, 1, 1, 1, 3, 3, 3, 6, 7, 6, 7, 2, 1, 3, 7, 6, 4, 4, 5, 5, 8]
# myQuiverWrels = A9_4
# for v in mutationVertexList:
#     myQuiverWrels = quiverMutationAtVertex(myQuiverWrels, v)
#     myQuiverWrels = reduceQuiver(myQuiverWrels)
#     printQuiver(myQuiverWrels)
# mutList = readMutationsFromFile('A9_4L6.txt')
# newMut = (myQuiverWrels, mutationVertexList, {1:6, 2:7, 3:3, 4:1, 5:2, 6:9, 7:8, 8:5, 9:4})
# mutList.append(newMut)
# mutationList = mutationListLineCleanup(mutList)
# open('A9mutClass1.txt', 'w+').close()
# for mut in mutationList:
#     saveLineQuiverMutation(mut[0], mut[1], mut[2], 'A9mutClass1.txt')

# mutList1 = readMutationsFromFile('A9mutClass1.txt')
# mutList2 = readMutationsFromFile('A9mutClass2.txt')
# for mut1 in mutList1:
#     for mut2 in mutList2:
#         if mut1[0]['rels'].edges == mut2[0]['rels'].edges:
#             print('Mutations: ', mut1[1])
#             printQuiver(mut1[0])


# mutationList = [8, 8, 6, 6, 5, 5, 5, 5] #mutList[randomMutIndex][1]
#
# myQuiverWrels = SA5
# quiverName = 'SA5'
#
# print([])
# printQuiver(myQuiverWrels)
# fileName = quiverName
# saveToFile = True
# showPlot = False
# rootFolder = '/Users/didrikfo/OneDrive - NTNU/PhD/QuiverMutation'
# folder = rootFolder + '/' + fileName + 'm' + ''.join(map(str, mutationList)) + '/'
# plotQuiver(myQuiverWrels, showPlot, saveToFile, fileName, folder)
# fileName = fileName + 'm'
# for i in range(0, len(mutationList)):
#     print(mutationList[:i+1])
#     myQuiverWrels = quiverMutationAtVertex(myQuiverWrels, mutationList[i])
#     myQuiverWrels = reduceQuiver(myQuiverWrels)
#     printQuiver(myQuiverWrels)
#     fileName = fileName + str(mutationList[i])
#     plotQuiver(myQuiverWrels, showPlot, saveToFile, fileName, folder)


# end = time.time()
# print('Runtime: ', end - start, 's')

# mutList = readMutationsFromFile('SA5L2.txt')
# open('SA5DualMutations.txt', 'w+').close()
# for mut1 in mutList:
#     mutIndex = mutList.index(mut1)
#     for mut2 in mutList[mutIndex:]:
#         if isRelationDualLineQuiver(mut1[0], mut2[0]):
#             with open('SA5DualMutations.txt', "a") as f:
#                 f.write('Mutations: {0}\n'.format(mut1[1]))
#                 f.write('Numbering: {0}\n'.format(mut1[2]))
#                 f.write('Vertices: {0}\n'.format(mut1[0]['quiver'].nodes))
#                 f.write('Arrows: {0}\n'.format(mut1[0]['quiver'].edges))
#                 f.write('Relations: {0}\n'.format(mut1[0]['rels'].edges))
#                 f.write('|\n')
#                 f.write('Mutations: {0}\n'.format(mut2[1]))
#                 f.write('Numbering: {0}\n'.format(mut2[2]))
#                 f.write('Vertices: {0}\n'.format(mut2[0]['quiver'].nodes))
#                 f.write('Arrows: {0}\n'.format(mut2[0]['quiver'].edges))
#                 f.write('Relations: {0}\n'.format(mut2[0]['rels'].edges))
#                 f.write('-\n')
#             f.close()