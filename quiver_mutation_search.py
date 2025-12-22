import copy
import itertools
import math
import time
import networkx as nx
import numpy as np
import pathAlgebraClass
from quiver_mutation_line import lineQuiverExample, saveLinePathAlgMutation, mutationListLineCleanup, createMutationClassCSV, importMutationClassCSV, saveLineRelationsAndMutationsToCSV, mutationSearchDepthFirst
from quiver_mutation_io import printPathAlgebra, readMutationsFromFile

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
