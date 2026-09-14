"""The pipeline `classification` replaced, kept until its last caller is gone.

These are the older entry points: a per-length search that round-tripped through
text files, the CSV combining steps that followed it, and the "easy relations"
expansion that `classification.expandClassByMoves` supersedes.  They are here
rather than deleted because `main.py` and some by-hand workflows still call
them, and because the text transcripts they write are still the readable record
of a search.

Nothing in the package depends on this module.  See NOTES.md idea 11.
"""

import copy
import glob
import itertools
import math
import time

import networkx as nx

from . import classification
from . import fileio
from . import invariants
from . import lines
from . import mutation
from . import pathAlgebra
from . import paths
from . import search



def findMutationClassesForLine(lineLength, lineName, n_max = 5, manualDFdepth = 1, importAlreadyDoneSearches = False):
    start = time.time()
    DFgrowthFactor = 1 / math.ceil(lineLength / 2)
    mutationClassFiles = []
    allPossibleRelSets = lines.generateAllPossibleLineRelations(lineLength)
    reachedRelSets = []
    allRelsReached = False
    firstUnreachedRelSet = []
    relSetNumber = '0'*(lineLength-2)
    edgeList = []
    for i in range(1, lineLength):
        edgeList.append((i, i+1))
    pathAlg = pathAlgebra.PathAlgebra()
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
        search.mutationSearchDepthFirst(pathAlg, DFdepth, [], quiverName, vertexRelabeling)
        print('\nFirst search for round done! \n')
        mutListDF = fileio.readMutationsFromFile('{0}DF.txt'.format(quiverName))
        mutList = lines.mutationListLineCleanup(mutListDF)
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
                    pathAlgebra.printPathAlgebra(mut[0])
                    DFdepthIncreaseByRound = max( -2, -math.floor(roundCount / math.sqrt(lineLength)))
                    DFdepthIncreaseByLap = 0
                    if roundCount > 4:
                        DFdepthIncreaseByLap = math.floor(n * DFgrowthFactor)
                    DFdepth = math.floor(math.sqrt(lineLength)) + 1 + DFdepthIncreaseByLap + DFdepthIncreaseByRound
                    search.mutationSearchDepthFirst(mut[0], DFdepth, mut[1], quiverName, mut[2])
            lapTimeEnd = time.time()
            print('Lap time for lap {0}: {1} s'.format(n, lapTimeEnd - lapTimeStart))
            lapTimeStart = time.time()
            mutListDF = fileio.readMutationsFromFile('{0}DF.txt'.format(quiverName))
            mutList = lines.mutationListLineCleanup(mutListDF)
            checkedAlready = []
            open('{0}.txt'.format(quiverName), 'w+').close()
            for mut in mutList:
                fileio.saveLinePathAlgMutation(mut[0], mut[1], mut[2], '{0}.txt'.format(quiverName))
                checkedAlready.append(False)
                if mut[0].rels in reachedRelations:
                    checkedAlready[mutList.index(mut)] = True
                else:
                    reachedRelations.append(mut[0].rels)
        fileio.saveLineRelationsToFile('{0}'.format(quiverName))
        mutationClassFiles.append('{0}Relations.txt'.format(quiverName))
        reachedRelSetsThisRound = fileio.readRelationsFromFile('{0}Relations.txt'.format(quiverName))
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
        listOfRelSetLists.append(fileio.readRelationsFromFile(relFile))
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
        listOfRelSetLists.append(fileio.readRelationsFromFile('{0}'.format(relFile)))
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
    allPossibleRelSets = lines.generateAllPossibleLineRelations(lineLength)
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
        oneMutList = fileio.readMutationsFromFile(name)
        for mut in oneMutList:
            with open('{0}Lines.txt'.format(lineName), "a") as f:
                f.write('Mutations: {0}\n'.format(mut[1]))
                f.write('Numbering: {0}\n'.format(mut[2]))
                f.write('Vertices: {0}\n'.format(mut[0].vertices()))
                f.write('Arrows: {0}\n'.format(mut[0].arrows()))
                f.write('Relations: {0}\n'.format(mut[0].rels))
                f.write('-\n')
    reachedAndMissingRelations = collectMutationClasses(lineLength)
    mutListNotReduced = fileio.readMutationsFromFile('{0}Lines.txt'.format(lineName))
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


def combineMutationClasses(lineLength):
    mutClasses = fileio.readMutationClassesFromFile(lineLength, 'A{0}MutationClasses.txt'.format(lineLength))
    print(len(mutClasses))
    coxPolsForClasses = []
    for clas in mutClasses:
        coxPol = invariants.coxeterPoly(clas[0])
        coxPolsForClasses.append(coxPol)
        for quiv in clas[1:]:
            if invariants.coxeterPoly(quiv) != coxPol:
                print('DANGER! Coxeter polynomial does not match!')
                print('Class coxeter polynomial: ', coxPol)
                print('Quiver coxeter polynomial: ', invariants.coxeterPoly(quiv))
                pathAlgebra.printPathAlgebra(quiv)
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
                    search.mutationSearchDepthFirst(quiv, 6, [], 'A{0}TempFile'.format(lineLength), [])
                mutListDF = fileio.readMutationsFromFile('A{0}TempFileDF.txt'.format(lineLength))
                mutList = lines.mutationListLineCleanup(mutListDF)
                isSameClass = False
                for mut in mutList:
                    quiversInThisClass.append(mut[0])
                    if mut[0] in quiversInThisCombinedClass:
                        isSameClass = True
                if isSameClass:
                    for quiv in quiversInThisClass:
                        if not quiv in quiversInThisCombinedClass:
                            quiversInThisCombinedClass.append(quiv)


def combineMutationClassesInCSVfile(filename, mutationDepth, lineLength):
    oldCSVdata = fileio.importMutationClassCSV(filename)
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

    mutationClassCSV = fileio.importMutationClassCSV(filename)
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
                    pathAlg = lines.lineQuiverExample(lineLength, lineRelList, currentRow[2])
                    pathAlgebra.printPathAlgebra(pathAlg)
                    open('{0}DF.txt'.format(quiverName), 'w').close()
                    search.mutationSearchDepthFirst(pathAlg, mutationDepth, [], quiverName)
                    mutList = fileio.readMutationsFromFile('{0}DF.txt'.format(quiverName))
                    cleanMutList = lines.mutationListLineCleanup(mutList)
                    open('{0}.txt'.format(quiverName), 'w').close()
                    for mut in cleanMutList:
                        print('Mutations: {0}'.format(mut[1]))
                        print('Relations: {0}'.format(mut[0].rels))
                        fileio.saveLinePathAlgMutation(mut[0], mut[1], mut[2], '{0}.txt'.format(quiverName))
                    print('csv data: ', mutationClassCSV)
                    mutationClassCSV = classification.saveLineRelationsAndMutationsToCSV('A_{0}_mutation_classes.csv'.format(lineLength), cleanMutList, mutationClassCSV, lineNumberString)
                    if mutationClassCSV[j][1] == baseClass:
                        baseQuiverRelations = [int(rel) in baseClass.split('')]
                        baseQuiver = lines.lineQuiverExample(len(baseQuiverRelations) + 2, baseQuiverRelations)
                        quiverWithRightNumbering = mutation.quiverMutationAtVertices(baseQuiver, )
                        mutDictKeys = list()
                        reverseMutationVertices = []
                        for k in range(numberOfCSVrows):
                            additionalRow = mutationClassCSV[k]
                            additionalClass = additionalRow[1]
                            if additionalClass == currentClass:
                                mutationClassCSV[k][1] = baseClass
                                mutationClassCSV[k][2] = mutationClassCSV[j][2] + ';'
    return newCSVdata


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
    relevantPowerSet = paths.powerset(possible2RelPositions)
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
        pathAlg = lines.lineQuiverExample(lineLength, lineRelList)
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
        pathAlg = lines.lineQuiverExample(lineLength, lineRelListToAdd)
        if not (pathAlg, [], {}) in expandingRelPairMutList:
            expandingRelPairMutList.append((pathAlg, [], {}))
    return expandingRelPairMutList


def expandAllClassesWithEasyRels(lineLength, startRow = 0):
    mutationClassCSV = fileio.importMutationClassCSV('A_{0}_mutation_classes.csv'.format(lineLength))
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
                mutationClassCSV = classification.saveLineRelationsAndMutationsToCSV('A_{0}_mutation_classes.csv'.format(lineLength),
                                                                  mutListWith2RelsAndPairs, mutationClassCSV, lineNumberString)
    return
