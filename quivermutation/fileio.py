"""Reading and writing: the CSV tables and the older text transcripts.

The classification table is a CSV with a parquet alongside it
(`mutationClassTable`).  The text transcripts predate it -- a search used to
write `<name>DF.txt` and the caller parsed it back by string slicing -- and are
kept only for inspecting a search by hand; the pipeline collects in memory now.
"""

import ast
import copy
import csv

import networkx as nx

from . import lines
from . import pathAlgebra



def saveLinePathAlgMutation(pathAlg, mutationVertices = None, vertexRelabeling = None, fileName = 'lineQuiver.txt'):
    mutationVertices = [] if mutationVertices is None else mutationVertices
    vertexRelabeling = {} if vertexRelabeling is None else vertexRelabeling
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
    rPathAlg = pathAlgebra.PathAlgebra()
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
    combinedMutationListClean = lines.mutationListLineCleanup(combinedMutationList)
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


def createMutationClassCSV(lineLength):
    csvData = []
    allLineRels = lines.generateAllPossibleLineRelations(lineLength)
    relSetsAsStrings = []
    for relSet in allLineRels:
        relSetsAsStrings.append(lines.relSetToString(relSet))
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


def saveQuipusToCSV( quipusList, fileName, overwrirteFile = True ):
    if overwrirteFile:
        open(fileName, 'w+').close()
    with open(fileName, 'a') as f:
        for Q in quipusList:
            f.write('{0}\n'.format(str(Q.edges)))
        f.close()
    return


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
                quiverList.append(lines.makeStandardLineQuiver(lineLength, rRels))
            line = f.readline()
    return mutationClasses
