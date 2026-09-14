"""The linear quiver 1 -> 2 -> ... -> n, and the names its algebras go by.

Linearly oriented Nakayama algebras are the objects being classified, and they
get written three ways throughout the repo -- per-vertex relation lengths
`[2,2,3,0,0]`, the class name `'22300'`, and the relation string
`'1;2;3|2;3;4|3;4;5;6'`.  The conversions live here, along with building the
algebra, enumerating all of them for a length, and the standard renumbering that
recognises a mutated quiver as a line again.

`nakayama.LinearNakayamaAlgebra` is the object-oriented face of the same thing
and is what new code should use.
"""

import copy

import networkx as nx

from . import pathAlgebra



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


def makeStandardLineQuiver(lineLength, relationSet):
    pathAlg = pathAlgebra.PathAlgebra()
    pathAlg.add_vertices_from(range(1, lineLength + 1))
    pathAlg.add_arrows_from([[i, i+1] for i in range(1, lineLength)])
    pathAlg.add_rels_from(relationSet)
    return pathAlg


def generateAllLineQuiversWithRelations(lineLength):
    allPossibleLineRelations = generateAllPossibleLineRelations(lineLength)
    allLineQuiversWithRelations = []
    for i in range(len(allPossibleLineRelations)):
        lineQuiver = makeStandardLineQuiver(lineLength, allPossibleLineRelations[i])
        allLineQuiversWithRelations.append(lineQuiver)
    return allLineQuiversWithRelations


def lineQuiverExample(lineLength, relationList, vertexRelabeling = None):
    vertexRelabeling = {} if vertexRelabeling is None else vertexRelabeling
    pathAlg = pathAlgebra.PathAlgebra()
    if len(relationList) != lineLength - 2:
        print('len(relationList) = ', len(relationList))
        print(lineLength - 2)
        print('Error! Invalid relation set.')
    elif bool(relationList) and (max(relationList) > lineLength - 1):
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


def relSetToString(relSet):
    stringList = []
    for i in range(len(relSet)):
        stringInts = [str(int) for int in relSet[i][0]]
        stringOfInts = ";".join(stringInts)
        stringList.append(stringOfInts)
    joinedString = "|".join(stringList)
    return joinedString


def relationStringToLineRelLengths(lineLength, relationString):
    """'1;2;3|3;4;5;6' -> [2, 0, 3, 0] for a line of the given length.

    Entry i is the number of arrows in the relation starting at vertex i + 1.
    """
    lineRelList = [0] * (lineLength - 2)
    for pathString in relationString.split('|'):
        vertexStrings = pathString.split(';')
        if not bool(vertexStrings[0]):
            continue
        path = [int(v) for v in vertexStrings]
        lineRelList[path[0] - 1] = len(path) - 1
    return lineRelList


def lineRelLengthsToClassName(lineRelLengths):
    return ''.join(str(n) for n in lineRelLengths)


def relabelLineAlgebra(pathAlg, currentRelabeling = None):
    currentRelabeling = {} if currentRelabeling is None else dict(currentRelabeling)
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


def convertLineFromCSVnotation( lineLength, lineInCSVnotation ):
    #lineInCSVnotation'1;2;3;4;5|2;3;4;5;6|4;5;6;7;8;9;10'
    listOfRels = []
    listOfRelsAsStr = lineInCSVnotation.split('|')
    for relString in listOfRelsAsStr:
        relAsList = [ int(i) for i in relString.split(';') ]
        listOfRels.append([relAsList])
    lineQuiver = makeStandardLineQuiver(lineLength, listOfRels)
    return lineQuiver


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
                    pathAlgebra.printPathAlgebra(mut[0])
                (pathAlg, newVertexRelabeling) = relabelLineAlgebra(mut[0], vertexRelabeling)
                if printOutput:
                    print('Renumbering: ', newVertexRelabeling)
                    pathAlgebra.printPathAlgebra(pathAlg)
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
