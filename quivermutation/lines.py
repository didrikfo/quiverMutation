"""The linear quiver 1 -> 2 -> ... -> n, and the names its algebras go by.

Linearly oriented Nakayama algebras are the objects being classified, and they
get written three ways throughout the repo -- per-vertex relation lengths
`[2,2,3,0,0]`, the class name `'22300'`, and the relation string
`'1;2;3|2;3;4|3;4;5;6'`.  The conversions live here, along with enumerating
every LNA of a length and the standard renumbering that recognises a mutated
quiver as a line again.

Building one is `nakayama.LinearNakayamaAlgebra`, which is also where these
conversions appear as methods; what is left here is what has no algebra to be a
method of -- a relation set from a mutated quiver, a relation string read out of
the table, a bare list of relation lengths.
"""

import copy

import networkx as nx

from . import pathAlgebra



def generateAllPossibleLineRelations(lineLength):
    """Every admissible relation set on the line of that length, as path sets.

    There are Catalan(lineLength - 1) of them.  Built by recursion on the
    length: every relation set of the shorter line is still one here, and each
    of them also admits a new relation ending at the new last vertex, starting
    anywhere after the previous relation's start.

    The result is sorted, which is what makes the table's row order stable.
    """
    lineStop = lineLength
    if lineLength <= 2:
        return [[]]
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


def relSetToString(relSet):
    """A relation set -> '1;2;3|3;4;5;6', the key every table row is found by.

    The inverse of `relationStringToLineRelLengths` only up to the relation
    lengths: this writes out the vertices a relation passes through, so it works
    on the relations of any quiver, not only a line.  Each relation contributes
    its first path, which for an LNA is its only one.
    """
    return "|".join(";".join(str(v) for v in rel[0]) for rel in relSet)


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


def className(relLengths):
    """[2, 0, 3, 0] -> '2030', the name an LNA's class goes by.

    The third of the three forms, and the one the class names in the tables are
    written in.  `LinearNakayamaAlgebra.className` is this on an algebra.
    """
    return ''.join(str(n) for n in relLengths)


def relabelLineAlgebra(pathAlg, currentRelabeling = None):
    """Renumber a quiver that is a line so that its arrows run 1 -> 2 -> ... -> n.

    A mutation keeps every vertex label, so a quiver that comes back as a line
    is a line with the labels in some other order; renumbering it is what lets
    it be recognised as an LNA and looked up in the table.  Returns the
    renumbered algebra and the map from new labels back to whatever the caller
    was already tracking, so a chain of mutations can be read in the original
    numbering.

    **Renumbers in place** and returns the same object, which is why
    `lnaMoves` copies before calling it.  A quiver that is not a line is
    returned untouched -- as its `quiver`, not as the algebra, which is a wart
    the caller has to know about.
    """
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


def mutationListLineCleanup(mutationList, relabelNodes = True, printOutput = True):
    """Reduce what a search collected to one entry per LNA, shortest path first.

    A depth-first search reaches the same LNA many times, by paths of different
    lengths and in different numberings.  This keeps the quivers that are lines,
    renumbers them to the standard 1 -> ... -> n, and then keeps one entry per
    relation set: the one with the fewest mutations.  The result is sorted by
    relation set, so two searches of the same class produce the same list.
    """
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
                newVertexRelabeling = vertexRelabeling
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
