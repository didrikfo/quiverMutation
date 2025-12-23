import ast
import copy
import os
import matplotlib.pyplot as plt
import networkx as nx
from . import path_algebra_class

def print_path_algebra(pathAlg):
    print('Vertices: ', pathAlg.quiver.nodes)
    print('Arrows: ', pathAlg.quiver.edges)
    print('Relations: ', pathAlg.rels, '\n')
    return

def plot_quiver(pathAlg, showPlot = True, saveToFile = False, fileName = 'quiverPlot.png', folder = ''):
    try:
        os.mkdir(folder)
    except OSError:
        print("Creation of the directory %s failed" % folder)
    else:
        print("Successfully created the directory %s " % folder)
    savePath = folder + fileName
    nx.draw_networkx(pathAlg.quiver)
    if saveToFile:
        plt.savefig(savePath)
        plt.close()
    if showPlot:
        plt.show()
    return

def read_mutations_from_file(fileName):
    rPathAlg = path_algebra_class.PathAlgebra()
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

def read_relations_from_file(fileName):
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
