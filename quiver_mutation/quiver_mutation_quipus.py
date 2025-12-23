import csv
import copy
import itertools
import math
import networkx as nx
import numpy as np
import sympy
from sympy.matrices import Matrix

def generate_all_quipus_up_to_length( length ):
    timeStart = time.time()
    quipusOfAllLengths = [[nx.path_graph(1)]]
    pathGraphs = []
    for i in range(length):
        pathGraphs.append(nx.path_graph(i+1))
        quipusOfThisLength = []
        for j in range(int(np.floor((i+3)/2)), i+1 ):
            pathGraph = pathGraphs[j]
            for vertexSet in itertools.combinations(range(1,j), i + 1 - j):
                heightOneQuipu = pathGraph.copy()
                k = 1
                for v in vertexSet:
                    heightOneQuipu.add_edge(v, j+k)
                    k += 1
                isNewQuipu = True
                for Q in quipusOfThisLength:
                    if nx.is_isomorphic(heightOneQuipu, Q):
                        isNewQuipu = False
                        break
                if isNewQuipu:
                    quipusOfThisLength.append(heightOneQuipu)
        for quipu in quipusOfAllLengths[i]:
            for v in quipu.nodes():
                if quipu.degree(v) <= 1:
                    longerQuipu = quipu.copy()
                    longerQuipu.add_edge(v, i+1)
                    isNewQuipu = True
                    for Q in quipusOfThisLength:
                        if nx.is_isomorphic(longerQuipu, Q):
                            isNewQuipu = False
                            break
                    if isNewQuipu:
                        quipusOfThisLength.append(longerQuipu)
        quipusOfAllLengths.append([])
        for newQuipu in quipusOfThisLength:
            quipusOfAllLengths[i+1].append(newQuipu.copy())
        timeEnd = time.time()
        print('Generated ', len(quipusOfThisLength), ' quipus of length ', i+2, ' in ',  timeEnd - timeStart, 's.')
    return quipusOfAllLengths

def generate_all_height_one_quipus( mainStringLength, mainStringLengthStart = 1 ):
    startTime = time.time()
    heightOneQuipusWithDupes = [nx.path_graph(1)]
    for i in range(mainStringLengthStart,mainStringLength):
        heightOneQuipusWithDupes.append(nx.path_graph(i+1))
        vertexPowerset = powerset(range(1,i))
        j = 0
        while j < 2**(i-1):
            vertexSet = vertexPowerset[j]
            setToRemove = {}
            for v in vertexSet:
                setToRemove.append(i-v)
            if setToRemove:
                print('before', vertexPowerset)
                vertexPowerset = vertexPowerset - setToRemove
                print('after', vertexPowerset)
                j += 1
            j += 1
        for vertexSet in vertexPowerset:
            heightOneQuipu = nx.path_graph(i+1)
            if bool(vertexSet):
                k = 1
                for v in vertexSet:
                    heightOneQuipu.add_edge(v, i+k)
                    k += 1
                heightOneQuipusWithDupes.append(heightOneQuipu)
        print('Done with ', i+1, ' in ', time.time() - startTime, 's')
    heightOneQuipus = []
    print('Found ', len(heightOneQuipusWithDupes), ' height one quipus with dupes.')
    for Q1 in heightOneQuipusWithDupes:
        isNewQuipu = True
        for Q2 in heightOneQuipus:
            if nx.is_isomorphic(Q1, Q2):
                isNewQuipu = False
                break
        if isNewQuipu:
            heightOneQuipus.append(Q1)
    return heightOneQuipus

def save_quipus_to_csv( quipusList, fileName, overwrirteFile = True ):
    if overwrirteFile:
        open(fileName, 'w+').close()
    with open(fileName, 'a') as f:
        for Q in quipusList:
            f.write('{0}\n'.format(str(Q.edges)))
        f.close()
    return

def generate_all_quipus( length ):
    trees = list(nx.generators.nonisomorphic_trees(length))
    quipus = []
    for tree in trees:
        isQuipu = False
        maxDeg = 0
        for deg in tree.degree:
            if deg[1] > maxDeg:
                maxDeg = deg[1]
                if deg[1] > 3:
                    break
        if maxDeg <= 3:
            isQuipu = True
            deg3SubGraph = nx.Graph()
            deg3Vertices = []
            for v in tree.nodes:
                if tree.degree(v) == 3:
                    deg3Vertices.append(v)
            if len(deg3Vertices) > 3:
                nx.add_path(deg3SubGraph, nx.shortest_path(tree, deg3Vertices[0], deg3Vertices[1]))
                for v1 in deg3Vertices[2:]:
                    pathToAddAsGraph = nx.Graph()
                    nx.add_path(pathToAddAsGraph, nx.shortest_path(tree, deg3Vertices[0], v1))
                    deg3SubGraph = nx.compose(deg3SubGraph, pathToAddAsGraph)
                maxDegForDeg3SubGraph = max(deg3SubGraph.degree,key=itemgetter(1))[1]
                if maxDegForDeg3SubGraph >= 3:
                    isQuipu = False
        if isQuipu:
            quipus.append(tree)
    return quipus

def generate_all_quipus_gpt(length):
    quipus = []
    for tree in nx.generators.nonisomorphic_trees(length):
        maxDeg = max(dict(tree.degree()).values())
        if maxDeg <= 3:
            isQuipu = True
            deg3Vertices = [v for v, d in tree.degree() if d == 3]
            if len(deg3Vertices) > 3:
                deg3SubGraph = nx.Graph()
                nx.add_path(deg3SubGraph, nx.shortest_path(tree, deg3Vertices[0], deg3Vertices[1]))
                for v1 in deg3Vertices[2:]:
                    pathToAddAsGraph = nx.Graph()
                    nx.add_path(pathToAddAsGraph, nx.shortest_path(tree, deg3Vertices[0], v1))
                    deg3SubGraph = nx.compose(deg3SubGraph, pathToAddAsGraph)
                maxDegForDeg3SubGraph = max(dict(deg3SubGraph.degree()).values())
                if maxDegForDeg3SubGraph >= 3:
                    isQuipu = False
        else:
            isQuipu = False
        if isQuipu:
            quipus.append(tree)
    return quipus

def generate_quipus(length):
    for tree in nx.generators.nonisomorphic_trees(length):
        is_quipu = False
        max_deg = 0
        for deg in tree.degree:
            if deg[1] > max_deg:
                max_deg = deg[1]
                if deg[1] > 3:
                    break
        if max_deg <= 3:
            is_quipu = True
            deg3_subgraph = nx.Graph()
            deg3_vertices = []
            for v in tree.nodes:
                if tree.degree(v) == 3:
                    deg3_vertices.append(v)
            if len(deg3_vertices) > 3:
                nx.add_path(deg3_subgraph, nx.shortest_path(tree, deg3_vertices[0], deg3_vertices[1]))
                for v1 in deg3_vertices[2:]:
                    path_to_add_as_graph = nx.Graph()
                    nx.add_path(path_to_add_as_graph, nx.shortest_path(tree, deg3_vertices[0], v1))
                    deg3_subgraph = nx.compose(deg3_subgraph, path_to_add_as_graph)
                max_deg_for_deg3_subgraph = max(deg3_subgraph.degree, key=itemgetter(1))[1]
                if max_deg_for_deg3_subgraph >= 3:
                    is_quipu = False
        if is_quipu:
            yield tree

def count_quipus_v1(length):
    num_quipus = 0
    for tree in nx.generators.nonisomorphic_trees(length):
        is_quipu = False
        max_deg = 0
        for deg in tree.degree:
            if deg[1] > max_deg:
                max_deg = deg[1]
                if deg[1] > 3:
                    break
        if max_deg <= 3:
            is_quipu = True
            deg3_vertices = []
            for v in tree.nodes:
                if tree.degree(v) == 3:
                    deg3_vertices.append(v)
            if len(deg3_vertices) > 3:
                #deg3_subgraph = nx.Graph()
                #nx.add_path(deg3_subgraph, dfs_shortest_path(tree, deg3_vertices[0], deg3_vertices[1]))
                deg3_path = dfs_shortest_path(tree, deg3_vertices[0], deg3_vertices[-1])
                #deg3_subgraph_edges = [(deg3_path[x], deg3_path[x+1]) for x in range(len(deg3_path)-1)]
                ##for v1 in deg3_vertices[2:]:
                ##    path_to_add = nx.Graph()
                ##    #nx.add_path(path_to_add, dfs_shortest_path(tree, deg3_vertices[0], v1))
                ##    nx.add_path(path_to_add, bfs_shortest_path_to_subgraph(tree, v1, deg3_subgraph))
                ##    deg3_subgraph = nx.compose(deg3_subgraph, path_to_add)

                deg3_vertices_not_checked = [x for x in deg3_vertices if x not in deg3_path]
                for i in range(1, len(deg3_vertices)-1):
                    v1 = deg3_vertices[i]
                    if v1 in deg3_vertices_not_checked:
                        ##path_as_list = bfs_shortest_path_to_subgraph(tree, v1, deg3_subgraph)
                        #path_as_list = bfs_shortest_path_to_subgraph_edges(tree, v1, deg3_subgraph_edges)
                        path_as_list = bfs_shortest_path_to_subgraph_path(tree, v1, deg3_path)
                        if len(path_as_list) > 1:
                            if path_as_list[-1] == deg3_path[0]:
                                #print(path_as_list, deg3_path)
                                deg3_path = path_as_list[:-1] + deg3_path
                                #print(path_as_list, deg3_path)
                            elif path_as_list[-1] == deg3_path[-1]:
                                deg3_path.extend(reversed(path_as_list[:-1]))
                            else:
                                #print(deg3_path)
                                #print(path_as_list)
                                is_quipu = False
                                break
                            #connection_point_degree = sum(edge.count(path_as_list[-1]) for edge in deg3_subgraph_edges)
                            #if connection_point_degree == 2:
                            #    is_quipu = False
                            #    break
                            ##path_to_add = nx.Graph()
                            ##nx.add_path(path_to_add, path_as_list)
                            ##deg3_subgraph = nx.compose(deg3_subgraph, path_to_add)
                            #deg3_subgraph_edges.extend([(path_as_list[x], path_as_list[x+1]) for x in range(len(path_as_list)-1)])
                            path_as_set = set(path_as_list)
                            deg3_vertices_not_checked = [x for x in deg3_vertices_not_checked if not x in path_as_set]

                ##max_deg_for_deg3_subgraph = max(deg3_subgraph.degree, key=itemgetter(1))[1]
                ##if max_deg_for_deg3_subgraph >= 3:
                ##    is_quipu = False
        if is_quipu:
            num_quipus += 1
    return num_quipus

def dfs_shortest_path(tree, start, end):
    visited = set()
    stack = [(start, [start])]
    while stack:
        (node, path) = stack.pop()
        if node == end:
            return path
        if node not in visited:
            visited.add(node)
            for neighbor in tree.neighbors(node):
                if neighbor not in visited:
                    stack.append((neighbor, path + [neighbor]))
    return None

def bfs_shortest_path_to_subgraph(tree, start_node, subgraph):
    queue = [(start_node, [start_node])]
    visited = set()

    while queue:
        node, path = queue.pop(0)

        if node in subgraph:
            return path

        visited.add(node)

        for neighbor in tree.neighbors(node):
            if neighbor not in visited:
                queue.append((neighbor, path + [neighbor]))

    return []

def bfs_shortest_path_to_subgraph_edges(tree, start_node, subgraph_edges):
    queue = [(start_node, [start_node])]
    visited = set()

    while queue:
        node, path = queue.pop(0)

        if any(node in edge for edge in subgraph_edges):
            return path

        visited.add(node)

        for neighbor in tree.neighbors(node):
            if neighbor not in visited:
                edge = (node, neighbor)
                if edge not in subgraph_edges and (edge[1], edge[0]) not in subgraph_edges:
                    queue.append((neighbor, path + [neighbor]))

    return None

def bfs_shortest_path_to_subgraph_path(tree, start_node, subgraph_path):
    queue = [(start_node, [start_node])]
    visited = set()
    while queue:
        node, path = queue.pop(0)

        if node in subgraph_path:
            return path

        visited.add(node)
        for neighbor in tree.neighbors(node):
            #print(neighbor)
            if neighbor not in visited:
                queue.append((neighbor, path + [neighbor]))
    return None

def count_quipus(n):
    count = [0] * (n + 1)  # Initialize count list with zeros

    # Special case for n=2
    count[2] = 1

    # Iterate over all possible lengths of quipus
    for length in range(3, n + 1):
        # Generate all non-isomorphic trees with the given length
        trees = []
        if length == 3:
            G = nx.path_graph(length)
            trees.append(G)
        else:
            for tree in nx.nonisomorphic_trees(length):
                trees.append(tree)

        # Count the number of quipus with the given length
        for tree in trees:
            degree_3_nodes = [node for node, degree in tree.degree() if degree == 3]
            if len(degree_3_nodes) < 4:
                continue
            for path in nx.all_simple_paths(tree, degree_3_nodes[0], degree_3_nodes[-1]):
                if all(tree.degree(node) <= 3 for node in path):
                    count[length] += 1
                    break

    return count
