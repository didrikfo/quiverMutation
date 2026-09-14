"""The container: a quiver and a list of relations, and nothing more.

A relation is a list of paths and a path is a list of vertices, so
`[[1,2,4],[1,3,4]]` is the commutativity relation between the two paths from 1
to 4 and `[[1,2,3]]` is the zero relation on `1 -> 2 -> 3`.  By convention
`rel[0][0]` is a relation's source and `rel[0][-1]` its target, and every path
in a relation shares both.  See NOTES.md, "The model of a path algebra", for
what this deliberately cannot express.

Every operation on a path algebra is a free function elsewhere in the package
taking one as its first argument; `nakayama` holds the two subclasses that do
carry their own behaviour.
"""

import networkx as nx


class PathAlgebra():

    def __init__(self):
        self.quiver = nx.MultiDiGraph()
        self.rels = []

    def vertices(self):
        """ returns the vertices of a quiver """
        return list(self.quiver.nodes)

    def arrows(self):
        """ returns the arrows of a quiver """
        return self.quiver.edges

    def add_vertex(self, vertex):
        if vertex not in self.quiver.nodes:
            self.quiver.add_node(vertex)

    def add_vertices_from(self, vertices):
        for vertex in vertices:
            if vertex not in self.quiver:
                self.quiver.add_node(vertex)

    def add_arrow(self, arrowStart, arrowEnd):
        self.quiver.add_edge(arrowStart, arrowEnd)

    def add_arrows_from(self, arrows):
        self.quiver.add_edges_from(arrows)

    def add_path(self, path):
        for i in range(len(path) - 1):
            self.quiver.add_edge(path[i], path[i + 1])

    def add_paths_from(self, pathList):
        for path in pathList:
            for i in range(len(path) - 1):
                self.quiver.add_edge(path[i], path[i + 1])

    def add_rel(self, rel):
        self.rels.append(sorted(rel)[:])

    def add_rels_from(self, rels):
        for rel in rels:
            rel.sort()
        self.rels.extend(rels[:])

    def update_quiver(self, newQuiver, newRels = None):
        newRels = [] if newRels is None else newRels
        self.quiver = newQuiver
        self.rels = newRels

    def out_arrows(self, vertex):
        return self.quiver.out_edges(vertex)

    def out_rels(self, vertex):
        outRels = []
        for rel in self.rels:
            if rel[0][0] == vertex:
                outRels.append(rel)
        return outRels

    def in_rels(self, vertex):
        inRels = []
        for rel in self.rels:
            if rel[0][-1] == vertex:
                inRels.append(rel)
        return inRels

    def rels_between(self, inVertex, outVertex):
        betweenRels = []
        for rel in self.rels:
            if rel[0][0] == inVertex and rel[0][-1] == outVertex:
                betweenRels.append(rel)
        return betweenRels

    def clear_rels(self):
        self.rels = []
        return

def printPathAlgebra(pathAlg):
    print('Vertices: ', pathAlg.quiver.nodes)
    print('Arrows: ', pathAlg.quiver.edges)
    print('Relations: ', pathAlg.rels, '\n')
    return


def dualPathAlgebra( pathAlg ):
    dualPathAlg = PathAlgebra()
    dualPathAlg.add_vertices_from(pathAlg.vertices())
    for arrow in pathAlg.arrows():
        dualPathAlg.add_arrow(arrow[1], arrow[0])
    for rel in pathAlg.rels:
        dualRel = []
        for relPath in rel:
            dualRel.append(list(reversed(relPath)))
        dualPathAlg.add_rel(dualRel)
    return dualPathAlg
