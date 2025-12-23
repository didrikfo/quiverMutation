import networkx as nx


class PathAlgebra:
    def __init__(self):
        self.quiver = nx.MultiDiGraph()
        self.rels = []

    def vertices(self):
        """Returns the vertices of a quiver."""
        return list(self.quiver.nodes)

    def arrows(self):
        """Returns the arrows of a quiver."""
        return self.quiver.edges

    def add_vertex(self, vertex):
        if vertex not in self.quiver.nodes:
            self.quiver.add_node(vertex)

    def add_vertices_from(self, vertices):
        for vertex in vertices:
            if vertex not in self.quiver:
                self.quiver.add_node(vertex)

    def add_arrow(self, arrow_start, arrow_end):
        self.quiver.add_edge(arrow_start, arrow_end)

    def add_arrows_from(self, arrows):
        self.quiver.add_edges_from(arrows)

    def add_path(self, path):
        for i in range(len(path) - 1):
            self.quiver.add_edge(path[i], path[i + 1])

    def add_paths_from(self, path_list):
        for path in path_list:
            for i in range(len(path) - 1):
                self.quiver.add_edge(path[i], path[i + 1])

    def add_rel(self, rel):
        self.rels.append(sorted(rel)[:])

    def add_rels_from(self, rels):
        for rel in rels:
            rel.sort()
        self.rels.extend(rels[:])

    def update_quiver(self, new_quiver, new_rels=None):
        if new_rels is None:
            new_rels = []
        self.quiver = new_quiver
        self.rels = new_rels

    def out_arrows(self, vertex):
        return self.quiver.out_edges(vertex)

    def out_rels(self, vertex):
        out_rels = []
        for rel in self.rels:
            if rel[0][0] == vertex:
                out_rels.append(rel)
        return out_rels

    def in_rels(self, vertex):
        in_rels = []
        for rel in self.rels:
            if rel[0][-1] == vertex:
                in_rels.append(rel)
        return in_rels

    def rels_between(self, in_vertex, out_vertex):
        between_rels = []
        for rel in self.rels:
            if rel[0][0] == in_vertex and rel[0][-1] == out_vertex:
                between_rels.append(rel)
        return between_rels

    def clear_rels(self):
        self.rels = []
        return
