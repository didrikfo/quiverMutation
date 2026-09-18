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
        # The procedure works on linear combinations of paths, `rels` records
        # only the sets.  When the coefficients are known they are kept here,
        # in the same order as `rels`, so a walk of several mutations does not
        # have to guess them back at every step.  It is a cache, not a second
        # source of truth: `procedure.relationsFrom` checks it still describes
        # `rels` and falls back to the guess if anything has edited `rels`
        # behind its back.
        self.relCombinations = None
        # The arrows of a mutated quiver are not named by their endpoints: the
        # procedure produces parallel arrows, and then a path as a sequence of
        # vertices no longer says which arrow it runs along.  So the relations
        # are also kept as `arrowPaths` combinations, over paths that are tuples
        # of `(tail, head, key)` arrows, and *that* is the faithful record --
        # `rels` is its projection and the table's key.  `procedure.relationsFrom`
        # prefers this and checks it still describes `rels` and still names
        # arrows of the quiver before it does.  See NOTES.md, "Parallel arrows".
        self.arrowRels = None

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
        return self.quiver.add_edge(arrowStart, arrowEnd)

    def hasParallelArrows(self):
        """Whether two distinct arrows share both endpoints.

        Where this is true, `rels` does not determine the algebra and
        `arrowRels` is the only faithful reading of it.
        """
        from . import arrowPaths
        return arrowPaths.hasParallelArrows(self.quiver)

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

    # -- the operations, as methods ---------------------------------------
    #
    # The procedure, the reduction and the invariants are free functions in
    # the modules that own them, and this is their object-oriented face: the
    # methods delegate, they do not reimplement.  New code should read as
    # `algebra.mutateAt(3).reduce()`; the free functions remain the
    # implementation, and what the internals call.
    #
    # The imports are inside the methods on purpose.  `pathAlgebra` is the root
    # of the package's dependency graph -- `mutation`, `reduction` and
    # `invariants` all import it -- so importing them here at module level
    # would make that a cycle.  A deferred import is the price of the container
    # being at the bottom and still having behaviour.

    def canMutateAt(self, vertex):
        """Whether the mutation procedure may be applied at `vertex`.

        Applying it anyway still returns a quiver, just not a derived
        equivalent one -- research R-005 -- so this is not an optional check.
        """
        from . import mutation
        return mutation.mutationIsPossibleAtVertex(self, vertex)

    def mutateAt(self, vertex):
        """Mutate at one vertex, without reducing.  Negative means left."""
        from . import mutation
        if vertex < 0:
            return mutation.leftQuiverMutationAtVertex(self, -vertex)
        return mutation.quiverMutationAtVertex(self, vertex)

    def mutateAtVertices(self, vertices, printMutationSteps = False):
        """Mutate at each vertex in turn, reducing after each."""
        from . import mutation
        return mutation.quiverMutationAtVertices(self, vertices, printMutationSteps)

    def reduce(self):
        """The cleanup the procedure's steps 1-7 leave to the caller."""
        from . import reduction
        return reduction.reducePathAlgebra(self)

    def opposite(self):
        """The opposite algebra: every arrow and every relation reversed."""
        return dualPathAlgebra(self)

    def cartanMatrix(self, exact = True):
        from . import invariants
        return invariants.cartanMatrix(self, exact)

    def coxeterPolynomial(self, exact = True):
        from . import invariants
        return invariants.coxeterPoly(self, exact).as_expr()


def printPathAlgebra(pathAlg):
    print('Vertices: ', pathAlg.quiver.nodes)
    print('Arrows: ', pathAlg.quiver.edges)
    print('Relations: ', pathAlg.rels, '\n')
    return


def dualPathAlgebra( pathAlg ):
    """The opposite algebra: every arrow and every relation reversed.

    Arrow names are carried across as well, `(tail, head, key)` going to
    `(head, tail, key)`, which is injective -- so a parallel pair stays a
    parallel pair and the arrow relations reverse with it.  Without that, the
    dual of a quiver with parallel arrows could not be read back at all, and
    left mutation goes through the dual.
    """
    dualPathAlg = PathAlgebra()
    dualPathAlg.add_vertices_from(pathAlg.vertices())
    for tail, head, key in pathAlg.quiver.edges(keys = True):
        dualPathAlg.quiver.add_edge(head, tail, key = key)
    for rel in pathAlg.rels:
        dualRel = []
        for relPath in rel:
            dualRel.append(list(reversed(relPath)))
        dualPathAlg.add_rel(dualRel)
    arrowRels = getattr(pathAlg, "arrowRels", None)
    if arrowRels is not None:
        from . import arrowPaths
        dualArrowRels = [
            arrowPaths.combination(
                (tuple((head, tail, key) for tail, head, key in reversed(path)), coefficient)
                for path, coefficient in relation.items())
            for relation in arrowRels
        ]
        # `add_rel` sorts each relation's paths, and reversing a relation set
        # reverses the order they were written in, so the two lists have to be
        # matched up again rather than zipped as they stand.
        byProjection = {}
        for relation in dualArrowRels:
            key = tuple(sorted(tuple(arrowPaths.projectPath(p)) for p in relation))
            byProjection.setdefault(key, []).append(relation)
        matched = []
        for rel in dualPathAlg.rels:
            key = tuple(sorted(tuple(p) for p in rel))
            candidates = byProjection.get(key)
            if not candidates:
                matched = None
                break
            matched.append(candidates.pop(0))
        dualPathAlg.arrowRels = matched
    return dualPathAlg
