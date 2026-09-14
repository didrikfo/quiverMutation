"""Path algebras that know what kind of algebra they are.

`pathAlgebraClass.PathAlgebra` is a quiver and a list of relations and nothing
more, and every operation on it is a free function in `quiverMutation` taking it
as the first argument.  The two shapes this project actually works with -- the
linearly oriented Nakayama algebras being classified, and the quipu algebras
that turn out to classify them -- carry much more structure than that, and
carrying it explicitly removes a lot of the conversion between per-vertex
relation lengths, relation strings, class names and Kupisch series that is
currently spread across the module.
"""

import networkx as nx

import pathAlgebraClass
import quipuForms
import quiverMutation as qm
import relationAlgebra


class LinearNakayamaAlgebra(pathAlgebraClass.PathAlgebra):
    """kA_n / I: the linear quiver 1 -> 2 -> ... -> n with an admissible ideal.

    `relLengths[i]` is the number of arrows in the relation starting at vertex
    i + 1, or 0 if no relation starts there, so it has n - 2 entries.  This is
    the same per-vertex form the rest of the repo uses, and the digits of the
    class names: A7_22300 is the line on 7 vertices with relations of 2, 2 and 3
    arrows starting at vertices 1, 2 and 3.

    There are Catalan(n - 1) of these for each n.
    """

    def __init__(self, length, relLengths):
        super().__init__()
        if isinstance(relLengths, str):
            relLengths = [int(c) for c in relLengths]
        relLengths = list(relLengths)
        if len(relLengths) != length - 2:
            raise ValueError(
                "a line of length {0} takes {1} relation lengths, got {2}".format(
                    length, length - 2, len(relLengths)))
        for start, arrows in enumerate(relLengths, start=1):
            if arrows and (arrows < 2 or start + arrows > length):
                raise ValueError(
                    "the relation at vertex {0} has {1} arrows, which does not fit "
                    "in a line of length {2}".format(start, arrows, length))
        self.length = length
        self.relLengths = relLengths
        self.add_vertices_from(range(1, length + 1))
        self.add_arrows_from([[i, i + 1] for i in range(1, length)])
        self.add_rels_from([
            [list(range(start, start + arrows + 1))]
            for start, arrows in enumerate(relLengths, start=1)
            if arrows
        ])

    # -- naming ------------------------------------------------------------

    @classmethod
    def fromRelationString(cls, length, relationString):
        """From the 'a;b;c|d;e;f' form used as the table's key."""
        return cls(length, qm.relationStringToLineRelLengths(length, relationString))

    @classmethod
    def fromClassName(cls, className):
        """From '22300', which also fixes the length: n = len(className) + 2."""
        return cls(len(className) + 2, className)

    def relationString(self):
        return "|".join(
            ";".join(str(v) for v in range(start, start + arrows + 1))
            for start, arrows in enumerate(self.relLengths, start=1)
            if arrows
        )

    def className(self):
        return "".join(str(n) for n in self.relLengths)

    def __repr__(self):
        return "LinearNakayamaAlgebra({0}, {1!r})".format(self.length, self.className())

    def __eq__(self, other):
        return (isinstance(other, LinearNakayamaAlgebra)
                and self.length == other.length
                and self.relLengths == other.relLengths)

    def __hash__(self):
        return hash((self.length, tuple(self.relLengths)))

    # -- structure ---------------------------------------------------------

    def relations(self):
        """(start vertex, number of arrows) for each relation, left to right."""
        return [(start, arrows) for start, arrows in enumerate(self.relLengths, start=1) if arrows]

    def kupischSeries(self):
        """(dim P_1, ..., dim P_n), the classical invariant of a Nakayama algebra.

        dim P_i is the number of vertices on the longest nonzero path out of i.
        A path from i onwards dies at the first relation that starts at or after
        i, so dim P_i = min(n, (that relation's end) - 1) - i + 1.
        """
        series = []
        for i in range(1, self.length + 1):
            firstEnd = min(
                (start + arrows for start, arrows in self.relations() if start >= i),
                default=self.length + 1,
            )
            series.append(min(self.length, firstEnd - 1) - i + 1)
        return tuple(series)

    def hasAlmostSeparateRelations(self):
        """Whether consecutive relations overlap in at most one arrow.

        The condition n_{i+1} >= n_i + l_i - 1 of arXiv:2305.06642, which is what
        the published classification covers.
        """
        return quipuForms._hasAlmostSeparateRelations(self.length, self.relations())

    def relationDual(self):
        """The opposite algebra, renumbered back to 1 -> ... -> n.

        Reversing every arrow is one of the three class-preserving operations of
        arXiv:2305.06642, so this is always derived equivalent to self.
        """
        dual = [0] * (self.length - 2)
        for start, arrows in self.relations():
            dual[self.length - start - arrows] = arrows
        return LinearNakayamaAlgebra(self.length, dual)

    # -- invariants --------------------------------------------------------

    def quipu(self):
        """The quipu this algebra is derived equivalent to, or None.

        None means the algebra does not have almost separate relations, in which
        case theorem `thm:QuipuToAn` of arXiv:2305.06642 says nothing about it
        and the class has to be found by mutation instead.
        """
        return quipuForms.quipuForAlmostSeparateLNA(self.length, self.relLengths)

    def quipuName(self):
        return quipuForms.formatQuipu(self.quipu())

    def cartanMatrix(self):
        return relationAlgebra.cartanMatrixExact(self)

    def coxeterPolynomial(self):
        return qm.coxeterPoly(self).as_expr()

    # -- enumeration -------------------------------------------------------

    @classmethod
    def allOfLength(cls, length):
        """Every LNA of the given length, in the table's order."""
        return [
            cls.fromRelationString(length, qm.relSetToString(relSet))
            for relSet in qm.generateAllPossibleLineRelations(length)
        ]


class QuipuAlgebra(pathAlgebraClass.PathAlgebra):
    """The path algebra of a quipu quiver, which has no relations.

    A quipu is a tree of maximum degree 3 whose degree-3 vertices all lie on one
    path; the quipu quiver orients the main string left to right and each cord
    away from the main string.  Since the algebra is hereditary, every
    orientation of the same quipu gives a derived equivalent algebra, so the
    underlying undirected tree is the whole of its derived equivalence class.
    """

    def __init__(self, k, m):
        super().__init__()
        if len(k) != len(m) + 1:
            raise ValueError("a quipu needs one more k than m, got {0} and {1}".format(k, m))
        self.k = tuple(k)
        self.m = tuple(m)
        graph = quipuForms.graphFromQuipuParameters(self.k, self.m)
        self.mainString, self.cords = self._orient(graph)
        self.add_vertices_from(sorted(graph.nodes))
        for path in [self.mainString] + self.cords:
            for earlier, later in zip(path, path[1:]):
                self.add_arrow(earlier, later)

    def _orient(self, graph):
        """Recover the main string and cords, oriented as the paper's D quiver."""
        degrees = dict(graph.degree())
        branchVertices = {v for v, d in degrees.items() if d == 3}
        leaves = sorted(v for v, d in degrees.items() if d <= 1)
        for first in leaves:
            for last in leaves:
                if first == last:
                    continue
                candidate = nx.shortest_path(graph, first, last)
                if not branchVertices.issubset(candidate):
                    continue
                onMain = set(candidate)
                cords = []
                ok = True
                for vertex in candidate:
                    if vertex not in branchVertices:
                        continue
                    off = [w for w in graph.neighbors(vertex) if w not in onMain]
                    if len(off) != 1:
                        ok = False
                        break
                    cord = [vertex, off[0]]
                    while True:
                        onward = [w for w in graph.neighbors(cord[-1]) if w != cord[-2]]
                        if not onward:
                            break
                        cord.append(onward[0])
                    cords.append(cord)
                if ok and len(cords) == len(branchVertices):
                    return candidate, cords
        return list(graph.nodes), []

    @classmethod
    def fromLNA(cls, lna):
        """The quipu algebra an LNA with almost separate relations comes from."""
        parameters = lna.quipu()
        if parameters is None:
            raise ValueError(
                "{0!r} does not have almost separate relations, so the quipu "
                "theorem does not apply to it".format(lna))
        return cls(*parameters)

    def quipuParameters(self):
        return (self.k, self.m)

    def quipuName(self):
        return quipuForms.formatQuipu((self.k, self.m))

    def underlyingGraph(self):
        return quipuForms.underlyingGraph(self)

    def canonicalForm(self):
        return quipuForms.canonicalTreeForm(self.underlyingGraph())

    def correspondingLNA(self):
        """The LNA this quipu is derived equivalent to, by `thm:QuipuToAn`.

        n_i = k_0 + sum_{j=1..i} (m_{j-1} + k_j + 1); the algebra is
        A_{n_{r+1}, (k_0, n_1, ..., n_r)}^{(m_0 + 2, ..., m_r + 2)}.
        """
        starts = [self.k[0]]
        for i in range(1, len(self.k)):
            starts.append(starts[-1] + self.m[i - 1] + self.k[i] + 1)
        length = starts[-1]
        relLengths = [0] * (length - 2)
        for index, cord in enumerate(self.m):
            relLengths[starts[index] - 1] = cord + 2
        return LinearNakayamaAlgebra(length, relLengths)

    def coxeterPolynomial(self):
        return qm.coxeterPoly(self).as_expr()

    def __repr__(self):
        return "QuipuAlgebra({0}, {1})".format(self.k, self.m)

    def __eq__(self, other):
        return isinstance(other, QuipuAlgebra) and (self.k, self.m) == (other.k, other.m)

    def __hash__(self):
        return hash((self.k, self.m))
