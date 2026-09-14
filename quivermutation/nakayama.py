"""Path algebras that know what kind of algebra they are.

`pathAlgebra.PathAlgebra` is a quiver and a list of relations and nothing more,
and every operation on it is a free function elsewhere in the package taking it
as the first argument.  The two shapes this project actually works with -- the
linearly oriented Nakayama algebras being classified, and the quipu algebras
that turn out to classify them -- carry much more structure than that, and
carrying it explicitly removes a lot of the conversion between per-vertex
relation lengths, relation strings, class names and Kupisch series that `lines`
still does by hand.
"""

import networkx as nx

from . import lines
from . import pathAlgebra
from . import quipuForms



class LinearNakayamaAlgebra(pathAlgebra.PathAlgebra):
    """kA_n / I: the linear quiver 1 -> 2 -> ... -> n with an admissible ideal.

    `relLengths[i]` is the number of arrows in the relation starting at vertex
    i + 1, or 0 if no relation starts there, so it has max(0, n - 2) entries --
    A_1 and A_2 admit no relation and so take none.  This is
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
        if length < 1:
            raise ValueError("a line needs at least one vertex, got {0}".format(length))
        # A_1 and A_2 admit no relation at all, so they take no relation lengths
        # -- length - 2 would be negative for A_1.
        expected = max(0, length - 2)
        if len(relLengths) != expected:
            raise ValueError(
                "a line of length {0} takes {1} relation lengths, got {2}".format(
                    length, expected, len(relLengths)))
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
        return cls(length, lines.relationStringToLineRelLengths(length, relationString))

    @classmethod
    def fromClassName(cls, className):
        """From '22300', which also fixes the length: n = len(className) + 2."""
        return cls(len(className) + 2, className)

    def relationString(self):
        return lines.relSetToString(self.rels)

    def className(self):
        return lines.className(self.relLengths)

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

    # -- the class-preserving operations of the paper ----------------------
    #
    # `cor:EquivNakayamaAlgebras` of arXiv:2305.06642 lists three operations
    # that move an LNA with almost separate relations to another one in the same
    # derived equivalence class, and says there are at most eight such algebras
    # per class once every relation has length >= 3.  They are here as methods
    # because the class of an LNA was otherwise only ever computed through the
    # quipu, where the same symmetry appears indirectly: two parameter pairs that
    # name the same tree.  Having both routes means each can check the other, and
    # tests/test_quipu_symmetry.py does exactly that.

    def withoutShortRelations(self):
        """The same class, with every relation of length 2 dropped.

        Operation 2: a relation of two arrows never changes the derived
        equivalence class of an LNA with almost separate relations, which is why
        `quipu` drops them before naming the class.  The other two operations are
        stated for algebras whose relations all have length >= 3, so they apply
        to this form of the algebra.
        """
        return LinearNakayamaAlgebra(
            self.length, [0 if arrows == 2 else arrows for arrows in self.relLengths])

    def swapFirstRelation(self):
        """Operation 1 at the first relation, or None where it does not apply.

        In the quipu, the first relation's start vertex is the length k_0 of the
        main string before the first cord foot, and the relation's own length is
        m_0 + 2, where m_0 is that cord.  The two branches at the foot are the
        k_0 path and the m_0 cord, so exchanging them is an isomorphism of the
        tree -- which on the LNA side replaces the first relation
        (n_0, l_0) by (l_0 - 2, n_0 + 2):

            the relation's END vertex n_0 + l_0 does not move; its start vertex
            and its length trade places.

        None when there is no relation of length >= 3 to swap, or when the swap
        would need a start vertex before the first (l_0 = 2 after dropping the
        short relations means there is nothing to exchange).
        """
        algebra = self.withoutShortRelations()
        relations = algebra.relations()
        if not relations:
            return None
        (start, arrows), rest = relations[0], relations[1:]
        if arrows - 2 < 1:
            return None
        return algebra._withRelations([(arrows - 2, start + 2)] + rest)

    def swapLastRelation(self):
        """Operation 1 at the last relation, or None where it does not apply.

        The mirror of `swapFirstRelation`, at the other end of the main string:
        k_{r+1} = n - n_r - l_r + 1 is the main string after the last foot and
        m_r = l_r - 2 is the cord there, so exchanging them replaces the last
        relation (n_r, l_r) by (n_r, n - n_r - l_r + 3):

            the relation's START vertex does not move; its length becomes what
            was left of the line beyond it.

        Note that it is the *start* that is fixed here and the *end* in
        `swapFirstRelation`, not the other way round -- the operation is the
        exchange of two branches at a foot, and the foot is what stays put.
        """
        algebra = self.withoutShortRelations()
        relations = algebra.relations()
        if not relations:
            return None
        rest, (start, arrows) = relations[:-1], relations[-1]
        swapped = algebra.length - start - arrows + 3
        if swapped < 3:
            return None
        return algebra._withRelations(rest + [(start, swapped)])

    def classPreservingOrbit(self):
        """Every LNA the paper's operations reach from this one, as a frozenset.

        The closure of `withoutShortRelations` under `swapFirstRelation`,
        `swapLastRelation` and `relationDual`.  For an algebra with almost
        separate relations this is the paper's list of at most eight algebras in
        the class, and it must coincide with the set of LNAs `quipu` sends to the
        same quipu -- which is checked exhaustively for n <= 10 in
        tests/test_quipu_symmetry.py.

        For an algebra outside the theorem the operations still return something,
        but the paper claims nothing about it, so neither does this.
        """
        orbit = {self.withoutShortRelations()}
        frontier = list(orbit)
        while frontier:
            algebra = frontier.pop()
            for operation in ('swapFirstRelation', 'swapLastRelation', 'relationDual'):
                reached = getattr(algebra, operation)()
                if reached is not None and reached not in orbit:
                    orbit.add(reached)
                    frontier.append(reached)
        return frozenset(orbit)

    def _withRelations(self, relations):
        """A sibling algebra carrying the given (start, arrows) relations.

        None when that list is not an admissible ideal on this many vertices --
        two relations starting at the same vertex, one running off the end, or a
        pair violating the paper's standing assumption n_i + l_i < n_{i+1} +
        l_{i+1}.  The swaps can produce any of those at the edge of their range,
        and a rejected swap is not a derived equivalence.
        """
        relations = sorted(relations)
        starts = [start for start, _ in relations]
        if len(set(starts)) != len(starts):
            return None
        if any(start < 1 or arrows < 2 or start + arrows > self.length
               for start, arrows in relations):
            return None
        if any(start + arrows >= nextStart + nextArrows
               for (start, arrows), (nextStart, nextArrows) in zip(relations, relations[1:])):
            return None
        lengths = [0] * (self.length - 2)
        for start, arrows in relations:
            lengths[start - 1] = arrows
        return LinearNakayamaAlgebra(self.length, lengths)

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

    # -- enumeration -------------------------------------------------------

    @classmethod
    def allOfLength(cls, length):
        """Every LNA of the given length, in the table's order."""
        return [
            cls.fromRelationString(length, lines.relSetToString(relSet))
            for relSet in lines.generateAllPossibleLineRelations(length)
        ]


class QuipuAlgebra(pathAlgebra.PathAlgebra):
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

    # -- the notation's ambiguity, as operations -----------------------

    def exchangeAtFirstFoot(self):
        """The same quipu with k_0 and m_0 exchanged.

        A different name for the same tree, so a derived equivalent algebra.
        `LinearNakayamaAlgebra.swapFirstRelation` is what this does to the LNA
        the quipu corresponds to.
        """
        return QuipuAlgebra(*quipuForms.exchangeAtFirstFoot(self.k, self.m))

    def exchangeAtLastFoot(self):
        """The same quipu with k_{r+1} and m_r exchanged."""
        return QuipuAlgebra(*quipuForms.exchangeAtLastFoot(self.k, self.m))

    def endExchanges(self):
        """Both end exchanges, trivial ones included.

        Every one of these has the same underlying tree as `self`, hence the
        same canonical parameters; the point of having them is that the LNA-side
        operations can be checked against them.
        """
        return [self.exchangeAtFirstFoot(), self.exchangeAtLastFoot()]

    def __repr__(self):
        return "QuipuAlgebra({0}, {1})".format(self.k, self.m)

    def __eq__(self, other):
        return isinstance(other, QuipuAlgebra) and (self.k, self.m) == (other.k, other.m)

    def __hash__(self):
        return hash((self.k, self.m))
