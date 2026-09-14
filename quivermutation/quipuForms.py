"""Canonical forms for the hereditary algebras a mutation search reaches.

When a mutation search leaves a quiver with no relations, the algebra is
hereditary, and for hereditary algebras derived equivalence is well understood:
two path algebras of trees are derived equivalent exactly when the trees are
isomorphic as undirected graphs, since the orientations are related by BGP
reflections.  So the underlying undirected graph of any relation-free quiver a
search reaches is a complete derived invariant of the class, and a far sharper
one than the Coxeter polynomial.

That gives a way to separate two classes that share a Coxeter polynomial: if
they reach relation-free quivers whose underlying graphs are not isomorphic,
they are not derived equivalent, whatever the polynomial says.

For the linear Nakayama algebras the trees in question are quipus -- trees of
maximum degree 3 whose degree-3 vertices all lie on one path -- which is the
subject of arXiv:2305.06642.  `quipuParameters` recovers the paper's
P^(m_0..m_r)_(k_0..k_{r+1}) notation where it applies, in a canonical form,
since that notation does not determine the quipu on its own.
"""

import itertools

import networkx as nx


def underlyingGraph(pathAlg):
    """The quiver's underlying undirected simple graph."""
    graph = nx.Graph()
    graph.add_nodes_from(pathAlg.quiver.nodes)
    graph.add_edges_from((a[0], a[1]) for a in pathAlg.quiver.edges)
    return graph


def _rootedTreeForm(graph, root, parent):
    """The AHU canonical string of the subtree at root, coming from parent."""
    children = sorted(
        _rootedTreeForm(graph, child, root)
        for child in graph.neighbors(root)
        if child != parent
    )
    return "(" + "".join(children) + ")"


def canonicalTreeForm(graph):
    """A canonical string for an undirected tree, exact up to isomorphism.

    This is the AHU encoding rooted at the tree's centre.  A tree has one or two
    centres; with two, the smaller of the two encodings is taken, so the result
    depends only on the isomorphism class.
    """
    return min(_rootedTreeForm(graph, centre, None) for centre in nx.center(graph))


def canonicalUndirectedForm(graph):
    """A canonical key for the underlying graph of a hereditary quiver.

    Trees get their exact AHU encoding.  Anything else -- a quiver whose
    underlying graph has a cycle, which no LNA search has produced so far --
    falls back to a Weisfeiler-Leman hash, which is prefixed so the two can
    never be confused.  The hash is complete for trees but not in general, so a
    'wl:' key is evidence and not proof.
    """
    if graph.number_of_nodes() == 0:
        return "()"
    if nx.is_tree(graph):
        return canonicalTreeForm(graph)
    return "wl:" + nx.weisfeiler_lehman_graph_hash(graph)


def isQuipu(graph):
    """Whether the graph is an open quipu: a tree of maximum degree 3 whose
    degree-3 vertices all lie on a single path."""
    return quipuParameters(graph) is not None


def quipuParameters(graph):
    """The quipu's (k, m) parameters in the paper's notation, or None.

    Returns a pair of tuples ((k_0, ..., k_{r+1}), (m_0, ..., m_r)) for the
    quipu P^(m_0,...,m_r)_(k_0,...,k_{r+1}): m_i is the length of cord i, and
    k_i the number of main-string vertices between the foot of cord i-1 and the
    foot of cord i (with k_0 before the first foot and k_{r+1} after the last).

    The notation does not determine the quipu -- reading the main string from
    the other end, or swapping an end segment of the main string with the cord
    at the outermost foot, gives another valid parameter pair for the same graph
    -- so the lexicographically smallest pair over every valid reading is
    returned, making it canonical.
    """
    if graph.number_of_nodes() == 0 or not nx.is_tree(graph):
        return None
    degrees = dict(graph.degree())
    if max(degrees.values(), default=0) > 3:
        return None

    branchVertices = {v for v, d in degrees.items() if d == 3}
    leaves = [v for v, d in degrees.items() if d <= 1]
    if len(leaves) < 2:
        # A single vertex: the quipu with no cord and no main string beyond it.
        return ((0, 0), (0,))

    candidates = []
    # Both directions of each main string, since reading it from the other end
    # reverses k and m and can give the smaller parameter pair.  Using
    # combinations here made the result depend on how the vertices happened to
    # be labelled.
    for first, last in itertools.permutations(leaves, 2):
        mainString = nx.shortest_path(graph, first, last)
        if not branchVertices.issubset(mainString):
            continue
        parameters = _parametersAlongMainString(graph, mainString, branchVertices)
        if parameters is not None:
            candidates.append(parameters)
    return min(candidates) if candidates else None


def _parametersAlongMainString(graph, mainString, branchVertices):
    """(k, m) for one choice of main string, or None if the cords are not paths."""
    onMainString = set(mainString)
    cordLengths = []
    footIndices = []
    for index, vertex in enumerate(mainString):
        if vertex not in branchVertices:
            continue
        offMainString = [w for w in graph.neighbors(vertex) if w not in onMainString]
        if len(offMainString) != 1:
            return None
        length = _cordLength(graph, vertex, offMainString[0], onMainString)
        if length is None:
            return None
        cordLengths.append(length)
        footIndices.append(index)

    if not footIndices:
        # No cords at all: a path, written with one cord of length zero.
        return ((0, len(mainString) - 1), (0,))

    k = [footIndices[0]]
    for earlier, later in zip(footIndices, footIndices[1:]):
        k.append(later - earlier - 1)
    k.append(len(mainString) - 1 - footIndices[-1])
    return (tuple(k), tuple(cordLengths))


def _cordLength(graph, foot, first, onMainString):
    """The number of vertices hanging off `foot` through `first`, or None if
    that branch is not a bare path."""
    length = 0
    previous, current = foot, first
    while True:
        if current in onMainString:
            return None
        length += 1
        onward = [w for w in graph.neighbors(current) if w != previous]
        if not onward:
            return length
        if len(onward) > 1:
            return None
        previous, current = current, onward[0]


def formatQuipu(parameters):
    """'P^(m_0,...,m_r)_(k_0,...,k_{r+1})' for a parameter pair."""
    if parameters is None:
        return ""
    k, m = parameters
    return "P^({0})_({1})".format(",".join(map(str, m)), ",".join(map(str, k)))


def graphFromQuipuParameters(k, m):
    """Build the quipu P^(m_0,...,m_r)_(k_0,...,k_{r+1}) as an undirected graph.

    The inverse of quipuParameters up to isomorphism, so that a quipu named in
    the paper's notation can be canonicalised and compared with one a search
    found.
    """
    if len(k) != len(m) + 1:
        raise ValueError("a quipu needs one more k than m, got {0} and {1}".format(k, m))
    graph = nx.Graph()
    nextLabel = itertools.count(1)
    mainString = [next(nextLabel) for _ in range(k[0])]
    for cordIndex, cordLength in enumerate(m):
        foot = next(nextLabel)
        mainString.append(foot)
        previous = foot
        for _ in range(cordLength):
            current = next(nextLabel)
            graph.add_edge(previous, current)
            previous = current
        mainString.extend(next(nextLabel) for _ in range(k[cordIndex + 1]))
    graph.add_nodes_from(mainString)
    for earlier, later in zip(mainString, mainString[1:]):
        graph.add_edge(earlier, later)
    return graph


def canonicalQuipuParameters(k, m):
    """The canonical parameter pair for the quipu named by (k, m)."""
    return quipuParameters(graphFromQuipuParameters(k, m))


def quipuForAlmostSeparateLNA(lineLength, relLengths):
    """The quipu an LNA with almost separate relations is derived equivalent to.

    This inverts theorem `thm:QuipuToAn` of arXiv:2305.06642, which states that
    for a quipu P^(m_0,...,m_r)_(k_0,...,k_{r+1}) and

        n_i = k_0 + sum_{j=1..i} (m_{j-1} + k_j + 1)   for 1 <= i <= r+1,

    the quipu algebra is derived equivalent to

        A_{n_{r+1}, (k_0, n_1, ..., n_r)}^{(m_0+2, m_1+2, ..., m_r+2)}

    and to no quipu of any other shape.  Reading it backwards, an LNA
    A_{n,(n_0,...,n_r)}^{(l_0,...,l_r)} with almost separate relations, all of
    length >= 3, comes from

        m_i = l_i - 2,   k_0 = n_0,
        k_i = n_i - n_{i-1} - m_{i-1} - 1   for 1 <= i <= r,
        k_{r+1} = n - n_r - m_r - 1.

    Relations of length 2 do not change the derived equivalence class of such an
    algebra, so they are dropped first.

    `relLengths` is the per-vertex form used throughout the repo: entry i is the
    number of arrows in the relation starting at vertex i + 1, or 0.

    Returns canonical quipu parameters, or None when the algebra does not have
    almost separate relations, in which case the theorem says nothing about it.
    """
    relations = [(start + 1, length) for start, length in enumerate(relLengths) if length]
    if not _hasAlmostSeparateRelations(lineLength, relations):
        return None

    longRelations = [(start, length) for start, length in relations if length >= 3]
    if not longRelations:
        # Every relation has length 2, so the algebra is derived equivalent to
        # the path algebra of A_n itself.
        return canonicalQuipuParameters((0, lineLength - 1), (0,))

    starts = [start for start, _ in longRelations]
    m = [length - 2 for _, length in longRelations]
    k = [starts[0]]
    for i in range(1, len(starts)):
        k.append(starts[i] - starts[i - 1] - m[i - 1] - 1)
    k.append(lineLength - starts[-1] - m[-1] - 1)
    if any(value < 0 for value in k):
        return None
    return canonicalQuipuParameters(tuple(k), tuple(m))


def _hasAlmostSeparateRelations(lineLength, relations):
    """Whether consecutive relations overlap in at most one arrow.

    relations is a list of (start vertex, number of arrows), in increasing order
    of start vertex.  The paper's condition is n_{i+1} >= n_i + l_i - 1, along
    with the standing assumptions n_i < n_{i+1}, n_i + l_i < n_{i+1} + l_{i+1}
    and n_r + l_r <= n.
    """
    for (start, length), (nextStart, nextLength) in zip(relations, relations[1:]):
        if nextStart < start + length - 1:
            return False
        if start >= nextStart or start + length >= nextStart + nextLength:
            return False
    if relations and relations[-1][0] + relations[-1][1] > lineLength:
        return False
    return True


def isQuipuByDegrees(graph):
    """Whether the graph is a quipu, decided from its degrees.

    A tree of maximum degree 3 whose degree-3 vertices all lie on one path --
    read straight off the definition, since every path in a tree extends to one
    between two leaves.

    `isQuipu` answers the same question the other way round, by asking
    `quipuParameters` for a main string that accounts for every branch vertex
    and leaves bare paths hanging off it.  Two routes to one definition, so the
    two disagreeing would mean one of them is wrong -- which is what
    `quipusByTreeEnumeration` is for.
    """
    if graph.number_of_nodes() == 0 or not nx.is_tree(graph):
        return False
    degrees = dict(graph.degree())
    if max(degrees.values(), default=0) > 3:
        return False
    branchVertices = {v for v, d in degrees.items() if d == 3}
    if len(branchVertices) <= 1:
        return True
    leaves = [v for v, d in degrees.items() if d <= 1]
    return any(branchVertices.issubset(nx.shortest_path(graph, first, last))
               for first, last in itertools.combinations(leaves, 2))


def quipusByTreeEnumeration(order):
    """Every quipu of an order, reached by enumerating trees, not parameters.

    The independent route to `allQuipusOfOrder`: that one enumerates the
    P^(m)_(k) parameter pairs and canonicalises them, this one enumerates the
    non-isomorphic trees of the order and keeps the ones `isQuipuByDegrees`
    accepts.  Nothing about the notation enters the choice of which trees to
    keep, so agreement between the two checks both the parameter enumeration and
    the two readings of the definition.  They agree for orders 1 to 12 --
    research E-013.

    It is the more expensive route by far, since the number of trees grows much
    faster than the number of quipus, which is why it is the check and not the
    implementation.
    """
    found = set()
    for tree in nx.nonisomorphic_trees(order) if order >= 2 else [nx.empty_graph(1)]:
        if isQuipuByDegrees(tree):
            found.add(quipuParameters(tree))
    return sorted(found)


def allQuipusOfOrder(order):
    """Every quipu on `order` vertices, as canonical parameter pairs.

    A quipu P^(m_0,...,m_r)_(k_0,...,k_{r+1}) has r + 1 + sum(k) + sum(m)
    vertices, so enumerating r, then the k and m that fit, and canonicalising,
    gives each quipu exactly once.

    Counts: 1, 1, 1, 2, 2, 4, 6, 11, 18, 36, 64, 127 for orders 1 to 12, and the
    orders 6, 7 and 8 agree with the table in arXiv:2305.06642.
    """
    found = set()
    for cords in range(1, order + 1):
        # cords = r + 1, so there are `cords` cord lengths and cords + 1 gaps.
        budget = order - cords
        if budget < 0:
            continue
        for m in _compositions(cords, budget):
            for k in _compositions(cords + 1, budget - sum(m)):
                if sum(k) + sum(m) + cords != order:
                    continue
                canonical = quipuParameters(graphFromQuipuParameters(k, m))
                if canonical is not None:
                    found.add(canonical)
    return sorted(found)


def _compositions(parts, total):
    """Every tuple of `parts` non-negative integers summing to at most `total`."""
    if parts == 0:
        yield ()
        return
    for first in range(total + 1):
        for rest in _compositions(parts - 1, total - first):
            yield (first,) + rest


def adjacencySpectrumPolynomial(graph):
    """The characteristic polynomial of the undirected adjacency matrix.

    Two non-isomorphic trees with the same one are *cospectral*, and cospectral
    quipus are exactly where the Coxeter polynomial stops separating derived
    equivalence classes: the Coxeter polynomial of the path algebra of a tree is
    determined by the tree's spectrum, so cospectral quipus give algebras that
    are not derived equivalent yet share a Coxeter polynomial.
    """
    import sympy

    nodes = sorted(graph.nodes)
    matrix = sympy.Matrix(
        [[1 if graph.has_edge(u, v) else 0 for v in nodes] for u in nodes])
    return matrix.charpoly().as_expr()


def cospectralQuipuGroups(order):
    """Groups of two or more distinct quipus of an order that are cospectral.

    Each group is a set of derived equivalence classes that no Coxeter polynomial
    can tell apart, found without running a single mutation.  An empty result
    means the Coxeter polynomial separates every class of that order.
    """
    groups = {}
    for parameters in allQuipusOfOrder(order):
        graph = graphFromQuipuParameters(*parameters)
        groups.setdefault(adjacencySpectrumPolynomial(graph), []).append(parameters)
    return {
        polynomial: quipus
        for polynomial, quipus in groups.items()
        if len(quipus) > 1
    }
