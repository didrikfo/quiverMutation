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
    nextLabel = itertools.count()
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
