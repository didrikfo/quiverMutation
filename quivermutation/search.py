"""Walking the mutation graph.

`mutationSearchDepthFirst` is the one search primitive: descend through
admissible right mutations to a bounded depth, recording every quiver reached
that is again a line, and every quiver reached with no relations left.  The
second kind identifies the class completely, since for hereditary algebras of
tree type the underlying undirected tree is the whole of the class -- as a
derived class, and, because the orientations are joined by reflections that are
themselves mutations (`reflections`, F-036), as a mutation class too.  That is
what `hereditaryFormsReachedFrom` collects.

Reachability here is one-way: the search only walks right mutations, so A can
reach B at a depth where B reaches nothing.  Searching from the relation dual as
well is what covers the other direction -- see `memberAndItsDual`.
"""

import contextlib
import copy

import networkx as nx

from . import arrowPaths
from . import invariants
from . import lines
from . import mutation
from . import nakayama
from . import pathAlgebra
from . import procedure
from . import quipuForms
from . import reduction


def _coxeterKeyOrNone(pathAlg):
    """`coxeterKey`, or None where there is no invariant to be had.

    `coxeterCoefficients` raises when the Cartan matrix is not unimodular, which
    is what a quiver with an oriented cycle gives -- there are infinitely many
    paths and the matrix does not mean what the identity needs it to mean.
    """
    try:
        return invariants.coxeterKey(pathAlg)
    except (ValueError, ZeroDivisionError):
        return None


# -- conditional deeper probing ------------------------------------------
#
# A depth-bounded search spends its budget uniformly: every branch gets the
# same number of mutations, whether it is passing through quivers that look
# like nothing in particular or through the few that look like something.  The
# quivers with parallel arrows are the current example -- the procedure only
# learned to mutate into them in F-039, they are rare, and the one exit traced
# by hand out of `3030` takes six further mutations to come back to a quiver
# without them, which a run that stops at depth 8 will never see.  Raising the
# depth for *every* branch to reach them costs about fivefold per level.
#
# A `DeeperWhen` raises it only where a condition holds.  The same object with
# `extraDepth = 0` is a pure recorder, which is the other way of asking the
# question -- note down every interesting quiver a pass goes through and search
# from those afterwards -- so the two are one mechanism here, not two.  They are
# not equivalent, and the difference is the budget: see `DeeperWhen`.


def hasParallelArrows(pathAlg):
    """Two distinct arrows sharing both endpoints."""
    return arrowPaths.hasParallelArrows(pathAlg.quiver)


def parallelArrowsAtLeast(count):
    """At least `count` arrows beyond the ones a simple quiver would have."""
    def condition(pathAlg):
        return countParallelArrows(pathAlg.quiver) >= count
    condition.__name__ = 'parallelArrowsAtLeast({0})'.format(count)
    return condition


def hasNoRelations(pathAlg):
    """The algebra is hereditary -- its derived class is settled outright."""
    return not bool(pathAlg.rels)


def hasOrientedCycle(pathAlg):
    """The quiver has an oriented cycle.

    **This one cannot deepen anything.**  `mutationSearchDepthFirst` does not
    descend from a cyclic quiver at all, so a node where this fires has no
    children to spend extra depth on whatever the grant says.  It is here to be
    recorded with, and is registered so that a run can ask how often the search
    walks into one.
    """
    # Lazily, unlike the `noCycles` the walk computes: this is asked at every
    # node and the first cycle is the whole answer.
    return next(nx.simple_cycles(pathAlg.quiver), None) is not None


def anyOf(*conditions):
    """True where any of the conditions is."""
    def condition(pathAlg):
        return any(each(pathAlg) for each in conditions)
    condition.__name__ = 'anyOf({0})'.format(
        ', '.join(getattr(each, '__name__', '?') for each in conditions))
    return condition


def allOf(*conditions):
    """True where all of the conditions are."""
    def condition(pathAlg):
        return all(each(pathAlg) for each in conditions)
    condition.__name__ = 'allOf({0})'.format(
        ', '.join(getattr(each, '__name__', '?') for each in conditions))
    return condition


DEEPER_CONDITIONS = {
    'parallel-arrows': hasParallelArrows,
    'no-relations': hasNoRelations,
    'oriented-cycle': hasOrientedCycle,
}


def countParallelArrows(quiver):
    """How many arrows the underlying simple *directed* graph would lose.

    Not what `describeRelationFreeQuiver` counts: that one subtracts the edges
    of the underlying *undirected* graph, which also merges a pair of opposite
    arrows.  Here a 2-cycle is two arrows, as it should be -- they are not
    parallel.
    """
    return quiver.number_of_edges() - len({(tail, head)
                                           for tail, head, _ in quiver.edges(keys = True)})


def describeNode(pathAlg, mutationVertices = None):
    """The cheap summary of a node, for recording a probe firing.

    Deliberately much less than `describeRelationFreeQuiver`: this runs at every
    firing of a condition that may hold at thousands of nodes, and canonical
    forms are not free.
    """
    return {
        'path': list(mutationVertices or []),
        'vertices': pathAlg.quiver.number_of_nodes(),
        'arrows': pathAlg.quiver.number_of_edges(),
        'parallelArrows': countParallelArrows(pathAlg.quiver),
        'relations': len(pathAlg.rels),
    }


class DeeperWhen:
    """Extra depth for the branches that reach a quiver worth looking at.

    `condition(pathAlg)` is asked at every node the search reaches, before the
    node's children are walked.  Where it holds, the branch is given up to
    `extraDepth` more mutations -- and `budget` is what keeps that finite.

    **The budget is the whole safety argument.**  A grant that renewed at every
    node where the condition holds would not terminate: parallel arrows tend to
    beget parallel arrows, so a branch inside that region would refill its depth
    faster than it spent it.  `budget` caps the total extra depth one branch may
    ever accumulate, so no branch runs longer than `depth + budget` mutations
    and the search is as finite as it was.  It defaults to `extraDepth`, which
    means the grant is made once: the first quiver on a branch that meets the
    condition buys the extra depth, and re-entering the region later on the same
    branch buys nothing.  Pass a larger budget to allow that.

    `limit` caps the number of *grants* in the whole run, across branches, as a
    valve for an overnight job: past it the condition still records but buys
    nothing.  Firings are recorded either way, in `firings`, in visit order,
    with `keepQuivers` deciding whether each record carries a copy of the
    algebra as well as its summary -- copies are what makes a second pass
    possible and are also what makes a long run large, so they are off by
    default.

    With `extraDepth = 0` nothing is granted and the object is a recorder; see
    `recordOnly`.  Recording and re-searching afterwards is not quite the same
    search, and it is the wider of the two:

        one pass at depth d with budget b  <=  union over the firings a plain
        pass at depth d records of a plain search to that firing's remaining
        depth + b

    because every recorded firing starts its own fresh budget, where one pass
    spends a single budget along the whole branch.  The containment can only be
    strict for a branch that *leaves* the region and comes back: a firing below
    the one that bought the depth is already inside the subtree the grant paid
    for, at exactly the remaining depth re-searching it would give.  At every
    size measured the two come out equal (E-036), so the difference is not what
    to choose between them on -- the second's advantage is that the count of
    firings is visible before the deeper round is paid for.
    """

    def __init__(self, condition, extraDepth = 2, budget = None, limit = None,
                 keepQuivers = False, name = None):
        if extraDepth < 0:
            raise ValueError("extraDepth is a grant, not a penalty: {0}".format(extraDepth))
        self.condition = condition
        self.extraDepth = extraDepth
        self.budget = extraDepth if budget is None else budget
        if self.budget < 0:
            raise ValueError("a negative budget is not a budget: {0}".format(self.budget))
        self.limit = limit
        self.keepQuivers = keepQuivers
        self.name = name or getattr(condition, '__name__', 'condition')
        self.firings = []
        self.grants = 0

    def grant(self, pathAlg, mutationVertices, depth, spent):
        """(depth, spent) for a node the search has reached, after any grant."""
        if not self.condition(pathAlg):
            return depth, spent
        allowed = max(0, min(self.extraDepth, self.budget - spent))
        if self.limit is not None and self.grants >= self.limit:
            allowed = 0
        firing = describeNode(pathAlg, mutationVertices)
        firing['depth'] = depth
        firing['granted'] = allowed
        if self.keepQuivers:
            firing['pathAlg'] = copy.deepcopy(pathAlg)
        self.firings.append(firing)
        if allowed:
            self.grants += 1
        return depth + allowed, spent + allowed

    def spec(self):
        """The `deeperWhenFromSpec` string for this probe.

        What a checkpoint records, so that a run resumed under a *different*
        condition does not skip the work the condition would have changed.  It
        round-trips only for a registered condition -- a probe built around a
        function of one's own has that function's name here and
        `deeperWhenFromSpec` will not know it, which is a reason to keep
        `DEEPER_CONDITIONS` the vocabulary a long run is asked for in.
        """
        fields = [self.name, str(self.extraDepth), str(self.budget)]
        if self.limit is not None:
            fields.append(str(self.limit))
        return ':'.join(fields)

    def summarise(self):
        """What a run's firings came to."""
        return {
            'condition': self.name,
            'extraDepth': self.extraDepth,
            'budget': self.budget,
            'firings': len(self.firings),
            'grants': self.grants,
            'depthGranted': sum(firing['granted'] for firing in self.firings),
            'shallowest': min((len(f['path']) for f in self.firings), default = None),
            'deepest': max((len(f['path']) for f in self.firings), default = None),
        }


def recordOnly(condition, keepQuivers = True):
    """A `DeeperWhen` that grants nothing and only writes down what it sees.

    This is the other half of the mechanism: run a pass at the depth you were
    going to run anyway, see how many interesting quivers it went through and
    what they are, and only then decide what a deeper pass from those is worth.
    `keepQuivers` is on here because a record you cannot search from again is
    the one thing this is for.
    """
    return DeeperWhen(condition, extraDepth = 0, keepQuivers = keepQuivers)


def deeperWhenFromSpec(spec):
    """Parse 'parallel-arrows', 'parallel-arrows:3', 'parallel-arrows:3:6:100'.

    The fields are condition, extraDepth, budget, limit, in that order, and any
    tail of them may be left off.  This is what lets a run be asked for from a
    command line; `DEEPER_CONDITIONS` is the set of names.
    """
    fields = spec.split(':')
    name = fields[0]
    if name not in DEEPER_CONDITIONS:
        raise ValueError("no such condition {0!r}; known: {1}".format(
            name, ', '.join(sorted(DEEPER_CONDITIONS))))
    numbers = []
    for field in fields[1:]:
        try:
            numbers.append(int(field))
        except ValueError:
            raise ValueError("{0!r} is not a number, in {1!r}".format(field, spec))
    if len(numbers) > 3:
        raise ValueError("at most condition:extraDepth:budget:limit, in {0!r}".format(spec))
    extraDepth = numbers[0] if len(numbers) > 0 else 2
    budget = numbers[1] if len(numbers) > 1 else None
    limit = numbers[2] if len(numbers) > 2 else None
    return DeeperWhen(DEEPER_CONDITIONS[name], extraDepth = extraDepth,
                      budget = budget, limit = limit, name = name)


def mutationSearchDepthFirst(pathAlg, depth, mutationVertices = None, quiverName = 'quiver', vertexRelabeling = None, printOutput = True, collected = None, collectedHereditary = None, visitor = None, coxeterGuard = True, baseKey = None, deeperWhen = None, deeperSpent = 0):
    """Walk mutations of pathAlg to the given depth, recording the lines found.

    Every quiver reached that is again a line is recorded as a triple
    (path algebra, mutation path, vertex numbering).  Pass a list as `collected`
    to receive those triples in memory, in the order the search visits them.

    Pass a list as `collectedHereditary` to also receive, for every quiver
    reached that has no relations left, a triple (canonical form of the
    underlying undirected graph, the quipu notation for it where it applies,
    mutation path).  Those are the hereditary algebras in the class, and they
    identify it completely.

    Pass `visitor` to see *every* quiver the search reaches, line or not: it is
    called as `visitor(pathAlg, mutationVertices)` at each node, before the
    node's children are walked.  The two collectors above are the two questions
    asked often enough to have been built in; a visitor is for the rest, such as
    asking which quivers of a given shape a class passes through.

    `deeperWhen` is a `DeeperWhen`: a condition that buys the branches reaching
    a quiver worth a longer look more depth than the rest of the walk gets, and
    records every node it fires at.  `deeperSpent` carries how much of its
    budget the branch has already taken down the recursion, and is not for
    callers to pass.  See the class for why a budget is what makes the deeper
    walk terminate.

    `quiverName` only labels the progress output.  It used to name a
    '<quiverName>DF.txt' transcript that the caller parsed back by string
    slicing; that round trip is gone -- see NOTES.md idea 11.

    **`coxeterGuard` is what makes the walk a walk in one derived class.**
    `mutationIsPossibleAtVertex` rules mutation *out*, not in -- the paper's own
    hypothesis is on the algebra and it says plainly that this is in general not
    equivalent to a condition on the quiver -- so a step it admits can still fail
    to be a derived equivalence.  R-005 recorded that for rule discovery and made
    `lnaMoves.verifyMove` require three things: the predicted result, every step
    admissible, and the Coxeter polynomial unchanged.  This search asked only for
    the second, and F-038 is what that let through: 97 quivers at `n = 7` alone,
    acyclic and with no parallel arrows, whose Coxeter polynomial has moved and
    which the search then walks straight on from.  The guard is the third
    requirement, applied per step: a mutation whose `coxeterKey` differs from the
    start's is not taken.  `baseKey` carries the start's key down the recursion
    and is computed here when the caller does not supply it.

    Passing `coxeterGuard = False` restores the old behaviour, and is for
    measuring what the guard changes, not for producing answers.
    """
    # These used to default to [] and {}, which Python evaluates once at
    # definition time.  The relabeling dict is filled in below and so leaked
    # between searches: a second search in the same process inherited the
    # first one's numbering, and crashed as soon as the quiver was longer.
    mutationVertices = [] if mutationVertices is None else mutationVertices
    vertexRelabeling = {} if vertexRelabeling is None else dict(vertexRelabeling)
    if coxeterGuard and baseKey is None:
        baseKey = _coxeterKeyOrNone(pathAlg)
        if baseKey is None:
            # No usable invariant to compare against -- a cyclic quiver has no
            # unimodular Cartan matrix.  Nothing to guard with, so do not.
            coxeterGuard = False
    vertices = list(pathAlg.vertices())
    baseQuiver = copy.deepcopy(pathAlg.quiver)
    quiverAtThisDepth = copy.deepcopy(pathAlg.quiver)
    rels = copy.deepcopy(pathAlg.rels)
    relsAtThisDepth = copy.deepcopy(pathAlg.rels)
    if not bool(vertexRelabeling):
        for vertex in vertices:
            vertexRelabeling[vertex] = vertex
    longestPathLength = 0
    noCycles = not bool(list(nx.simple_cycles(baseQuiver)))
    if noCycles:
        longestPathLength = nx.dag_longest_path_length(baseQuiver)
    if printOutput:
        print('Quiver name: ', quiverName)
        print("Mutations: ", mutationVertices)
        print('Numbering: {0}'.format(vertexRelabeling))
        print("Longest path: ", longestPathLength)
        pathAlgebra.printPathAlgebra(pathAlg)
    if not bool(rels) and (collectedHereditary is not None or _SIGHTING_SINKS):
        # No relations left: the algebra is hereditary, and the underlying
        # undirected graph of its quiver is a complete derived invariant.
        graph = quipuForms.underlyingGraph(pathAlg)
        if collectedHereditary is not None:
            collectedHereditary.append((
                quipuForms.canonicalUndirectedForm(graph),
                quipuForms.formatQuipu(quipuForms.quipuParameters(graph)),
                mutationVertices[:],
            ))
        for sink in _SIGHTING_SINKS:
            sink.append(describeRelationFreeQuiver(pathAlg, mutationVertices, graph))
    if visitor is not None:
        visitor(pathAlg, mutationVertices)
    isLine = (longestPathLength == len(vertices) - 1) and (len(baseQuiver.edges) == len(vertices) - 1)
    if isLine and collected is not None:
        foundPathAlg = pathAlgebra.PathAlgebra()
        foundPathAlg.quiver = baseQuiver
        foundPathAlg.rels = rels
        collected.append((foundPathAlg, mutationVertices[:], dict(vertexRelabeling)))
    # debugVertexList = [1, 1, 2, 1, 2, 3, 5, 3, 4, 4, 5, 2, 2, 3, 6, 1, 4, 1, 2, 3, 1, 4]
    # for i in range(7, len(debugVertexList)):
    #      if mutationVertices == debugVertexList[:i]:
    #          input('Press enter to continue...')
    if deeperWhen is not None:
        # Before the depth test, not after: the node this is for is typically
        # the one the walk has just run out of budget at.
        depth, deeperSpent = deeperWhen.grant(pathAlg, mutationVertices, depth, deeperSpent)
    if depth > 0 and noCycles:
        depth = depth - 1
        for vertex in reversed(vertices):
            discardMutation = False
            pathAlg.quiver = copy.deepcopy(quiverAtThisDepth)
            pathAlg.rels = copy.deepcopy(relsAtThisDepth)
            mutationVerticesAtDepth = mutationVertices[:]
            # `mutationIsPossibleAtVertex` is the whole gate.  It used to be
            # followed here by a second test of the same idea, counting the
            # paths into the vertex against the paths through each arrow out of
            # it and refusing the vertex if any one arrow lost a path.  That is
            # the *strict* reading -- every arrow must keep every path -- where
            # the paper's theorem rules mutation out only when a path dies
            # against all of them, so the search was refusing mutations the
            # paper allows, twice over.  F-016.
            if mutation.mutationIsPossibleAtVertex(pathAlg, vertex):
                mutationVerticesAtDepth.append(vertexRelabeling[vertex])
                mutPathAlg = mutation.quiverMutationAtVertex(pathAlg, vertex)
                # The check is over the *arrow* relations, not `rels`.  Two
                # parallel paths write down as the same vertex sequence, so
                # `paths.isIllegalRelation` called a commutativity relation
                # between them a repeated path and discarded a mutation that is
                # perfectly legal -- which is how the parallel-arrow branches
                # used to die even before the gate refused them.
                for rel in procedure.relationsFrom(mutPathAlg):
                    if arrowPaths.isIllegalRelation(mutPathAlg.quiver, rel):
                        print('ILLEGAL RELATION!')
                        print('The relation {0}'.format(arrowPaths.projectToPathSets([rel])))
                        print('is illegal in the following path algebra:')
                        print(arrowPaths.describeQuiver(mutPathAlg.quiver,
                                                        procedure.relationsFrom(mutPathAlg)))
                        discardMutation = True
                        break
                if discardMutation:
                    # `continue`, not `break`: only this vertex is discarded.
                    # It used to break the loop over vertices, so one illegal
                    # relation abandoned every vertex still to be tried at this
                    # node -- and since the loop runs in reverse, that was every
                    # lower-numbered one.  E-033.
                    continue
                mutPathAlg = reduction.reducePathAlgebra(mutPathAlg)
                if coxeterGuard:
                    movedKey = _coxeterKeyOrNone(mutPathAlg)
                    if movedKey is not None and movedKey != baseKey:
                        # Admissible but not a derived equivalence.  See the note
                        # on `coxeterGuard` above, and F-038.
                        continue
                    # movedKey is None for a cyclic quiver, which the search does
                    # not descend from anyway; it is let through so that cycles
                    # end a branch exactly as they did before the guard.
                mutationSearchDepthFirst(copy.deepcopy(mutPathAlg), depth, mutationVerticesAtDepth, quiverName, vertexRelabeling, printOutput, collected, collectedHereditary, visitor, coxeterGuard, baseKey, deeperWhen, deeperSpent)
    return


def hereditaryFormsReachedFrom(pathAlg, depth, deeperWhen = None):
    """The hereditary algebras reachable from pathAlg within `depth` mutations.

    Returns a dict mapping the canonical form of the underlying undirected graph
    to (quipu notation, shortest mutation path found to it).  An empty result
    means no relation-free quiver was reached at this depth, not that none
    exists.
    """
    found = []
    mutationSearchDepthFirst(pathAlg, depth, [], 'hereditary', printOutput = False,
                             collected = None, collectedHereditary = found,
                             deeperWhen = deeperWhen)
    forms = {}
    for canonical, quipu, path in found:
        if canonical not in forms or len(path) < len(forms[canonical][1]):
            forms[canonical] = (quipu, path)
    return forms


def findHereditaryFormForClass(table, lineLength, className, maxDepth = 8, printOutput = False,
                               deeperWhen = None):
    """Search the members of one class for a relation-free quiver.

    Iterative deepening from each member in turn, returning as soon as any
    member reaches one.  Members are tried shortest-relation-string first, on
    the observation that an LNA with fewer and shorter relations tends to need
    fewer mutations to shed them all.

    Each member is searched from itself and from its relation dual.  The search
    only walks right mutations, so reachability is one-way; reversing every arrow
    is one of the class-preserving operations of arXiv:2305.06642, and
    rightMutate(dual(P)) = dual(leftMutate(P)), so a right-mutation path out of
    the dual is a left-mutation path out of the member and everything it reaches
    is still in the class.  The dual of a tree quiver is the same tree, so a form
    reached from the dual is the class' form unchanged.

    Returns the hereditary form as it should be written into the table -- the
    quipu notation where the graph is a quipu, otherwise the canonical tree
    encoding -- or '' if nothing was reached.
    """
    members = sorted(table.membersOfClass(className), key = lambda r: (len(r), r))
    for depth in range(2, maxDepth + 1):
        for relationString in members:
            for startPoint in memberAndItsDual(lineLength, relationString):
                forms = hereditaryFormsReachedFrom(startPoint, depth, deeperWhen)
                if forms:
                    if printOutput:
                        print('class {0} reaches {1} at depth {2} from {3!r}'.format(
                            className, sorted(forms), depth, relationString))
                    return formatHereditaryForms(forms)
    return ''


def memberAndItsDual(lineLength, relationString):
    """An LNA and its relation dual, both as path algebras, without repeats."""
    algebra = nakayama.LinearNakayamaAlgebra.fromRelationString(lineLength, relationString)
    dual = algebra.relationDual()
    return [algebra] if dual == algebra else [algebra, dual]


def formatHereditaryForms(forms):
    """Render the result of hereditaryFormsReachedFrom for the table.

    All the hereditary algebras in one derived equivalence class have isomorphic
    underlying graphs, so this is normally one value.  More than one would mean
    either a bug in the mutation procedure or a non-tree in the mix, so they are
    all reported, joined by '|', rather than silently reduced to one.
    """
    return '|'.join(sorted(quipu or canonical for canonical, (quipu, _path) in forms.items()))


def hereditaryFormFromTheorem(lineLength, relationString):
    """The hereditary form of one LNA, straight from the quipu theorem.

    Returns '' when the LNA does not have almost separate relations, since
    theorem `thm:QuipuToAn` of arXiv:2305.06642 says nothing about those.  For
    the ones it does cover this is O(1), where reaching the same answer by
    mutation search costs a depth-4-to-9 traversal.
    """
    algebra = nakayama.LinearNakayamaAlgebra.fromRelationString(lineLength, relationString)
    return algebra.quipuName()


def linesReachedFrom(pathAlg, depth, alsoDual = True, deeperWhen = None):
    """The LNAs a bounded mutation search out of `pathAlg` reaches.

    Returns a dict from relation string to the shortest mutation path found to
    it.  This is what settles a Coxeter-polynomial lead: the polynomial says two
    algebras *could* be derived equivalent, and a mutation path from one to the
    other says they are.

    That last step used to be justified here by "every step of the procedure is a
    tilting mutation", which is true of the procedure and was **not** true of the
    walk this function performs -- R-012.  It holds now because
    `mutationSearchDepthFirst` guards every step on the Coxeter polynomial as
    well as on admissibility (F-038).  What the guard gives is still a necessary
    condition rather than a proof; H-015 is whether it is also sufficient, and a
    link that matters is worth replaying and checking rather than trusting.

    With `alsoDual`, the search is run from the opposite algebra as well, whose
    lines are read back through the dual.  Two algebras are derived equivalent
    exactly when their opposites are, and the opposite of an LNA is an LNA, so a
    line reached from the opposite is as good a witness as one reached directly
    -- and it is the only way to see what a *left* mutation path would reach,
    since the search walks right mutations only.

    `deeperWhen` is passed straight through, and is shared by both searches:
    what it records is the whole run, and its per-branch budget starts afresh at
    each of the two start points, as a branch of one is not a branch of the
    other.
    """
    startPoints = [pathAlg]
    if alsoDual:
        startPoints.append(pathAlgebra.dualPathAlgebra(pathAlg))
    reached = {}
    for index, startPoint in enumerate(startPoints):
        collected = []
        mutationSearchDepthFirst(copy.deepcopy(startPoint), depth, [], 'lines',
                                 printOutput = False, collected = collected,
                                 deeperWhen = deeperWhen)
        for found in lines.mutationListLineCleanup(collected, printOutput = False):
            algebra = found[0] if index == 0 else pathAlgebra.dualPathAlgebra(found[0])
            # Sorted, because a relation set read off the dual comes out in the
            # reverse order and `relSetToString` writes it down as it stands --
            # which would make one LNA look like two.
            relationString = lines.relSetToString(sorted(_renumberedLine(algebra).rels))
            path = found[1]
            if relationString not in reached or len(path) < len(reached[relationString]):
                reached[relationString] = path
    return reached


def _renumberedLine(pathAlg):
    """A line algebra with its vertices renumbered 1 -> ... -> n along the line.

    `mutationListLineCleanup` already does this for what it collects; taking the
    opposite afterwards reverses the numbering, so it has to be done again.
    """
    order = nx.topological_sort(pathAlg.quiver)
    relabeling = {vertex: position for position, vertex in enumerate(order, start = 1)}
    renumbered = pathAlgebra.PathAlgebra()
    renumbered.add_vertices_from(sorted(relabeling.values()))
    for tail, head in pathAlg.quiver.edges():
        renumbered.add_arrow(relabeling[tail], relabeling[head])
    for rel in pathAlg.rels:
        renumbered.add_rel([[relabeling[vertex] for vertex in path] for path in rel])
    return renumbered


# -- relation-free sightings ---------------------------------------------
#
# Relations are what the procedure spends its time on, so a mutation that leaves
# a quiver with *none* is an event: the algebra is hereditary, and the underlying
# graph of the quiver settles its derived equivalence class outright.  The
# searches already use that -- `hereditaryFormsReachedFrom` is nothing else --
# but they use it locally and throw the rest away, which means nobody has ever
# looked at what those quivers are.  The expectation is that every one is a tree
# and almost every one a quipu; a relation-free quiver whose graph has a *cycle*
# would be a hereditary algebra of a kind no LNA class has produced, and would be
# worth stopping for.  Neither is known, and the cost of finding out is a few
# lines, since every search already passes through the place where it could be
# recorded.
#
# Sightings are opt-in: with no sink open, the branch above does no more work
# than it did before.

_SIGHTING_SINKS = []


@contextlib.contextmanager
def relationFreeSightings():
    """Record every relation-free quiver the searches inside the block reach.

    Yields the list the sightings accumulate in, one dict per sighting, in the
    order the searches visit them.  Sinks nest, so an inner block does not stop
    an outer one from seeing what it sees.

        with search.relationFreeSightings() as sightings:
            classification.classifyLength(9)
        print(search.summariseSightings(sightings))
    """
    sink = []
    _SIGHTING_SINKS.append(sink)
    try:
        yield sink
    finally:
        _SIGHTING_SINKS.remove(sink)


def describeRelationFreeQuiver(pathAlg, mutationVertices = None, graph = None):
    """What a relation-free quiver is, in the terms worth counting.

    `isTree` and `isQuipu` are about the *underlying undirected* graph, since
    that is what decides a hereditary algebra's derived equivalence class.
    `hasOrientedCycle` is about the quiver, and `parallelArrows` counts arrows
    the underlying simple graph merges -- either of those would mean the algebra
    is not the path algebra of a tree, whatever the undirected picture says.
    """
    graph = quipuForms.underlyingGraph(pathAlg) if graph is None else graph
    isTree = nx.is_tree(graph) if graph.number_of_nodes() else False
    return {
        'canonical': quipuForms.canonicalUndirectedForm(graph),
        'quipu': quipuForms.formatQuipu(quipuForms.quipuParameters(graph)),
        'vertices': graph.number_of_nodes(),
        'arrows': pathAlg.quiver.number_of_edges(),
        'isTree': isTree,
        'isQuipu': bool(isTree and quipuForms.isQuipuByDegrees(graph)),
        'isConnected': bool(graph.number_of_nodes()) and nx.is_connected(graph),
        'hasOrientedCycle': bool(list(nx.simple_cycles(pathAlg.quiver))),
        'parallelArrows': pathAlg.quiver.number_of_edges() - graph.number_of_edges(),
        'maxDegree': max((degree for _, degree in graph.degree()), default = 0),
        'path': list(mutationVertices or []),
    }


def summariseSightings(sightings):
    """Counts of the kinds of relation-free quiver a run saw.

    `quipus`, `otherTrees` and `notTrees` partition the sightings; `distinct` is
    how many isomorphism classes of underlying graph they came to, which is the
    number that says whether a run saw one thing many times or many things.
    `oddities` is the sightings that are not trees, kept in full, since those are
    the ones there is no reason to expect.
    """
    oddities = [sighting for sighting in sightings if not sighting['isTree']]
    return {
        'sightings': len(sightings),
        'distinct': len({sighting['canonical'] for sighting in sightings}),
        'quipus': sum(1 for sighting in sightings if sighting['isQuipu']),
        'otherTrees': sum(1 for sighting in sightings
                          if sighting['isTree'] and not sighting['isQuipu']),
        'notTrees': len(oddities),
        'oddities': oddities,
    }
