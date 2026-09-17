"""The mutation procedure in the set-of-paths model.

The procedure itself is `procedure`, which runs on linear combinations of
paths; this module is its face for everything that speaks `PathAlgebra` and
`rels` -- the search, the tables, the naming.  `quiverMutationAtVertex` and
`leftQuiverMutationAtVertex` delegate,  `quiverMutationAtVertices` composes
them with `reduction.reducePathAlgebra` after each step, and results carry
their coefficients so a chain of mutations does not lose them.

`mutationIsPossibleAtVertex` is the admissibility condition, and is the one
thing here that is *not* delegated.  It is stricter than the criterion the paper
gives when the vertex has more than one arrow out of it: it rejects as soon as
one arrow out of the vertex kills a nonzero path, where the paper's theorem rules
mutation out only when *every* arrow does.  Since the paper is explicit that the
real condition, `Hom(P_i*[1], Lambda) = 0`, is in general **not** equivalent to
any condition on the quiver, neither reading is exact, and refusing too much is
the safe direction -- it loses reachability where being too permissive risks
performing a rewrite that is not a derived equivalence (research R-005).  So the
strict one stays the search's gate, and `procedure.isMutable` is the paper's
criterion read exactly, for when that is what is wanted.  F-015 measures the
difference.

A mutation applied outside the condition still returns a quiver, just not a
derived equivalent one -- research R-005 -- so callers that care about the class
must test admissibility first.
"""

import numpy as np

from . import invariants
from . import pathAlgebra
from . import plotting
from . import procedure
from . import reduction



def quiverMutationAtVertex(pathAlg, vertex):
    """Steps 1 to 7 of the procedure at `vertex`, without the cleanup after.

    The work is `procedure.mutateAtVertex`, which runs on linear combinations
    of paths; this is its face in the set-of-paths model the tables and the
    naming use.  The result carries its coefficients (`relCombinations`), so a
    chain of mutations does not lose them between steps.

    Does not check admissibility -- `mutationIsPossibleAtVertex` is that, and
    applying the rewrite without it returns a quiver that is not derived
    equivalent (research R-005).
    """
    return procedure.toPathAlgebra(*procedure.mutateAtVertex(
        pathAlg.quiver, procedure.relationsFrom(pathAlg), vertex))


def quiverMutationAtVertices(pathAlg, vertices : list, printMutationSteps = False):
    for i in range(len(vertices)):
        vertex = vertices[i]
        if vertex > 0:
            pathAlg = quiverMutationAtVertex(pathAlg, vertex)
        else:
            pathAlg = leftQuiverMutationAtVertex(pathAlg, -vertex)
        pathAlg = reduction.reducePathAlgebra(pathAlg)
        if printMutationSteps:
            print('Mutations: ', vertices[:i + 1])
            pathAlgebra.printPathAlgebra(pathAlg)
    return pathAlg


def leftQuiverMutationAtVertex(pathAlg, vertex):
    """Left mutation: the same procedure on the opposite algebra, read back."""
    return procedure.toPathAlgebra(*procedure.mutateLeftAtVertex(
        pathAlg.quiver, procedure.relationsFrom(pathAlg), vertex))


def mutationIsPossibleAtVertex(pathAlg, vertex):
    """Whether the mutation procedure may be applied to pathAlg at vertex.

    The paper's criterion, in `procedure.isMutable`; this is its face in the
    set-of-paths model.  It used to be a stricter reading of the same theorem --
    see F-016 for what the two differed by and what changing it changed.

    It rules mutation *out*, not in: the theorem's own hypothesis is on the
    algebra, and the paper says it is in general not equivalent to a condition
    on the quiver.  A rewrite performed on the strength of this can still fail
    to be a derived equivalence (research R-005), so a caller that cares about
    the class wants more than this alone.
    """
    return procedure.isMutable(pathAlg.quiver, procedure.relationsFrom(pathAlg), vertex)


def showMutationSteps(pathAlg, mutationVertexList, firstDisplayedStep = 0):
    """Mutate step by step, printing and plotting each one, for looking by hand.

    `quiverMutationAtVertices` is the same walk without the commentary.  This
    also flags a step where the Coxeter polynomial moves, which it must not do
    along an admissible path -- if it does, the sequence left the admissible
    region (research R-005).

    It was called `quiverMutation`, which now names the package.
    """
    baseCoxPol = invariants.coxeterPoly(pathAlg)
    print(invariants.coxeterPoly(pathAlg))
    pathAlgebra.printPathAlgebra(pathAlg)
    if firstDisplayedStep == 0:
        plotting.plotQuiver(pathAlg)
    for i in range(len(mutationVertexList)):
        if mutationVertexList[i] >= 0:
            pathAlg = quiverMutationAtVertex(pathAlg, mutationVertexList[i])
        else:
            pathAlg = leftQuiverMutationAtVertex(pathAlg, -mutationVertexList[i])
        pathAlg = reduction.reducePathAlgebra(pathAlg)
        currentCoxPol = invariants.coxeterPoly(pathAlg)
        cartMat = invariants.cartanMatrix(pathAlg)
        print(np.matrix(cartMat))
        if currentCoxPol != baseCoxPol:
            print('COXETER POLYNOMIAL HAS CHANGED!')
        print('Mutations: ', mutationVertexList[0:i + 1])
        print(currentCoxPol)
        pathAlgebra.printPathAlgebra(pathAlg)
        if i + 1 >= firstDisplayedStep:
            plotting.plotQuiver(pathAlg)
    return pathAlg


def reverseMutationSequence(mutationVertices, vertexNumbering):
    reverseMutationVertices = []
    for n in range(len(mutationVertices) - 1, -1, -1):
        reverseMutationVertices.append(-getVertexNumberingKeyFromValue(vertexNumbering, mutationVertices[n]))
    return reverseMutationVertices


def getVertexNumberingKeyFromValue(vertexNumbering, vertex):
    keyList = list(vertexNumbering.keys())
    valueList = list(vertexNumbering.values())
    if vertex > 0:
        position = valueList.index(vertex)
        key = keyList[position]
    else:
        position = valueList.index(-vertex)
        key = -keyList[position]
    return key
def leftQuiverMutationAtVertex(pathAlg, vertex):
    dualPathAlg = pathAlgebra.dualPathAlgebra(pathAlg)
    mutDualPathAlg = quiverMutationAtVertex(dualPathAlg, vertex)
    leftMutPathAlg = pathAlgebra.dualPathAlgebra(mutDualPathAlg)
    return leftMutPathAlg
