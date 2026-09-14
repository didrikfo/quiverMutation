"""Tools for mutating quivers with relations and working with their path algebras.

The mutation procedure is the combinatorial rule for tilting mutation of
D. Fosse, *A combinatorial procedure for tilting mutation*, arXiv:2112.08129,
and the main application is the classification of linearly oriented Nakayama
algebras up to derived equivalence, which produced the results in *Quipu quivers
and Nakayama algebras with almost separate relations*, arXiv:2305.06642.

The modules, in dependency order:

| module | holds |
|---|---|
| `pathAlgebra` | the container: a quiver and a list of relations |
| `paths` | paths and relations inside one, and what a relation set does to a path |
| `reduction` | the cleanup after a mutation, to a fixed point |
| `mutation` | the procedure itself, and its admissibility condition |
| `invariants` | the Cartan matrix and the Coxeter polynomial |
| `lines` | the linear quiver, and the names its algebras go by |
| `search` | walking the mutation graph, and the hereditary quivers it reaches |
| `classification` | classifying a whole length, end to end |
| `quipuForms` | canonical forms for quipus, and the quipu theorem inverted |
| `relationAlgebra` | relations as integer combinations of paths, and exact ideals |
| `mutationClassTable` | the classification table, as CSV and parquet |
| `nakayama` | `LinearNakayamaAlgebra` and `QuipuAlgebra`, the two shapes with structure |
| `lnaMoves` | verified mutation shortcuts between LNAs |
| `piecewiseHereditary` | certificates that an algebra is in no quipu class |
| `quiverExamples` | small quivers to try things on by hand |
| `plotting` | drawing a quiver |

The procedural modules are re-exported flat, so `import quivermutation as qm`
reaches `qm.classifyLength` and the rest directly.  The modules that carry their
own namespace -- `nakayama`, `quipuForms`, `lnaMoves`, `piecewiseHereditary`,
`relationAlgebra`, `mutationClassTable`, `quiverExamples` -- are imported as
names: `from quivermutation import nakayama as nk`.
"""

from . import (
    pathAlgebra,
    paths,
    reduction,
    mutation,
    invariants,
    lines,
    search,
    classification,
    quipuForms,
    relationAlgebra,
    mutationClassTable,
    nakayama,
    lnaMoves,
    piecewiseHereditary,
    quiverExamples,
    plotting,
)

from .pathAlgebra import (
    PathAlgebra,
    dualPathAlgebra,
    printPathAlgebra,
)

from .paths import (
    allMinimalRelsBetweenVertices,
    allRelsBetweenVertices,
    allRelsInPathAlgebra,
    applyRelSetToPath,
    extendRel,
    isIllegalRelation,
    isSubRelOf,
    listIntersection,
    numberOfPathsUpToRels,
    pathHasZeroRel,
    powerset,
    sublistExists,
    zeroizeRels,
)

from .reduction import (
    reducePathAlgebra,
    removeDuplicateRelPaths,
    removeDuplicateRels,
    removeExistingSubrelations,
    removeNonminimalZeroRels,
    removeRedundantRelations,
)

from .mutation import (
    getVertexNumberingKeyFromValue,
    leftQuiverMutationAtVertex,
    mutationIsPossibleAtVertex,
    quiverMutationAtVertex,
    quiverMutationAtVertices,
    reverseMutationSequence,
    showMutationSteps,
)

from .invariants import (
    cartanMatrix,
    coxeterPoly,
)

from .lines import (
    className,
    generateAllPossibleLineRelations,
    mutationListLineCleanup,
    relSetToString,
    relabelLineAlgebra,
    relationStringToLineRelLengths,
)

from .search import (
    findHereditaryFormForClass,
    formatHereditaryForms,
    hereditaryFormFromTheorem,
    hereditaryFormsReachedFrom,
    mutationSearchDepthFirst,
)

from .classification import (
    adoptClassesByMoves,
    annotateHereditaryForms,
    assignMutationClassInTable,
    classifyLength,
    expandClassByMoves,
    mergeReport,
    mutationSearch,
    resolveMergeCandidates,
    seedTableFromQuipuTheorem,
)

from .plotting import (
    plotQuiver,
)

__all__ = [
    "PathAlgebra",
    "adoptClassesByMoves",
    "allMinimalRelsBetweenVertices",
    "allRelsBetweenVertices",
    "allRelsInPathAlgebra",
    "annotateHereditaryForms",
    "applyRelSetToPath",
    "assignMutationClassInTable",
    "cartanMatrix",
    "className",
    "classification",
    "classifyLength",
    "coxeterPoly",
    "dualPathAlgebra",
    "expandClassByMoves",
    "extendRel",
    "findHereditaryFormForClass",
    "formatHereditaryForms",
    "generateAllPossibleLineRelations",
    "getVertexNumberingKeyFromValue",
    "hereditaryFormFromTheorem",
    "hereditaryFormsReachedFrom",
    "invariants",
    "isIllegalRelation",
    "isSubRelOf",
    "leftQuiverMutationAtVertex",
    "lines",
    "listIntersection",
    "lnaMoves",
    "mergeReport",
    "mutation",
    "mutationClassTable",
    "mutationIsPossibleAtVertex",
    "mutationListLineCleanup",
    "mutationSearch",
    "mutationSearchDepthFirst",
    "nakayama",
    "numberOfPathsUpToRels",
    "pathAlgebra",
    "pathHasZeroRel",
    "paths",
    "piecewiseHereditary",
    "plotQuiver",
    "plotting",
    "powerset",
    "printPathAlgebra",
    "quipuForms",
    "quiverExamples",
    "quiverMutationAtVertex",
    "quiverMutationAtVertices",
    "reducePathAlgebra",
    "reduction",
    "relSetToString",
    "relabelLineAlgebra",
    "relationAlgebra",
    "relationStringToLineRelLengths",
    "removeDuplicateRelPaths",
    "removeDuplicateRels",
    "removeExistingSubrelations",
    "removeNonminimalZeroRels",
    "removeRedundantRelations",
    "resolveMergeCandidates",
    "reverseMutationSequence",
    "search",
    "seedTableFromQuipuTheorem",
    "showMutationSteps",
    "sublistExists",
    "zeroizeRels",
]
