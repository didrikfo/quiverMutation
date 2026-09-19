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
| `relationAlgebra` | relations as integer combinations of vertex paths, and exact ideals |
| `arrowPaths` | the same over paths that name their arrows, so parallel arrows can be said |
| `procedure` | steps 1-7 and the cleanup, on those combinations |
| `reduction` | the cleanup, in the set-of-paths model |
| `mutation` | the procedure, in the set-of-paths model, and admissibility |
| `invariants` | the Cartan matrix and the Coxeter polynomial |
| `coxeterTables` | the polynomials of a whole length, as a table to match against |
| `lines` | the linear quiver, and the names its algebras go by |
| `search` | walking the mutation graph, and the hereditary quivers it reaches |
| `reflections` | reorienting a relation-free tree, which costs mutations and nothing else |
| `classification` | classifying a whole length, end to end |
| `quipuForms` | canonical forms for quipus, and the quipu theorem inverted |
| `treeSearch` | every tree as a hereditary algebra, against every LNA |
| `quipuRelations` | quipus that do carry relations, against every LNA |
| `mutationClassTable` | the classification table, as CSV and parquet |
| `classview` | reading a finished classification back, one row per class |
| `classpage` | the same, rendered as a page to browse |
| `nakayama` | `LinearNakayamaAlgebra` and `QuipuAlgebra`, the two shapes with structure |
| `lnaMoves` | verified mutation shortcuts between LNAs, floating and anchored |
| `endMoves` | the listed half of those: the rewrites that need an end of the quiver |
| `spectatorMoves` | the rewrites that hold with a relation in the window they do not touch |
| `overlap` | how much an LNA's relations overlap, and how far the rules reach |
| `edgeMoves` | the doubling at an end, which the window encoding cannot state |
| `freeMoves` | relations of two arrows, which cost no mutation at all |
| `piecewiseHereditary` | certificates that an algebra is in no quipu class |
| `quiverExamples` | small quivers to try things on by hand |
| `plotting` | drawing a quiver |

The procedural modules are re-exported flat, so `import quivermutation as qm`
reaches `qm.classifyLength` and the rest directly.  The modules that carry their
own namespace -- `nakayama`, `quipuForms`, `lnaMoves`, `overlap`, `edgeMoves`, `freeMoves`,
`piecewiseHereditary`, `coxeterTables`, `treeSearch`, `quipuRelations`,
`reflections`,
`relationAlgebra`, `arrowPaths`, `procedure`, `mutationClassTable`, `classview`,
`classpage`, `quiverExamples` -- are imported as names:
`from quivermutation import nakayama as nk`.
"""

from . import (
    pathAlgebra,
    paths,
    relationAlgebra,
    arrowPaths,
    procedure,
    reduction,
    mutation,
    invariants,
    coxeterTables,
    lines,
    search,
    reflections,
    classification,
    quipuForms,
    treeSearch,
    quipuRelations,
    mutationClassTable,
    nakayama,
    lnaMoves,
    endMoves,
    spectatorMoves,
    overlap,
    edgeMoves,
    freeMoves,
    piecewiseHereditary,
    classview,
    classpage,
    quiverExamples,
    plotting,
)

from .pathAlgebra import (
    PathAlgebra,
    dualPathAlgebra,
    printPathAlgebra,
)

from .paths import (
    allRelsBetweenVertices,
    allRelsInPathAlgebra,
    applyRelSetToPath,
    extendRel,
    isIllegalRelation,
    numberOfPathsUpToRels,
    pathHasZeroRel,
    powerset,
    sublistExists,
)

from .reduction import (
    reducePathAlgebra,
)

from .mutation import (
    getVertexNumberingKeyFromValue,
    leftQuiverMutationAtVertex,
    leftQuiverMutationAtVertex,
    mutationIsPossibleAtVertex,
    quiverMutationAtVertex,
    quiverMutationAtVertices,
    reverseMutationSequence,
    showMutationSteps,
)

from .invariants import (
    cartanMatrix,
    coxeterCoefficients,
    coxeterKey,
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
    DEEPER_CONDITIONS,
    DeeperWhen,
    deeperWhenFromSpec,
    describeNode,
    describeRelationFreeQuiver,
    findHereditaryFormForClass,
    formatHereditaryForms,
    hereditaryFormFromTheorem,
    hereditaryFormsReachedFrom,
    hasNoRelations,
    hasOrientedCycle,
    hasParallelArrows,
    linesReachedFrom,
    meetingPoints,
    mutationSearchDepthFirst,
    parallelArrowsAtLeast,
    quiverKey,
    quiversReachedFrom,
    recordOnly,
    relationFreeSightings,
    summariseSightings,
)

from .classification import (
    adoptClassesByMoves,
    annotateHereditaryForms,
    assignMutationClassInTable,
    classifyLength,
    expandClassByMoves,
    mergeReport,
    mutationSearch,
    nameClassesFromTheorem,
    nameRemainingClasses,
    progressPathFor,
    readProgress,
    resolveMergeCandidates,
    seedTableFromQuipuTheorem,
    writeProgress,
)

from .plotting import (
    plotQuiver,
)

__all__ = [
    "PathAlgebra",
    "adoptClassesByMoves",
    "arrowPaths",
    "allRelsBetweenVertices",
    "allRelsInPathAlgebra",
    "annotateHereditaryForms",
    "applyRelSetToPath",
    "assignMutationClassInTable",
    "cartanMatrix",
    "className",
    "classification",
    "classifyLength",
    "classpage",
    "classview",
    "coxeterCoefficients",
    "coxeterKey",
    "coxeterPoly",
    "coxeterTables",
    "DEEPER_CONDITIONS",
    "DeeperWhen",
    "deeperWhenFromSpec",
    "describeNode",
    "describeRelationFreeQuiver",
    "dualPathAlgebra",
    "edgeMoves",
    "expandClassByMoves",
    "extendRel",
    "findHereditaryFormForClass",
    "formatHereditaryForms",
    "freeMoves",
    "generateAllPossibleLineRelations",
    "getVertexNumberingKeyFromValue",
    "hasNoRelations",
    "hasOrientedCycle",
    "hasParallelArrows",
    "hereditaryFormFromTheorem",
    "hereditaryFormsReachedFrom",
    "invariants",
    "isIllegalRelation",
    "leftQuiverMutationAtVertex",
    "lines",
    "linesReachedFrom",
    "lnaMoves",
    "meetingPoints",
    "mergeReport",
    "mutation",
    "mutationClassTable",
    "mutationIsPossibleAtVertex",
    "mutationListLineCleanup",
    "mutationSearch",
    "mutationSearchDepthFirst",
    "nakayama",
    "nameClassesFromTheorem",
    "nameRemainingClasses",
    "numberOfPathsUpToRels",
    "parallelArrowsAtLeast",
    "pathAlgebra",
    "pathHasZeroRel",
    "paths",
    "piecewiseHereditary",
    "plotQuiver",
    "plotting",
    "powerset",
    "printPathAlgebra",
    "procedure",
    "progressPathFor",
    "quipuForms",
    "quipuRelations",
    "quiverExamples",
    "quiverKey",
    "quiverMutationAtVertex",
    "quiverMutationAtVertices",
    "quiversReachedFrom",
    "readProgress",
    "recordOnly",
    "reducePathAlgebra",
    "reduction",
    "reflections",
    "relSetToString",
    "relabelLineAlgebra",
    "relationAlgebra",
    "relationFreeSightings",
    "relationStringToLineRelLengths",
    "resolveMergeCandidates",
    "reverseMutationSequence",
    "search",
    "seedTableFromQuipuTheorem",
    "showMutationSteps",
    "sublistExists",
    "summariseSightings",
    "treeSearch",
    "writeProgress",
]
