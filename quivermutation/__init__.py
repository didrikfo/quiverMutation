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
| `relationAlgebra` | relations as integer combinations of paths, and exact ideals |
| `procedure` | steps 1-7 and the cleanup, on those combinations |
| `reduction` | the cleanup, in the set-of-paths model |
| `mutation` | the procedure, in the set-of-paths model, and admissibility |
| `invariants` | the Cartan matrix and the Coxeter polynomial |
| `lines` | the linear quiver, and the names its algebras go by |
| `search` | walking the mutation graph, and the hereditary quivers it reaches |
| `classification` | classifying a whole length, end to end |
| `quipuForms` | canonical forms for quipus, and the quipu theorem inverted |
| `mutationClassTable` | the classification table, as CSV and parquet |
| `classview` | reading a finished classification back, one row per class |
| `classpage` | the same, rendered as a page to browse |
| `nakayama` | `LinearNakayamaAlgebra` and `QuipuAlgebra`, the two shapes with structure |
| `lnaMoves` | verified mutation shortcuts between LNAs, floating and anchored |
| `endMoves` | the listed half of those: the rewrites that need an end of the quiver |
| `overlap` | how much an LNA's relations overlap, and how far the rules reach |
| `piecewiseHereditary` | certificates that an algebra is in no quipu class |
| `quiverExamples` | small quivers to try things on by hand |
| `plotting` | drawing a quiver |

The procedural modules are re-exported flat, so `import quivermutation as qm`
reaches `qm.classifyLength` and the rest directly.  The modules that carry their
own namespace -- `nakayama`, `quipuForms`, `lnaMoves`, `overlap`,
`piecewiseHereditary`,
`relationAlgebra`, `procedure`, `mutationClassTable`, `classview`, `classpage`,
`quiverExamples` -- are imported as names:
`from quivermutation import nakayama as nk`.
"""

from . import (
    pathAlgebra,
    paths,
    relationAlgebra,
    procedure,
    reduction,
    mutation,
    invariants,
    lines,
    search,
    classification,
    quipuForms,
    mutationClassTable,
    nakayama,
    lnaMoves,
    endMoves,
    overlap,
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
    nameClassesFromTheorem,
    nameRemainingClasses,
    resolveMergeCandidates,
    seedTableFromQuipuTheorem,
)

from .plotting import (
    plotQuiver,
)

__all__ = [
    "PathAlgebra",
    "adoptClassesByMoves",
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
    "leftQuiverMutationAtVertex",
    "lines",
    "lnaMoves",
    "mergeReport",
    "nameClassesFromTheorem",
    "nameRemainingClasses",
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
    "procedure",
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
    "resolveMergeCandidates",
    "reverseMutationSequence",
    "search",
    "seedTableFromQuipuTheorem",
    "showMutationSteps",
    "sublistExists",
]
