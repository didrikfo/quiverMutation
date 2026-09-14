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
| `lines` | the linear quiver, and the three names its algebras go by |
| `search` | walking the mutation graph, and the hereditary quivers it reaches |
| `classification` | classifying a whole length, end to end |
| `quipuForms` | canonical forms for quipus, and the quipu theorem inverted |
| `relationAlgebra` | relations as integer combinations of paths, and exact ideals |
| `mutationClassTable` | the classification table, as CSV and parquet |
| `nakayama` | `LinearNakayamaAlgebra` and `QuipuAlgebra`, the two shapes with structure |
| `lnaMoves` | verified mutation shortcuts between LNAs |
| `piecewiseHereditary` | certificates that an algebra is in no quipu class |
| `fileio` | the CSV tables and the older text transcripts |
| `plotting` | drawing a quiver |
| `legacy`, `legacyQuipus` | superseded; nothing in the package depends on them |

Everything the flat `quiverMutation` module used to expose is re-exported here,
so `import quivermutation as qm` reaches it all; new code is better off importing
the module it wants.
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
    fileio,
    plotting,
    quipuForms,
    relationAlgebra,
    mutationClassTable,
    nakayama,
    lnaMoves,
    piecewiseHereditary,
    quiverExamples,
    legacy,
    legacyQuipus,
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
    applyCommutativityRelSetToPath,
    applyRelSetToPath,
    extendRel,
    isIllegalRelation,
    isSubRelOf,
    listIntersection,
    nonMinimalOutRels,
    numberOfPathsUpToRels,
    pathHasZeroRel,
    powerset,
    replaceSubPath,
    sublistExists,
    zeroizeRels,
)

from .reduction import (
    minimizeCommutingRelation,
    reduceCommutativityRels,
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
    onePointExtension,
    quiverMutation,
    quiverMutationAtVertex,
    quiverMutationAtVertices,
    reverseMutationFromSequence,
    reverseMutationSequence,
)

from .invariants import (
    cartanMatrix,
    cartanMatrixForCanonicalAlgebra,
    coxPolyOfTree,
    coxeterPoly,
    coxeterPolyForCanonicalAlgebra,
    divisors,
    generateAllCoxeterPolynomials,
)

from .lines import (
    convertLineFromCSVnotation,
    generateAllKupischSeries,
    generateAllLineQuiversWithRelations,
    generateAllPossibleLineRelations,
    isRelationDualLineQuiver,
    lineQuiverExample,
    lineRelLengthsToClassName,
    makeStandardLineQuiver,
    mutationListLineCleanup,
    mutationListLineCleanupKeepDupes,
    relSetToString,
    relabelLineAlgebra,
    relationDualLineQuiver,
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
    saveLineRelationsAndMutationsToCSV,
    seedTableFromQuipuTheorem,
)

from .fileio import (
    createMutationClassCSV,
    generateListOfRelations,
    importMutationClassCSV,
    readMutationClassesFromFile,
    readMutationsFromFile,
    readRelationsFromFile,
    saveLinePathAlgMutation,
    saveLineRelationsAndMutationsToFile,
    saveLineRelationsToFile,
    saveQuipusToCSV,
)

from .plotting import (
    plotQuiver,
)

from .legacy import (
    collectMutationClasses,
    combineLineMutationFiles,
    combineMutationClasses,
    combineMutationClassesInCSVfile,
    expandAllClassesWithEasyRels,
    expandClassFurtherWithEqualRelPairs,
    expandClassWith2Rels,
    findMutationClassesForLine,
)

from .legacyQuipus import (
    bfs_shortest_path_to_subgraph,
    bfs_shortest_path_to_subgraph_edges,
    bfs_shortest_path_to_subgraph_path,
    count_quipus,
    count_quipusV1,
    dfs_shortest_path,
    generateAllHeightOneQuipus,
    generateAllQuipus,
    generateAllQuipusGPT,
    generateAllQuipusUpToLength,
    generate_quipus,
)

__all__ = [
    "PathAlgebra",
    "adoptClassesByMoves",
    "allMinimalRelsBetweenVertices",
    "allRelsBetweenVertices",
    "allRelsInPathAlgebra",
    "annotateHereditaryForms",
    "applyCommutativityRelSetToPath",
    "applyRelSetToPath",
    "assignMutationClassInTable",
    "bfs_shortest_path_to_subgraph",
    "bfs_shortest_path_to_subgraph_edges",
    "bfs_shortest_path_to_subgraph_path",
    "cartanMatrix",
    "cartanMatrixForCanonicalAlgebra",
    "classification",
    "classifyLength",
    "collectMutationClasses",
    "combineLineMutationFiles",
    "combineMutationClasses",
    "combineMutationClassesInCSVfile",
    "convertLineFromCSVnotation",
    "count_quipus",
    "count_quipusV1",
    "coxPolyOfTree",
    "coxeterPoly",
    "coxeterPolyForCanonicalAlgebra",
    "createMutationClassCSV",
    "dfs_shortest_path",
    "divisors",
    "dualPathAlgebra",
    "expandAllClassesWithEasyRels",
    "expandClassByMoves",
    "expandClassFurtherWithEqualRelPairs",
    "expandClassWith2Rels",
    "extendRel",
    "fileio",
    "findHereditaryFormForClass",
    "findMutationClassesForLine",
    "formatHereditaryForms",
    "generateAllCoxeterPolynomials",
    "generateAllHeightOneQuipus",
    "generateAllKupischSeries",
    "generateAllLineQuiversWithRelations",
    "generateAllPossibleLineRelations",
    "generateAllQuipus",
    "generateAllQuipusGPT",
    "generateAllQuipusUpToLength",
    "generateListOfRelations",
    "generate_quipus",
    "getVertexNumberingKeyFromValue",
    "hereditaryFormFromTheorem",
    "hereditaryFormsReachedFrom",
    "importMutationClassCSV",
    "invariants",
    "isIllegalRelation",
    "isRelationDualLineQuiver",
    "isSubRelOf",
    "leftQuiverMutationAtVertex",
    "legacy",
    "legacyQuipus",
    "lineQuiverExample",
    "lineRelLengthsToClassName",
    "lines",
    "listIntersection",
    "lnaMoves",
    "makeStandardLineQuiver",
    "mergeReport",
    "minimizeCommutingRelation",
    "mutation",
    "mutationClassTable",
    "mutationIsPossibleAtVertex",
    "mutationListLineCleanup",
    "mutationListLineCleanupKeepDupes",
    "mutationSearch",
    "mutationSearchDepthFirst",
    "nakayama",
    "nonMinimalOutRels",
    "numberOfPathsUpToRels",
    "onePointExtension",
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
    "quiverMutation",
    "quiverMutationAtVertex",
    "quiverMutationAtVertices",
    "readMutationClassesFromFile",
    "readMutationsFromFile",
    "readRelationsFromFile",
    "reduceCommutativityRels",
    "reducePathAlgebra",
    "reduction",
    "relSetToString",
    "relabelLineAlgebra",
    "relationAlgebra",
    "relationDualLineQuiver",
    "relationStringToLineRelLengths",
    "removeDuplicateRelPaths",
    "removeDuplicateRels",
    "removeExistingSubrelations",
    "removeNonminimalZeroRels",
    "removeRedundantRelations",
    "replaceSubPath",
    "resolveMergeCandidates",
    "reverseMutationFromSequence",
    "reverseMutationSequence",
    "saveLinePathAlgMutation",
    "saveLineRelationsAndMutationsToCSV",
    "saveLineRelationsAndMutationsToFile",
    "saveLineRelationsToFile",
    "saveQuipusToCSV",
    "search",
    "seedTableFromQuipuTheorem",
    "sublistExists",
    "zeroizeRels",
]
