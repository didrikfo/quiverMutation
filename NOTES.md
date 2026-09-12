# Working notes

State of the codebase, what the model can and cannot express, and the backlog.
Written while picking the project back up; keep it current.

## What the code does

The repo implements the combinatorial rule for tilting mutation of
[arXiv:2112.08129](https://arxiv.org/abs/2112.08129), and applies it to classify
linearly oriented Nakayama algebras (LNAs) up to derived equivalence, which is
how the data behind [arXiv:2305.06642](https://arxiv.org/abs/2305.06642) was
produced.

### The model of a path algebra

`pathAlgebraClass.PathAlgebra` holds

* `quiver`: a `networkx.MultiDiGraph` whose nodes are the vertices (integers,
  `1..n` for a line);
* `rels`: a list of relations.

A **relation** is a list of paths, and a **path** is a list of vertices. So
`[[1,2,4],[1,3,4]]` is the commutativity relation between the two paths from 1
to 4, and `[[1,2,3]]` is the zero relation on the single path `1 -> 2 -> 3`.
By convention `rel[0][0]` is the source and `rel[0][-1]` the target, and every
path in a relation shares both.

This is the deliberate simplification: a relation is a **set** of paths, not a
linear combination. There are no coefficients, so the model cannot tell
`p - q = 0` from `p + q = 0`, and cannot express `2p - 3q + r = 0` at all. A
relation with one path means that path is zero; a relation with two or more
means they are identified up to sign. See "Improving the model" below.

A vertex has an implied identity path, which nothing in `rels` represents; the
`+1` on the diagonal of `cartanMatrix` stands in for it.

### The mutation procedure

`quiverMutationAtVertex(pathAlg, vertex)` applies steps 1-7 of the paper's
procedure and returns a new `PathAlgebra`. It does not clean up after itself:
the result routinely contains inadmissible relations (ones containing a path of
length one) and redundant relations.

`reducePathAlgebra(pathAlg)` does the cleanup, in this order: drop illegal
relations, dedupe relation paths, cancel each length-one path in a relation
against the corresponding arrow (the "Note" after the paper's step 7),
substitute through the arrows that cancelled, drop non-minimal zero relations,
drop redundant relations and existing subrelations, iterate to a fixed point.

`quiverMutationAtVertices(pathAlg, vertices)` mutates at each vertex in turn,
reducing after each. A negative entry means left mutation, which
`leftQuiverMutationAtVertex` performs as dual -> right mutation -> dual.

`mutationIsPossibleAtVertex(pathAlg, vertex)` is the admissibility test: an
arrow out of the vertex must exist, the quiver must have no parallel arrows, and
`Hom(P_i*[1], Lambda)` must vanish.

### The LNA search

`mutationSearch(lineLength, mutationDepthStart, startRow, createNewCSVfile)` is
the entry point.

1. `generateAllPossibleLineRelations(n)` enumerates every admissible ideal on the
   linear quiver `1 -> ... -> n`. There are Catalan(n-1) of them.
2. `createMutationClassCSV(n)` writes one row per LNA to
   `A_<n>_mutation_classes.csv`, with empty class columns.
3. For each row still without a class, `mutationSearchDepthFirst` walks
   mutations to the given depth and appends every LNA it reaches to a text file.
4. `mutationListLineCleanup` normalises each of those back to the standard
   numbering `1..n` and dedupes.
5. `saveLineRelationsAndMutationsToCSV` writes the class name, the mutation path
   from the class representative, the Coxeter polynomial and the vertex
   numbering into every row it reached.

A class is named after the LNA that seeded it, in the per-vertex relation-length
notation: `A7_22300` is the line on 7 vertices with relations of 2, 2 and 3
arrows starting at vertices 1, 2 and 3. The same notation appears as the
`relationList` argument of `lineQuiverExample`.

The depth-first search is bounded by `depth`, which `mutationSearch` decays as
`max(depthStart - floor(log10(row)), 2)` to keep later rows affordable, and it
stops descending as soon as a cycle appears in the quiver.

Because the depth is bounded, **the search finds a lower bound on each class**:
two LNAs in the same class can end up in different search classes if no mutation
path between them fits in the depth. Merging those was the manual step. The
Coxeter polynomial identifies the candidates, since it is invariant under
derived equivalence — but it is not a complete invariant, so agreement is not
proof. It does happen to separate every class for n <= 8.

## Verified against the papers

* `generateAllPossibleLineRelations(n)` returns Catalan(n-1) relation sets for
  n = 2..11.
* The acyclic worked example of arXiv:2112.08129 (mutate the 7-vertex quiver at
  vertex 3) reproduces the paper's quiver and relations exactly.
* The Coxeter polynomial does not move along any legal mutation path of depth
  <= 3 out of 14 different LNAs of length 5 to 7.
* Every class in the published n <= 8 table has a single Coxeter polynomial, and
  distinct classes have distinct ones.
* n = 5: the search finds both classes outright, 8 + 6 = 14 LNAs.
* n = 6: the search finds 5 classes with 4 distinct Coxeter polynomials. Merging
  the two that share the D_6 polynomial gives the paper's 4 classes for n = 6
  (A_6, D_6, E_6, D~_5) with sizes 16 + 13 + 12 + 1 = 42.
* n = 7: the search finds 11 classes with 6 distinct Coxeter polynomials,
  matching the 6 quipus of order 7 in the paper, sizes
  32 + 29 + 54 + 7 + 6 + 4 = 132.

## Known gaps and limitations

* **Step 3's cyclic case is not implemented.** For a minimal relation
  `r: i --> i` the procedure calls for one arrow `(alpha r-bar): i* -> t(alpha)`
  per arrow `alpha` out of `i`; the code adds a single arrow from `r`'s source to
  its target, which on a cycle is a loop on `i*`. Pinned as a strict xfail in
  `tests/test_mutation_procedure.py`. The LNA search never reaches it because it
  stops descending at the first cycle, so this has never affected a published
  result.
* **`mutationIsPossibleAtVertex` is stricter than the paper** when the vertex has
  more than one arrow out of it: it rejects the vertex as soon as one arrow out
  of it kills a nonzero path, where the paper only asks that some arrow out of it
  does not. The two agree when the vertex has a single arrow out of it, which is
  the only case the LNA search meets.
* **Parallel arrows are rejected outright**, because a path is a vertex sequence
  and so cannot name which of two parallel arrows it uses.
* **The recursive loop at length 12.** A full n = 12 run crashed in a recursive
  loop, somewhere in the mutation of relations. The offending LNA is not
  recorded. Reproducing it is the first step; a long-running n = 12 search with
  the recursion limit lowered and the failing quiver dumped on `RecursionError`
  would find it.
* **Relations are compared by sorted vertex sequences**, so two relations that
  are equal as ideals but written differently are distinct objects. Several
  functions exist mainly to paper over this (`removeDuplicateRels`,
  `removeDuplicateRelPaths`, `removeRedundantRelations`,
  `removeExistingSubrelations`, `minimizeCommutingRelation`).
* **Results round-trip through text files**, parsed by string slicing at fixed
  offsets in `readMutationsFromFile`. Fragile and slow; the search should hand
  the mutation list to the CSV writer in memory.

## Backlog

Roughly in the order that unblocks the most.

### Performance and correctness

1. **Memoize the path enumeration.** With the graph deepcopy gone, the hot spot
   is `nx.all_simple_paths`, called 93k times in one depth-4 search on quivers
   that have not changed between calls. Needs a cheap structural key for a path
   algebra — which is also what a proper `__hash__`/`__eq__` on `PathAlgebra`
   would give.
2. **Find and fix the length-12 recursive loop**, per above.
3. **Canonical form for a path algebra.** A normalised, hashable representation
   would replace the dedupe-by-list-comparison machinery, let the search
   memoize on quivers it has already visited, and make `reducePathAlgebra`
   testable as "reduces to the canonical form".
4. **Avoid re-deriving relations from scratch after each mutation.** The search
   recomputes every relation between every pair of vertices at each node; most
   of that is unchanged by a mutation at a single vertex.

### The model

5. **Linear combinations in relations.** Give a relation coefficients over the
   base field (or over Z, or over {+1,-1} as a first step) so it can express
   `sum c_p p = 0`. This is the main expressiveness limit, and it is what makes
   the steps 5 and 7 bookkeeping so delicate. Worth scoping as: what breaks if a
   relation becomes `dict[tuple[int, ...], int]` instead of `list[list[int]]`?
6. **Name arrows.** Paths as vertex sequences cannot express parallel arrows or
   distinguish two arrows with the same endpoints. Arrow identities would lift
   that restriction and make the quiver a plain `DiGraph` of named arrows.

### OOP and structure

7. **Make `PathAlgebra` carry its own behaviour.** It is currently a thin
   container and every operation is a free function taking it as the first
   argument. Mutation, reduction, invariants and admissibility are all methods.
8. **Subclass for LNAs**, holding `n` and the per-vertex relation lengths, with
   the standard numbering, the relation-string form, the Kupisch series and the
   relation dual as its own operations. Most of the `...Line...` free functions
   collapse into it.
9. **Subclass for quipu quivers**, holding the `P^{(m)}_{(k)}` parameters, the
   orientation, and the CR-swap of arXiv:2305.06642 section "Cord/relation-swap"
   as a method. This is what makes idea 12 possible.
10. **Split the 2600-line module** along the seams that already exist: the
    procedure, relation algebra, invariants, the line search, quipus, IO.
11. **Replace the text-file round trip** with in-memory objects, and write one
    CSV (or parquet) at the end.

### Features

12. **Better CSV post-processing, to cut the manual work.** Currently: group the
    rows by Coxeter polynomial, and for each group of two or more search classes
    try to find a mutation path between their representatives with a deeper or
    targeted search; report what merged, what did not, and what still shares a
    polynomial with nothing found. This is the highest-value feature — it is the
    step that was done by hand for n <= 11.
13. **Seed the table from quipu quivers.** For large n, use theorem
    `thm:QuipuToAn` of arXiv:2305.06642 to write down the class of every LNA with
    almost separate relations directly from the quipus of order n, and let the
    search fill in only the rest. `generateAllQuipus`, `count_quipus` and
    `generateAllHeightOneQuipus` already exist.
14. **More invariants, to separate classes the Coxeter polynomial cannot.**
    Candidates: the determinant and elementary divisors of the Cartan matrix,
    the Euler form, the number of indecomposables / the shape of the AR quiver,
    Hochschild cohomology dimensions, the derived invariants of Avella-Alaminos
    and Geiss for gentle algebras (LNAs are gentle), and the silting/tilting
    quiver's local structure. Each wants to be a function
    `PathAlgebra -> hashable`, so the merge step can key on a tuple of them.
15. **Certificates both ways.** A merge should record the mutation path that
    proves the equivalence; a split should record the invariant that separates
    the two classes. Then a classification is checkable without rerunning it.
16. **Push past n = 11.** Catalan growth means n = 12 is 58786 LNAs and n = 15
    is 2674440, so the search has to get cheaper per LNA and the table has to
    stop being a CSV read into memory.
