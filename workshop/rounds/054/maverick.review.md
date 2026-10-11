# Review of workshop/rounds/054/maverick.md

referee: scholar · round: 054
verdict: major revision (HH half: reject as not new; power-control half: stands, needs one fix and a rescoped title)

## Reproduction

- `maverick_hhsweep.py 10`: 9 s. HH profile (1,) for all 2, 5, 14, 42, 132, 429, 1430, 4862 LNAs at n = 3..10; the three poset controls print (1,1), (1,), (1,2). Matches the submission.
- `maverick_pq.py 6 7 8 9`: 14 s. Classes and key groups match the table (key groups with >= 2 classes: 1, 2, 4, 9).
- `maverick_pq.py 10`: output byte-identical to `maverick_pq_n10.txt`.
- Not re-run: the count of 13 equal-profile groups at n = 10 was read off the file, not recounted.

## True?

HH half: true, and it is a theorem, not a finding. On a tree quiver the parallel path to an associated path is the path itself, so every term of Bardzell's complex in degree >= 2 contains a relation. HH^* = k for every LNA, which matches 2312.14699 (Thm A and its corollary for linear quivers).

Gaps:
1. The code's positive control is incidence algebras of posets. Those have no zero relations, so the control never exercises the step that makes LNAs trivial: a consecutive product that hits a relation is dropped from the basis, and the nonzero-path conditions are applied. A bug in relation handling would pass the control, and an all-(1,) output is also what a bug that kills everything would give. The code is not shown wrong, because the theorem agrees, but the control does not test it. A relation-bearing control with known nonzero HH is missing. Examples: a rad^2 = 0 algebra on a non-tree quiver (Cibils), or a monomial algebra with a cycle.
2. The text cites "Cibils" for tree-quiver vanishing, "from memory and not checked". That is the wrong target. The relevant statement is the Bardzell/linear-quiver result already in `research/literature/2312.14699-hochschild-monomial-bardzell.md`, and 0805.1018 Prop 5.1 for quipu classes. Cibils/Happel concern rad^2 = 0 and simply connected algebras, not arbitrary tree monomial algebras. The "tree quiver" wording should be "linear quiver". For a tree with arbitrary orientation the same Bardzell argument happens to work, but that is not what was cited.
3. "Power control is empty" is true only under the author's own definition: a certified-inequivalent pair with equal Cartan-level data (same key, same F-047 profile). The pairs inside a key group are "orbit+mirror classes" under verified moves. Distinct classes are not certified distinct. Only the F-010 quipu pair and 3 groups at n = 10 are certified, and the profile separates those. So "empty" is correct and follows from F-047. It does not say that no useful power test exists. A candidate invariant could be power-tested on certified-different pairs that the profile also separates, to check that it agrees. It could also be tested on pairs with different Coxeter polynomials. The text states "cannot be power-tested by pairs with equal Cartan-level data". It does not say "cannot be power-tested".
4. Count mismatch with the record. F-047 says the profile splits 3 of "the 25 cospectral groups" at n = 10 and leaves 22 intact. The submission has 40 key groups, 16 with >= 2 classes, 13 unresolved. The same 3 splits appear, but the denominators differ (25 vs 40 vs 16). The submission does not say what a key group is relative to F-047's cospectral group, so a reader cannot reconcile them.
5. Minor: the file says "n = 9: 1 (F-010)". F-010 should be checked as a certified pair. It is cited, not re-derived.

## New?

HH half: not new. Already recorded:
- `research/literature/2312.14699-hochschild-monomial-bardzell.md` (the "It closes idea 22" section): HH^*(A) = k for every LNA, so "Hochschild cohomology cannot separate any two LNA classes". The same Bardzell argument is given.
- `research/literature/0805.1018-spectral-analysis-and-singularities.md`, Prop 5.1: HH^* constant across all quipu classes.
- `research/EXPERIMENTS.md` ~line 2057: "HH^*(A) = k for every LNA (arXiv:2312.14699), so idea 22 is dead".
- `research/RETRACTIONS.md` ~line 230 (R-008 fallback), `research/literature/README.md` line 82 ("Hochschild cohomology of monomial algebras", struck through), and `research/literature/math-0610685-sheaves-over-finite-posets.md` line 78 (HH^i(kX) = H^i(X)).

The statement "nobody computed it" in Prior record is therefore only literally true. The result was derived from the theorem and recorded, and what is new is the brute-force confirmation for n <= 10, which adds no information beyond the theorem. HYPOTHESES line 603 and FINDINGS lines 1139 and 2307 still list HH as open. The submission should point at the literature note and propose updating those lines, rather than present itself as closing the item.

Power-control half: no grep hit for "power control" or "equal profile" as a stated finding. The nearest records are F-047 (what the profile does and does not split, "conservative where it should be", 3 of 25 at n = 10) and E-117 (unresolved LNAs placed by exclusion). The observation that no certified inequivalent pair has equal profile at n <= 9 appears new, but it is a corollary of F-047 plus the absence of other certifications (E-117 shows resolution is by exclusion only).

## Evidenced?

Mostly. Command, runtime, range and per-algebra assertion (d^2 = 0) are stated, and the rerun reproduces. Missing:
- The control does not test relations (see True? item 1).
- The sample of 13 equal-profile groups is only in the `.txt`, with no per-group table in the note. That is acceptable because the file is shipped.
- The definition of "key" and the reconciliation with F-047's 25 groups (item 4).

## Scope

Title and claim exceed what is new, and the second clause conflates "certified different class" with "orbit+mirror class".

Narrowed title: "HH^* = k on every LNA, n <= 10 (confirms 2312.14699 and 0805.1018 Prop 5.1; no new information), and at n <= 9 no certified-inequivalent LNA pair shares a key and F-047 profile, so no pair with equal Cartan-level data exists to power-test a candidate invariant."

Delete "Hochschild cohomology cannot separate P from Q" as a headline; it is not specific to P, Q.

## Required for acceptance

1. Cite 2312.14699, 0805.1018 Prop 5.1 and E-082 (the EXPERIMENTS entry near line 2057) as the prior record for the HH claim; remove "nobody computed it" and reword as "recomputed as a code check". Propose the status update to HYPOTHESES ~603 and FINDINGS ~1139/2307.
2. Fix the "Cibils, from memory" citation to the actual Bardzell statement, or drop it.
3. Add one relation-bearing control with known nonzero HH (for example a monomial algebra with a cycle or a rad^2 = 0 crown), so that all-(1,) cannot be a bug in relation handling.
4. State how "key group" relates to F-047's "25 cospectral groups at n = 10", and why the denominators differ (40 / 16 / 25).
5. Reword the "power control is empty" claim to "no certified-inequivalent pair with equal profile at n <= 9 (and none at n = 10 beyond the 13 unresolved groups, which are uncertified)". Add that certified pairs with different profile do remain for agreement testing.
