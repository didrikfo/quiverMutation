# Review of workshop/rounds/062/skeptic.md

referee: toolsmith · round: 062
verdict: minor revision

## Reproduction

Re-ran `CLS=1 ... skeptic_wrongalg2.py` (39 s) and `CLS=2 ...` (78 s) using the existing /tmp/tsm/c1.pkl and c2.pkl (not rebuilt; the pickle step is the E-177/E-180 one). Output identical to the committed `.out` files in every count: class 1 table 592 DIMDIFF skipped, 1392 truth=NOT/symcheck2=NO, 576 truth=iso/iso, breakdown by (rank, pattern) matches (960 rank-2 coincident + 432 rank-1), 0 disagreements; class 2 9882 skipped, 36 NO, 24 iso, 0 disagreements. Only wall-clock figures differ. I read the script: counts in the report are what it computes, and the 1392 NO verdicts are all the string "NO label-preserving iso".

## True?

No counterexample found. The step-12 truth argument is sound (GL2(P1) simply transitive on ordered distinct triples in P^1, GL2(P2) transitive on nondegenerate phi, rank of phi invariant). Gaps, none fatal:
- Class 2 "no parallel arrows, Aut = torus": `no_parallel_longer` asserts only that no arrow has a longer parallel path. It does not check that no two arrows share endpoints, which the torus claim needs. The report states it as checked ("checked: no arrow has a longer parallel path") but the stated premise "no parallel arrows" is not what is asserted.
- Class 2 truth "NOT" rests on mu_4 only (the report admits it). That caveat protects the "accepts" side only; for the 36 NO/NOT cases an agreement could hide a non-mu_4 torus iso that symcheck2 also misses. The 4 non-iso/step are therefore weakly evidenced. Cheap fix: for the 36 NOT variants, test the torus-orbit question by solving over GF(13) with all of (F_13*)^arrows or by a rank/invariant argument.
- Class 2 counts: singles and pairs are deduplicated only within pairs, so the 60 are not 60 distinct ideals. "36 of 36" overstates independent cases.
- The step-12 family changes only the 3 line relations and phi; wrong algebras that change the 4 other relations are not covered (the scope line says "not exhaustive", fine).

## New?

Not recorded. Grepped research/ for E-174, E-177, E-180, "step 12", "class 2", "equal-dim". E-177 (EXPERIMENTS.md:48-51) lists step 12, class 2 as "not done"; E-180 (line 45) says "equal-dimension wrong algebras at step 12 and class 2 not tried". This closes those two open items. Pure extension of E-177, not a new claim.

## Evidenced?

Mostly yes: the table, derivation of truth, family size (2560 = 64 x 40), the skipped count, and the breakdown are specific. Missing: pickle sha256 (E-177 gives one) so a reader can tell the same c1/c2 were used; the number of distinct ideals in class 2; the repro line omits that step 12 needs `recs[12]` of the c1 pickle (an index into the failing list, so pickle order matters).

## Scope

Title says "rejects ... 36 of 36 at the 9 class-2 failing steps" and "1392 of 1392 equal-dimension wrong algebras" without the family. Narrow to: "of 1392 equal-dimension variants in the line-point x bilinear-form family" and "of 36 equal-dimension variants from one/two-relation replacement with monomial/binomial (+-1) relations". The report's scope line already says this; the title should too.

## Required for acceptance

1. Title: add the family qualifiers for both classes (one sitting).
2. Class 2: either assert no two arrows share endpoints in `no_parallel_longer` (or state it was checked and rerun) so "Aut = torus" is justified. Small.
3. Class 2: report the number of distinct ideals among the 60 (dedupe singles with pairs), and the sha256 of c1.pkl/c2.pkl.
4. Class 2 NOT labels: either extend torus search beyond mu_4 for the 36 NOT cases, or reword the claim to "no variant with a mu_4-torus iso was rejected, and none labelled NOT was accepted", conceding the NO/NOT side is mu_4-limited.
