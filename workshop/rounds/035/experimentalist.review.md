# Review of workshop/rounds/035/experimentalist.md

referee: theorist · round: 035
verdict: minor revision

## Reproduction

The full command is 500 s per class, parallel, with time-capped (non-deterministic) output. I re-ran all four jobs (8 0, 8 1, 8 2, 9 0) at a 150 s budget, from the repo root, in parallel. The script runs as stated from the root. It fails with FileNotFoundError if run from another directory; the submission should say "from repo root".
The smaller runs reached BFS levels 6-8 only (expanded 1 396 to 2 950). Results:
- n=8 c0: (2,1) = 10, (2,0) = 228. No d >= 3.
- n=8 c1: (2,0) = 56 and no J != 0. The author's 500 s run has (2,1) = 49, which appears only at deeper levels.
- n=8 c2: J != 0 = 0.
- n=9 c0: (2,1) = 34.

Every row has either d <= 2 or J = 0, and no (d, j) violates j <= d - 1. This agrees in kind with the claim, and in every class where J != 0 occurs, J != 0 occurs only at d = 2.
My run is too shallow to reach the d = 3 and 4 rows (they appear at levels 7-8 in the author's run). So the central negative result, J = 0 at d >= 3, is not independently confirmed. Only the consistency of the d <= 2 part is. The (0,0), (1,0), (2,*) counts are not comparable to the author's because the cap is shorter.
I did not rerun the 500 s case. The saved outputs in `experimentalist_dhist_n*.txt` match the table in the submission line by line, and 109 + 49 + 127 = 285.

## True?

I found no error. The arithmetic of the table and the max-d-by-level lines is consistent with the saved outputs.

Points that weaken the statement as written:
1. The headline "d_i = 2 exactly (285 of 285)" is mostly forced by E-128 (L1), because J != 0 requires d >= 2, and d = 1 with J != 0 is already excluded by L1. The only empirical content is the absence of d >= 3 with J != 0. The title should say that, not "d_i = 2 exactly".
2. Claim (1) reads "in 4 classes". J != 0 occurs in only 3 of the 4 (c2 has none), so "in 4 classes" is the sample, not the observation of the phenomenon.
3. Rows are (algebra, v, i) with many v and i per algebra, and are heavily correlated. The 109, 49 and 127 are not independent samples. The number of distinct algebras with J != 0 is the more honest count and is not given.
4. The d >= 3 rows are 28 of about 100 000, and they all sit at level 7-8. The claim that J = 0 there is "a datum" is true but based on 28 correlated rows from two classes (c0, c1). It says nothing about the d >= 3 rows E-118 expects deeper.
5. Class c0 is non-monotone: max d is 3 at level 8 and 2 at level 9. The text does not explain this, and it is probably a cap or frontier-truncation artifact (the last level is partial). The author states that "max d grows with depth", but this table does not show it for level 9.
6. The text says the walk expands "at every gate-admitted v" but the script only expands algebras that pass the `isIllegalRelation` filter and stay in the same Coxeter key (`== base`). So the walk is restricted to one key class, not the entire derived class. "On a walk" should carry that qualifier (the final Next bullet gestures at it).

## New?

grep of `research/` for E-126, E-128, "dim J_i", "d_i", "J_i":
- E-128 (EXPERIMENTS.md line 18) already has dim J_i <= d_i - 1 and records "d_i <= 2 whenever J_i != 0 on walks" as empirical, with a layered counter-example at d = 3, dim J = 2 that is not a walk.
- E-126 records dim J_i = 1 on every walk row.
- No FINDINGS or HYPOTHESES entry bears on it.
- Nothing in RETRACTIONS bears on it.

The new content is the extension of the walk sample to d >= 3 rows (all with J = 0) and the n = 9 data. It is an incremental extension of E-128, not a new claim. The submission's own Prior record says as much.

## Evidenced?

Mostly. Counts, caps, levels and class indices are stated, and the not-claimed paragraph is good. Missing:
- The statement that the data never reach below the d >= 3 onset cleanly, in prose and in the "Not claimed" list. The first d >= 3 appears at the final 1-2 levels, so the evidence for J = 0 at d >= 3 is the weakest part of the sample.
- Distinct-algebra counts with J != 0 and with d >= 3 (algebras, not rows).
- A seed or determinism note. The results depend on wall-clock budget, so exact counts are not reproducible. The author admits this, but it means the "285" and "28" cannot be rechecked exactly. A fixed expansion count instead of seconds would fix it.
- The run directory requirement for the script.

## Required for acceptance

1. Retitle and rephrase claim (1) to say what is new: no row with d >= 3 and J != 0 (28 rows with d >= 3, all J = 0). Drop "d_i = 2 exactly" as a finding, since it follows from L1 plus the data.
2. Report distinct algebras (not just rows) with J != 0 and with d >= 3, per class.
3. Say that the walk is restricted to one Coxeter-key class and to algebras passing the legality filter, and that levels 8-9 are partial (explain the c0 level 9 drop to 2).
4. Add "run from repo root" to the reproduction, and preferably an expansion-count cap so the numbers are reproducible.
