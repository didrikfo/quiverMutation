# Review of workshop/rounds/038/experimentalist.md

referee: theorist · round: 038
verdict: minor revision

## Reproduction

Ran the command from the repo root (it opens `workshop/rounds/033/experimentalist_bothdie.py` by relative path, so it fails from any other cwd). 560 s cap, 9m25 wall. Mine stopped at 8 503 expansions, 8 481 algs (author: 8 465 / 8 443), 383 rows (author 378). (d, J) counts: identical for every d >= 3 cell: (3,0) 119, (4,0) 16, (3,1) 2, (4,1) 4, (5,1) 3. Only (2,1) differs: 239 vs 234 (the extra 38 expansions). The 9 d >= 3, J != 0 rows, 5 algebras, exp 7798, 7810, 7822 x2, 7831 x2, 7836 x3, all out(i) = 3: same. d = 2 out(i) split: 2:194, 1:35, 3:10 (author 191/35/8), consistent with extra rows. d >= 3, J = 0 out(i) split identical to the author's. The core claim reproduces; the d = 2 numbers are load-dependent, as the author says.

## True?

Statement as written is true of the sample. Gaps:
- "out(i) = 3 necessary" rests on 9 rows / 5 algebras, all at BFS level 9, all in the last ~670 expansions, and rows from 3 expansion-indices cluster (7822, 7831, 7836 are one or two algebras each). Effective independent n is about 5, with dim A 47/64/75 each a different (d, dimA) pairing. The text says this; the headline title still says "necessary in the sample", which is accurate but weak.
- Algebra count of 5 is by id (key or hash); the two with keys vs three with hashes means pairwise distinctness is unchecked (author says so).
- Prediction of the claim: d >= 3 with J != 0 and out(i) != 3 does not appear. Not tested past the cap; the J = 0, d >= 3 rows show out(i) from 1 to 7 on the same tail, so the out(i) = 3 restriction is not a tail artefact of out-degree range generally, but the author does not compare J != 0 d = 2 vs d >= 3 rows at equal level (all d >= 3 are level 9; 146 d = 2 J != 0 rows at level 9 for comparison would be the fair control, and the author's own out(i) distribution there, which I find as 2/1/3, has 3 only about 5 percent).
- Mechanism is not offered; E-131 already gives one (relation p1 b = p2 b on the single non-parallel out-arrow). The two new rows have three distinct targets and no parallel arrow, so that mechanism does not cover them: this is the interesting part and is correctly flagged.

## New?

E-129 (d = 2 on all 285 rows; d >= 3 only at J = 0), E-131 (three d = 4, 4, 5 rows, J != 0, filtered on out-degree >= 3 of v; not out(i)), E-132 (hand (3,1) example, off-walk). Nothing found in research/ for "out-degree(i)" or a (3,1) row on a walk. So the two (3,1) walk rows with distinct targets, and the out(i) = 3 tabulation, are new. The d = 2 count (234) restates E-129's kind of count at larger scale; not new in content.

## Evidenced?

Mostly. Counts, table, script names, pkl and the cap condition are stated. Missing:
- The d >= 3 rows should be listed in the review text (aid, dimA, v, i, out targets, Cartan row) so the claim can be checked without opening the pkl; only the aggregate table is given.
- The rerun command works only from the repo root and depends on round 033's script; say so.
- "8 of the d = 2 out-degree-3 rows": count moves with load (8 vs 10 mine); the title should not rely on it. The title's "9 rows / 5 algebras" is stable.
- "Not sufficient" is supported (d = 2 rows with out(i) = 3), but that is a statement about d = 2 rows; the sufficiency question that matters is: among rows with out(i) = 3 and J != 0, what fraction has d >= 3 (9 of 17 to 19 in my run). State that ratio.

## Required for acceptance

1. State the cwd requirement and the round-033 dependency in Reproduction.
2. Replace load-dependent d = 2 counts in the title/claim by "about" or give both runs (mine: 239 (2,1), 10 with out(i) = 3).
3. Add a short per-row listing of the 9 d >= 3, J != 0 rows.
4. Give the conditional ratio P(d >= 3 | J != 0, out(i) = 3) (9 of 17 here by author's counts, 9 of 19 by mine) rather than only "not sufficient".
5. Keep "Not claimed" for non-isomorphism; either check it (Cartan matrix plus quiver isomorphism) or leave "5" as "up to 5 ids".
