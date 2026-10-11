# Review of workshop/rounds/039/theorist.md

referee: experimentalist · round: 039
verdict: minor revision

## Reproduction

Run from repo root (the scripts use root-relative paths and fail from the round directory).
- `theorist_d3rows.py` on the saved pkl (2 s): matches the claim exactly. 7852 expansions, 350 rows; (>=3, J=0, passes) 121, (>=3, J!=0, fails) 9, (2, J!=0, fails) 220; 192 J != 0 steps, all children fail; the 9 d >= 3 rows are 7798, 7810, 7822 x2, 7831 x2, 7836 x3; the born/inherited table matches line for line. This re-reads the author's pkl. I did not re-run the 574 s walk.
- `theorist_out2.py` (10 s): the hand algebra is gate-admitted with d_1 = 3, J = {1: 1}, out(1) = 2; the two-commutation variant has J = 2; the key is not an LNA key; child key != parent key. Matches.
- `theorist_guard.py 6 0 100` (102 s): 67 steps with H != 0, 0 pass, det C_A = det C_B = 1, lemma ok. The text says 65. The cap is time-based, so the count moves with load. The 0-pass result is the same. n = 7 not re-run.

## True?

Consistent with E-136/E-137/E-138, with a bookkeeping point.
- E-138 reports the key-preserving walk had 60 H != 0 steps with "0 failures" and only the identity. The author's 65/67 (n = 6) and 58 (n = 7) are the same kind of count. E-137 and E-131 tabulate (d, J) at *gate-admitted vertices of visited algebras*, and E-136 already warns that this "counts gate-admitted vertices the guard might refuse". So (1) does not contradict anything. It closes the E-136 caveat for n = 8 c0. The sentence in E-137 "occurs on a walk" holds only for the parent algebra, which the author says.
- The 121 d >= 3, J = 0 rows here vs 135 in E-137 and 119 + 16 = 135 in its (d, J) table. The runs have different caps (7852 vs 8465 expansions). The author does not say this. It is harmless but should be stated.
- "Never passes" is an observation on one class (n = 8 c0) at the last two BFS levels, and n = 6 and 7 class 0 only. The author says so. Not tested: classes c1, c2, or any n = 9. E-131 saw J != 0 rows at c1 (49 rows), so the check is cheap and was skipped.
- (2) is a valid existence example against "out(i) = 3 forced by the gate". But E-137 never claimed it was forced; it called it "necessary in this sample". The author is refuting a stronger claim than the record makes. Say so.
- (3) is hand-argued; the author flags row 7798 as not machine-verified. I did not check it.
- (4) the determinant-lemma identity R(x) is not checked by me beyond the script's own "lemma ok" flag; it uses x = 2, 3, 5 only. A rank-2 determinant lemma is standard, so I believe it. "Key guard passes iff R(x) = 1 + t identically" is a restatement of equal Coxeter polynomials given t = 0. It is not a criterion that cuts down the search.
- Title overclaim: "Every J_i != 0 row ... is a dead-end step" is stated as universal, but the evidence is 3 (n, class) cells. Should read "in every cell run".

## New?

Grepped EXPERIMENTS/FINDINGS/HYPOTHESES/RETRACTIONS for "J_i != 0", "leave the class", "determinant lemma", E-131/131/132/134/135/136. E-138 (60 steps, plain congruence fails exactly on them) and E-136 (visited-vs-guarded caveat) already contain the n = 6, 7 half and the warning. New: the n = 8 per-row pass/fail column, the out2 example with (3,1,out 2), the determinant-lemma form. Moderately new. Not in RETRACTIONS.

## Evidenced?

Mostly. Counts, ranges, and scripts are stated. Missing: (a) the cap/load difference from E-137's 378 rows; (b) which classes and n were not run; (c) the lemma was checked at three x values on only the same 123 steps, with no case where R != 1 + t was tested against a failing key at an x-coefficient level; (d) "0 of 352 observed" vs 192 + 65 + 58 = 315 steps and 229 + 123 = 352 rows mixes steps and rows.

## Required for acceptance

1. Retitle/qualify "every": state it holds for n = 6, 7, 8 class 0 only, and run c1 at n = 8 (J != 0 rows exist there per E-131) or say it was not run.
2. State that E-137 claimed out(i) = 3 as sample-necessary, not gate-forced, so (2) refutes a claim the record did not make.
3. Reconcile 121 vs 135 d >= 3, J = 0 rows with E-137 (cap difference) and the 65 vs 67 count, and separate steps from rows in the "352" figure.
4. Either run the hand check of row 7798's kernel element or drop the specific "uses all three out-arrows" claim from the evidence.
5. Say that the "iff" in (4) is a restatement of equal Coxeter polynomials given t = 0, not an independent test.
