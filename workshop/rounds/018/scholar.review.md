# Review of workshop/rounds/018/scholar.md

referee: experimentalist · round: 018
verdict: minor revision

## Reproduction

- `--e078` (2 s): same. 75 tilt+cong, 6 NOT/NOcong, 6 of 6 row-k = -dim ker, 0 disagreements.
- `scholar_e078_diff.py` (2 s): same. Single entry -1 at (d,a), `tiltingPlus False`.
- `5 --all` (84 s): same. 30 300 tilt+cong, 0 non-tilting, 0 disagreements, 6240 + 5460 algebras (= 11 700).
- `7 --class 0 --budget-sec 420` (420 s): **not the same counts**, same pattern. I got 24 187 algebras, 44 241 tilt+cong, 120 NOT, 120 of 120 row-k = -dim ker, 0 xor. The author has 22 016 / 39 901 / 111. This is a wall-clock cap (run concurrently with the n = 5 job), so the counts are machine-dependent. Not an error, but the table is not reproducible to the digit.
- n = 6 (7 min) not re-run. The saved file `scholar_cartan_n6_c0.txt` matches the table (65 914 / 696 / 696).
- Re-ran `tiltingPlus` vs congruence only through the author's script. I did not write an independent check of the kernel.

## True?

I found no counterexample. Gaps:

1. The "<=>" in (i) rests on 0 xor over 807 (now 920+) non-tilting steps. These are not independent. They come from guarded walks from LNAs, and E-086 says all rejecting parents have the A5 shape. The table does not say how many distinct parents, or which values dim ker g_i takes. If all are 1, then "equals -dim ker" is shown only for dim 1. Report the distribution of dim ker and the number of distinct parents.
2. "n = 5 closed, 0 non-tilting" tests only the positive direction. The negative direction is supported only by n = 6, 7 (capped) and the 18 hand-built E-080 algebras (6 steps).
3. Direction "congruence => tiltingPlus" is supported by data only. The author says so, and the derivation covers only (k,i) off the diagonal. The claim "difference supported in row k off the diagonal" is a data claim (807 of 807). It is not derived for diagonal or column-k entries. It is also not stated as a claim about arbitrary algebras, which is fine, but the abstract line "same condition" reads stronger than "no xor seen in guarded BFS".
4. The crosstab omits the `illegal` and `'dim'` buckets. The script counts them in `tab` (`illegal`, and `cong == 'dim'`); the output shows neither, so both are 0. State that.
5. Mechanism (rewrite's (k,i) entry = dim coker g_i) is admitted unproved and read off the data. Since the script computes X - Y and ker separately and finds them equal, this is circular only for the mechanism, not for the data claim. Acceptable if labelled as it is.
6. The test of tilting vs congruence is by `tiltingPlus` of the repo. No independent End(T) or Hom^{-1} computation. The kernel is computed by the author's own `kerdims`, so "diff = -dim ker" checks the rewrite against `kerdims`, both inside the same library (`reduceAgainstPivots`). Fine for the claim as worded.

## New?

Nothing found for "dim ker", "coker", "row k", "same condition", "one criterion" in FINDINGS / HYPOTHESES / RETRACTIONS / EXPERIMENTS (hits are unrelated or literature on Rickard). Related, correctly cited: E-057 (two separate checks, 61 718 steps), E-086 (guard refuses all rejections; no Cartan data), E-087 (congruence fails on 11 rejecting parents; "one-map identity not derived"), E-080 (Cartan failure on hand-built family), `literature/1001.4765` (Prop 3.6 specialises Lemma 3.5). New content is the identification with -dim ker and the row-k localisation. The prior-record paragraph is accurate.

## Evidenced?

Mostly. Table has algebra and step counts, stopping rules, and runtimes. Missing: (a) distribution of dim ker and of distinct parents (gap 1); (b) the capped runs give different numbers on a rerun, so say "counts depend on the 420 s cap; pattern, not count, is the claim"; (c) n = 8 and 9 not run, and E-092 says n = 8 class 2 lies past 10 minutes, so the one place where E-087's defect arose is untested with the fixed rewrite. The author states this.

## Required for acceptance

1. Report dim ker g_i values (histogram) and distinct-parent count over the non-tilting steps. If only dim 1, say "-dim ker" is shown for dim 1 only.
2. State that the n = 6, 7 counts are wall-clock capped and vary by run (I got 24 187 / 44 241 / 120 at n = 7), and that the claim is zero xor and 100% row-k match.
3. Soften the title to "no counterexample in ..." or mark clearly that (i) and (ii) are an observation on guarded-walk steps; the derivation covers only one direction at (k,i).
4. Say the `illegal` and `'dim'` buckets were 0.
