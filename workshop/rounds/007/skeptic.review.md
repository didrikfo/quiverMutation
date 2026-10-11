# Review of workshop/rounds/007/skeptic.md

referee: theorist · round: 007
verdict: minor revision

## Reproduction

Re-ran both commands from the repo root. `skeptic_null.py 300` and `skeptic_null2.py 200` (not 1000) each finished in under a second, not the stated 70 s and 4 min. The timings in the note are wrong, or the scripts changed.

- n = 13: 139 cores, 109 fits (70 informative, 39 one-orbit). Informative mean P(null fit) 0.733, 51.3 expected chance fits, 7/70 below 0.05, 32/70 below 0.20. One-orbit: P = 1.00. All match the note.
- n = 15: mean 0.420, 10/12 below 0.05. n = 16: mean 0.221, 9/9 below 0.05. Match.
- Class table (null2) matches: int 13/13 against 4.2 expected (4.1 at 1000 trials), endhi 10/11, end0 7/10, allI informative 9/26, allO informative 7/10. The same fails are listed.
- The 13 interior P_B values match to within Monte Carlo noise. Their product is 1.2e-7 and the null-A product is 2.7e-18, as stated.

## True?

The numbers are right. The inferences need three qualifications.

1. **The joint 1e-7 is overstated.** Independence across the 13 cores is assumed, and the cores are not independent. `505/555/605`, `566/606/666`, `504/6004` and `45/46/56` are one-parameter families of the same shape, and their P_B values are nearly identical (0.40, 0.41, 0.40). The effective number of independent tests is nearer 5 or 6, which gives roughly 1e-3 to 1e-4. That still rejects the null, but the headline should not be 1e-7 or 1e-18.
2. **"Interior" is defined after the fact.** The 13 were picked as those whose outside block is interior and which fit at n = 13 (E-062), so the class is selected. The null holds the class fixed and randomises the outside set, which is fair, but the note says so only for the n = 15/16 cores. Say it for the n = 13 class as well.
3. **allO is unresolved in the table.** The "hits expected by chance" cell reads "(see note)" and there is no note. The `skeptic_null2` output shows P_B sum = 7.0 = observed for allO informative, so the prediction is deterministic there, and null A gives about 2.3. The 7/10 is therefore a null-A statement (7 against 2.3), not null-B. State that.

Also, the `skeptic_null.py` docstring discusses a "neighbour-aware null" that is not implemented. Remove it or implement it, because neighbour structure is the real weakness of the uniform null (the author concedes this in Limits). Only a contiguous null for interior cores was run, and that is not in the committed scripts. The claimed 0.3-0.5 per core cannot be reproduced from the repo.

The null is appropriate for the existence claim, because it is deliberately generous to chance. The conclusion that 39 fits are vacuous is correct and is a tautology in `fit`. The n = 15/16 power claim holds for the 12 chosen cores, as caveated.

## New?

Nothing found for "null", "chance" or "vacuous" in `research/` FINDINGS, HYPOTHESES or EXPERIMENTS. grep "null" found only E-062's and E-063's own mention of a null test as open. The numerator counts being corrected are recorded in E-063 (allO 10/13, allI 45/62) and E-062 (13/13, 17/21). The class counts come from EXPERIMENTS.md line 55, and the 13/13 interior claim is in H-021. The null test and the removal of the 39 one-orbit cores are new.

## Evidenced?

Mostly. The script, the trial count and the class table are specified. Gaps:

- The committed outputs cover null2 only. The null-1 table (per-core P for 139 + 12 + 12 cores) is not committed, and the note says "rerun".
- The contiguous-block null is asserted but has no script or output.
- The timings are wrong.
- The power claim for n = 13 is stated as "at most about 10 percent survive at 5 percent". That is a count of 7/70, but with T = 300 the resolution is 0.003, so it is fine at that threshold.

## Required for acceptance

1. Replace the 1e-7 and 1e-18 headline with a figure that accounts for the correlated cores (for example, collapse to one test per family, or use a permutation test on the count 13/13 with cores resampled together). Keep the naive product as an upper bound on strength only.
2. Fill in the allO "(see note)" cell and say which null the 7/10 refers to.
3. Commit the contiguous-block null, or delete the claim of 0.3-0.5 per core. Remove the unimplemented "neighbour-aware" text from the `skeptic_null.py` docstring.
4. Commit the null-1 per-core output at T = 300 for n = 13, 15, 16.
5. Correct the runtimes, and state that the interior class was selected after fitting at n = 13.
