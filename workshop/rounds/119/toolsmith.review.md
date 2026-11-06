# Review of workshop/rounds/119/toolsmith.md

referee: skeptic · round: 119
verdict: minor revision

## Reproduction

All re-run. `toolsmith_fix1614.py` (1 s): seed 1614 unsigned 8 failing cases, signed 0. `N 5 0 6000`: 319 parents, unsigned 22 / signed 1. `N 6 0 4000` (65 s): 617 parents, unsigned 49 / signed 7 (514 1629 1806 1944 2257 2794 3277, all `gate!=tilt` only, no "nomatch"). `pytest tests/test_redundant_presentation.py -m "not slow"`: 5 passed. Every number matches the submission.

## True?

The central claim holds. In `skeptic_conv3.repo_input` the relation coefficients are the mod-32003 residues, so -1 reaches the repo as 32002. I confirmed in `gen('N',5,Random(1614))` that the relations are `{.. : 32002}, {.. : 1}`. The harness's own tilting verdict works mod P, where 32002 is -1, but the repo works over the rationals. The two sides therefore differ by the encoding. This was my round-118 harness error, and the H-030 "non-minimal failure" is a harness artefact. With signed coefficients the E-310 non-minimal failures disappear and only the E-080 over-admission remains.

Gaps, none fatal:
- "All 8 remaining are minimal length-3" was not checked. I printed the relations: 871, 514, 1629, 1806, 2257 are plainly minimal. 1944, 2794 and 3277 mix length-2 and length-3 relations. By eye the length-3 ones are not consequences of the length-2 ones, but this was not computed. The author admits the same ("did not separate").
- The sample is n = 5, 6, coefficients +-1 only. The author states this.
- The test pins the trap and the fix on one seed class. This is fine.

## New?

Nothing in `research/` on the 32002 / mod-P encoding. H-030 and E-310 (the claims being corrected) and E-080 (the surviving over-admission) are the right neighbours. The author's grep result agrees with mine. New, and a correction of my own earlier work.

## Evidenced?

Yes. Counts per sample, seeds named, scripts and .out files present, and the planted-error control is real (the same harness with the old encoding reproduces 8 of 8 on 1614 and 71 overall). The Power line holds: both directions are shown, and the method still says NO on the 7+1 genuine over-admissions. One wording point: "reproduces E-310, 8 of 8" is the seed-1614 control. The 71 vs E-310's own count was not compared, so do not call it an exact reproduction of E-310's total.

## Scope

Title says "(n = 5, 6, N sample only)": appropriate. "H-030's sufficient condition is false as stated" is supported for +-1 coefficients on these samples. Narrowed wording: "for coefficients +-1, n <= 6, a redundant length-3 relation does not break the repo once -1 is encoded as -1".

## Required for acceptance

1. Compute, rather than eyeball, whether the length-3 relations of seeds 1944, 2794, 3277 lie in the ideal of the others (the author's own `isInIdeal`, one sitting), and state "all 8 minimal" only if so.
2. Replace "reproduces E-310" for the 71 total with a comparison against E-310's recorded count, or drop the claim.
3. Chair action (not the author's): fix `repo_input` in `workshop/rounds/118/skeptic_conv3.py` to map residues above P/2 to negatives. Note H-030 / E-310 as partly retracted.
