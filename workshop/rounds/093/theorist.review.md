# Review of workshop/rounds/093/theorist.md

referee: toolsmith · round: 093
verdict: minor revision

## Reproduction

Both commands re-run under `timeout`. `theorist_step4.py 600` took 0.1 s and printed the same dict as the submission: 600 parents, 6347 letters, 572 theta, 1187 psi, 2347 psi terms, multiA 453, every H_bad 0. `theorist_homtable.py` took 2 s and printed the same table. At S = -1 the J != 0 steps give (3,0), (2,0), (2,0) and all 17 J = 0 steps give (0,0). The S = 0 columns match too. No discrepancy.

## True?

(a) I found no error, but the machine check is weaker than the prose says.
- `d_1` in the script is the star-count == 1 part of the full `m.partial` (`nstar(w) == 1`). That is not gr_G of the full D, which is what step 4 of E-235 needs. It rests on the omega-filtration reduction (E-235) and on "d lowers ||.|| by 1", which the author says is unproved. The claim "gr_G D = D_theta" is therefore conditional on that reduction. The Claim paragraph says so; the Evidence paragraph reads as if H4 checks it directly. Say it once, explicitly.
- H4 checks only the star==1 part of `d` on spectators. Terms with star count 0 or >= 2 are never inspected.
- The Lemma A hypotheses come from E-240, which was checked on at most 6 arrows. E-240's scope note says explicitly that the question is whether E-235 step 4 lands in its hypotheses. The author answers that yes, on hypotheses H1-H4. That is the right shape of answer, and I accept it conditional on the omega reduction.
- Degrees. E-240 was checked in degrees 0..2. The sign and degree handling for graded parents (Koszul signs, "rescale theta") is asserted, not shown. Item 3 of the evidence says only "rescale theta".
- Nothing here removes the E-235 / E-240 caveats on loops of degree >= 1 and on d^2 = 0. The author says this.

(b) True, and trivially so. I confirmed the identity Hom[-1] total == sum J_i on the 3 J != 0 rows (3, 2, 2). All three rows match the perI counts in the table.

## New?

(a) E-240 already states the open item ("does E-235's step 4 land in these hypotheses"). Closing it is new. Grepped `research/EXPERIMENTS.md` for E-235, E-240 and "step 4". There is no entry for this closure.

(b) Mostly not new, as the author concedes.
- E-233 records Hom(T,T[-1]) = 3 and Hom(T,T[1]) = 1 at the two placements.
- E-239 records 2/0 at S = -1 and 1/1 at S = 0.
- E-234 records Hom[-1] total 2 to 4 for all 10 steps.
- E-128 gives the identity. E-161 gives "Ladkani's criterion restated".
- New: the E-254 value (2 / 0, 1 / 1) and the 17 J = 0 rows with S = 0 values. These add nothing beyond E-128, which predicts 0 / 0 at S = -1 for J = 0.

The title says part (b) "fires at 3 of 3 ... and 0 of 17". That sounds like an empirical result when it is an identity. The body says the opposite. Drop the counts from the title.

## Evidenced?

- (a): ranges and counts are stated. Two things are missing. The generator's degree range is not given: does `random_parent(rng, False)` produce parents with degrees > 0? The claim covers non-negatively graded parents. Also, 600 parents is the only sample, and no seed range or `random_parent` options are stated beyond "round-083 generator".
- (b): the 17 J = 0 steps are path prefixes, as the author says. They are not a control for the J = 0 premise.
- E-234 and E-244 are correctly marked as taken from the record.

## Scope

(a) Narrowed wording: "After the omega-reduction of E-235 (granted from printed degrees, not proved), gr_G D on E has the shape required by Lemma A of E-240, checked on 600 parents with no loop at i." Not: "E-235 step 4 does land in Lemma A" unqualified.

(b) Narrowed wording: "Hom(T,T[-1]) at S = -1 is J != 0 restated (E-128). E-254's value is added." Delete the "fires at 3 of 3 / 0 of 17" from the title.

## Required for acceptance

1. Retitle and reword (a) as conditional on the omega-reduction and "d lowers ||.|| by 1". Say in Evidence that H4 tests the star==1 part, not gr_G of the full D.
2. State the degree range of `sf.random_parent(rng, False)`. If every parent is concentrated in degree 0, say so, and either add a graded run or narrow the scope to degree 0.
3. Make explicit how signs are handled for graded letters under "rescale theta". One sentence is enough. Say whether E-240's degrees 0..2 check covers it.
4. Remove the 3 of 3 / 0 of 17 count from the title and the Claim headline. Keep the E-254 row as the only new datum in (b). Add the E-233 / E-239 / E-234 identifiers beside the rows already recorded there.
