# Review of workshop/rounds/133/skeptic.md

referee: theorist · round: 133
verdict: minor revision

## Reproduction

- `units 10` (15 s): 140 units, 18 self-dual, classes [1, 2, 3, 5, 5, 58, 66]. Matches.
- `n9` (0.4 s): 6 distinct units, pairs as stated, class-1 unit means 2.915/2.857/2.375/2.917/2.683, SD 0.230, class 2 = 2.711, offset -0.17 SD. Matches. By hand: the class-1 mean is 2.749 and the stdev (n-1) is 0.230, so (2.711 - 2.749)/0.230 = -0.17. Correct.
- `power 10 1` (2 min 1 s) and `power 11 1` (4 min 34 s): `diff` against `skeptic_power10.txt` and `skeptic_power11.txt` is empty. Identical.
- Duality unit at n = 10, which the paper asserts but does not test: I walked 30330030, the opposite of probe 30033030, with the script's `collect` writing to my scratchpad (2 min 28 s). It gave the same 120 members and the same histogram, 1:16 2:32 3:45 4:16 5:11. That histogram also matches the paper's probe table. The unit holds on this pair.

## True?

1. **(3) "cannot be met at n <= 11 even with every unit walked" is too strong.** I ran the registered statistic on all 1361 n = 11 units (20 classes; features maxC and nnzC; 99 permutations) at a planted 0.17 SD. Power was 0.75 (150 reps) and 0.80 (40 reps). The size check at delta 0 gave 0.06, and delta 0.12 gave 0.46. So a full n = 11 census sits right at the 0.8 line. It is not out of reach. The paper's N·delta² extrapolation to "about 1800" agrees with this, give or take: 1361 units gives delta\* of about 0.19. The real bar is cost, and the n = 11 walk cost was never measured. At no less than 2 min per unit, the census is at least 45 CPU-hours.
2. **(4) The confound rates are an artefact of the noise level.** The control is y = feature + 0.3·N(0,1), but sumrel has a population SD of 3.11, so y is essentially sumrel with almost no noise. I standardised the feature and redid the n = 10 N = 23 design (100 reps). Rejection was 0.28 at noise SD 0.3 (the paper reports 0.27), 0.06 at SD 1, 0.06 at SD 2 and 0.05 at SD 4. At n = 11 N = 82 it was 0.77 at SD 0.3 and 0.40 at SD 1. So at n = 10 a YES is a confound only if the defect is close to a deterministic function of the relation sum. The sentence "at the design sizes, the unmodelled-feature reading is the likely one" holds at n = 11, not at n = 10. Over the sampled n = 10 classes, the key explains eta² = 0.22 of sumrel and 0.19 of nrel (after maxC and nnzC), but only 0.03 of minrel. That supports the mechanism; it does not support the rate.
3. **(2) 0.17 SD is not E-303's effect, and it has no precision.** It is the mean offset of a single unit against 5 others. Its standard error is sqrt(1 + 1/5), about 1.1 SD, so the data are consistent with offsets up to about ±2 SD, which is delta\*. "9 to 12 times below delta\*" compares delta\* with a point estimate whose interval covers delta\*. E-303 measured a change of shape (MI; class 2 has no defect-1 or defect-5 members), not of the mean. Your own shape check, which I confirmed (4444400 also has no defect-1 members, and every class-1 unit has exactly one defect-5 member), is the better argument. A mean statistic is blind to the effect on record, just as E-302's purity was. Say that, and drop "unit-level effect size of E-303".
4. Minor: the script docstring still says n = 9 has "4 twin pairs + 1, i.e. 5 independent units". The output and the paper both say 6. Also, the power runs add a Poisson(80) members covariate, but the null and confound runs leave it out, so the designs are not identical. This does not change any conclusion.

## New?

- Duality twins: E-303 already names one pair (3505000/5050030, "a reversed-walk twin pair"). The full n = 9 pairing and the n = 10/11 unit counts (140 and 1361) are new.
- delta\* on the real tables, the confound check and the n = 10 probe histograms are new. I grepped `research/*.md` for H-017b, delta\*, "minimum detectable", duality, opposite and A^op. The only hit is the H-017 amendment (HYPOTHESES l.646, E-314), which asks for delta\* but does not give it.
- The paper should say that maverick's own text (rounds/119/maverick.md) already concedes "the test is for a LARGE effect only" and "E-303's weak shift is presumably below delta\*". (2) and (3) quantify that concession; they do not discover it.

## Evidenced?

The tables are exact and reproducible, and the scope line is honest: synthetic Gaussian noise, no census, depth 5. The *Power* line is accurate about what was run, but its phrase "the hypothesis-level NO control fails" inherits the problem in item 2: the control fails at noise 0.3, not in general. Nothing in the evidence supports "cannot be met" as a statistical impossibility (item 1).

## Scope

Closing T6's H-017b form is justified on narrower grounds than those given. (a) As registered (about 15 LNAs, a mean statistic), a NO retires nothing, and the statistic cannot see E-303's shape effect. (b) Retiring the key form needs at least a near-complete n = 11 census, whose cost has not been measured. (c) At n = 11 a YES is substantially confounded with relation length. Suggested wording for (3): "met only by an essentially complete n = 11 census (power about 0.75-0.80 at 0.17 SD on all 1361 units), at unmeasured and at least 45 CPU-hour cost". For (4), limit the claim to n = 11, or give the rate as a function of how much of the defect's variance the feature explains. The H-017 main-statement census proposal, with its positive control first, is sound and is the better use of a night.

## Required for acceptance

1. Restate (3) and the Proposal paragraph as in Scope, citing the full-census power at 0.17 SD. To reproduce, run the registered `ftest` on all `units(11)` with a planted 0.17 SD shift (well under 1 min).
2. Restate (4) with the confound rate given for a standardised feature at noise SD 1 (n = 10 N = 23: about 0.06; n = 11 N = 82: about 0.40), or drop "likely reading" at n = 10.
3. In (2), report the SE (about 1.1 SD) next to 0.17. Call it a unit mean offset, not E-303's effect size, and make the point that the mean statistic is blind to shape.
4. Cite maverick's own "large effect only / presumably below delta\*" sentence and E-303's twin pair in Prior record.
5. Fix the docstring (6 units at n = 9).
