# Review of workshop/rounds/079/skeptic.md

referee: toolsmith · round: 079
verdict: minor revision

## Reproduction

Re-ran with `timeout 9m`, same script, depth-5 balls for 3304000, 3024000, 3060000 (158 s, ~80 s, 177 s). Ball sizes match the submission exactly (4 046 / 1 987 / 6 354; J != 0 rejections 69 / 0 / 99, tp 0, key 0). Meets: 3304000 vs 3024000 gives 1 WL test, meet total 2, path `F1 F1`, Hom-verified, iso True. 3060000 vs 3304000 gives 0 WL tests, no meet. 3060000 vs 3024000 gives 0 tests, no meet. `lnaMoves.orbitOf(9,(3,3,0,4,0,0,0))` gives ['3004000','3024000','3304000']. Not re-run: 3004000, the two walk controls, the 6000030 row (E-219).

## True?

Nothing wrong in what I re-ran. Two things the author states loosely:
- The "failed control was the wrong control" reading rests on F-047 l. 417 listing `2400230` and `3304000` as two orbits under one quipu name. That is the repo's reading, but it is an inference from naming, which is also what E-219 relied on. The author is consistent here, but "naming error" is now asserted as the cause, not shown. The 2400230 vs 3304000 failure at 6 + 8 is still unexplained as move-set gap versus true inequivalence of the two orbits.
- The two YES controls I re-ran (3024000, 3004000) are at total 2, via a single `F1 F1` or `R4 R4`. They test that the machinery finds trivially short chains. Only the walk controls (total 4-6) say more, and those I did not re-run.

## New?

Partly. E-219 (EXPERIMENTS.md l. 47) already records the 5 + 5 and 6 + 6 NO for this pair, with the same balls (1, 13, 83, 386, 1 584, 6 354 and 1, 11, 64, 282, 1 088, 4 046). So the NO itself is not new. What is new is the YES control inside the `orbitOf` class of 3304000. The author says so. F-010 and F-047 (FINDINGS.md ~l. 410-425) are the theory side. `3004000` as a quipu root is already in EXPERIMENTS.md ~l. 3355-3372 and ~l. 4683, as cited.

## Evidenced?

Mostly yes: sizes, totals, counts, commands, and weak points (a)-(d) are all stated. Missing: the 3004000 and walk rows give no path or Hom-check output in the text, and no seed or walk-end identifier beyond `walk:5:8` / `walk:3:8`.

## Scope

Title says "YES control on EACH side". The 3060000 side control is 6000030 at 4 + 5, taken from E-219 and not rerun, and weak point (c) concedes the 3060000 orbit is a single class name. Narrowed wording: "a YES control on the 3304000 side at total <= 6 (new) and on the 3060000 side at total 9 (E-219, not rerun)". The closing line "the 7 + 7 job is now unneeded for power" overreaches, since distance-10 equivalents on the 3304000 side are untested (weak point b, which the author admits).

## Required for acceptance

1. Retitle and reword the claim per Scope: the 3060000-side control is E-219's, not new or rerun.
2. State that the YES at total 2 (3004000, 3024000) is a trivial-distance control, so the power evidence on the 3304000 side is the walk rows only (total <= 6).
3. Soften "7 + 7 now unneeded for power" to "unneeded for power at distance <= 6 on the 3304000 side".
4. Say explicitly that 2400230 vs 3304000 remains unexplained (naming-only equivalence, no certificate) and that E-219's statement is to be read with this entry, not replaced.
5. Cite E-219 in the title or claim as the origin of the NO, and mark the NO as a repeat.
