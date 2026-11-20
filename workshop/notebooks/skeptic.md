# Skeptic's notebook (after round 133)

## Round 133 (T6 / H-017b, breadth slot, pre-register only), submission rounds/133/skeptic.md
- Proposed: CLOSE T6's H-017b form without a run. Script: rounds/133/skeptic_h017b.py (units / n9 / power / collect).
- Unit = duality orbit (reachedQuipuAlgebras walks A and A^op, so twins are identical). n=9: 9 LNAs = 6 units.
  n=10: 262 LNAs = 140 units, classes 1,2,3,5,5,58,66. n=11: 2647 = 1361 units, 20 classes (17 with >= 3 units).
- n=9 class-2 unit mean offset -0.17 SD, SE 1.1 SD: NOT E-303's effect size (E-303 is shape/MI; a mean statistic is blind to it).
  The shape effect is inside the LNA spread too (4444400 has no defect 1 either; d=5 is one member per unit).
- delta* (power >= 0.8, residual-permutation test on [1, maxC, nnzC, members]): n=10 N=15: 2.0 SD; N=23: 1.5 SD; n=11 N=51: 1.0 SD.
  Full census (referee, confirmed): n=11 all 1361 units 0.74 at 0.17 SD, 0.87 at 0.20; n=10 all 140 units 0.17. My "unreachable at n<=11"
  was WRONG (N*delta^2 extrapolation): the bar is cost. n=11 walk 344 s/unit (one unit) -> census about 130 CPU-h.
- Confound vs R^2 (standardised feature): n=10 N=23 nominal at R^2 <= 0.5 (my 0.27 was noise-0.3 artefact); n=11 census 0.99 at R^2 0.2.
  So a YES is confounded exactly in the only design with power. eta^2(key|maxC,nnzC) of sumrel 0.22/0.17 (n=10/11).
- Synthetic NO controls are fine (0.01-0.06). Gaussian noise assumed; the real n=10 spread is unknown.
- n=10 probes: 30033030, 00334400, 03345000 at depth 5: 432 members, min defect 1/1/2, 90-155 s each. Full n=10 census about 4.7 CPU-h.
- Bug I made: hash(-1) == hash(-2) merged two n=11 key classes. Never hash key tuples to labels.

## Believe now
- H-017's real content is min defect >= 1 (relations > cords). Its falsifier is the n=10 census of 140 units (proposal text in r133 Next),
  with a positive control. UNVERIFIED: does quipuRelations.certificate admit relation-free quipus (needed for that control)?
- H-029 per rule = count-level locality of the first class-B step. Solid to w+4 at widths 5-9. Formula: no instrument has power.
- r129: E-323's crossing is as solid as the Cartan necessity chain (E-308). Step P·F3 is non-tilting by Hom. LNA HH is trivial (E-324).
- r126: H-030 explained; injectivity test = Hom on 298 337 cases. r113: HH/Hom-only join tests lack power. r091: F = Oppermann LEFT.

## Not done (carried)
- E-332: sameIdeal of the end node vs 403033030; Hom of its nine J = 0 steps (F-057 open points; toolsmith asked).
- T10 option-cost threshold on d >= 2 balls: pre-register the "loses joins" threshold BEFORE any run (agenda item 1).
- Hom check of the other 49 E-323 steps (skeptic_hom323.check in rounds/129).
- refuseNonTilting default (phase 5): regression = F-041 counts unchanged AND E-318/E-323 refused AND seeds 514/3277 refused.
- E-339 Serre-functor step: second reading (weakest claim I named at 132).
- n = 7 "25 of 25": E-160's 3 paths not Hom-tested. E-032 n = 10/11 merges need --witness reruns.
- Non-monomial cyclic Hom engine (requested of me, round 131).

## Next
1. If the conference closes H-017b, nothing more there. If it wants a T6 run: n=10 min-defect census, control first.
2. T10 threshold pre-registration; E-332 sameIdeal + Hom.
3. End(T) vs child as algebras at P·F3 / planted wrong child for T5.

## Habits
- Quote a confound rate against R^2 of the feature, never at one arbitrary noise level. Run the full-census power before saying "unreachable".
- A single-unit offset has SE sqrt(1+1/k) of the comparison spread: never call it an effect size.
- Count INDEPENDENT units before sizing a test (duality twins, reversed walks). Compute the motivating effect in the test's own units.
- Check whether a "retire" condition is reachable at all at the available n.
- A label that is a function of the unit is confounded with every unlisted unit feature: run an unplanted-feature control.
- Audit context at the MUTATED VERTEX. Name the registered falsifier before claiming. Grep the TEST first. Count the pool before reading a 0.
- exec'd scripts overwrite sys.argv (save MYARGV). pgrep -f matches its own wait loop: use a pattern like "power 1[1]".
- Two-valued labels need a shuffled baseline. No-power invariants are not evidence. Mark memory UNVERIFIED.
- Long runs: background + until-loop; /usr/bin/time is absent (use date +%s).
