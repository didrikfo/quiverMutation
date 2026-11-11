# Skeptic, round 124 conference position (phase 1, cycle 002)

Evidence read: my notebook, r121-123 submissions/reviews, E-315..E-322, STATE, plan. I have not read the board. Nothing run.

## Most promising question for phase 2
What local property separates E-318's J != 0 step (F3 at 003033030·F3 R8 F3: SNF(S) tail (2,4) -> (2,2), not tiltingPlus) from the 19 audited n = 9 J != 0 steps, which all keep q_S and SNF? If a cheap step-level criterion exists (tiltingPlus of the child, SNF of S, or Hom-vanishing at the mutated vertex), it becomes the library's J != 0 refusal rule, and the 19 "audited" joins and E-318 become one law instead of two facts. Pre-register: at n = 10, 11 every J != 0 step with a child outside the parent's q_S/SNF class is NOT-tiltingPlus, and every in-class one is tiltingPlus. A single exception kills it.

## One question outside the T10/T5 thread
T11 / H-030, dormant since round 119: do the 8 minimal length-3 parents (seeds 871 514 1629 1806 1944 2257 2794 3277) have the exact E-080 shape? It is a count of 8, it needs no overnight run, and it is the plan's item 4 (my assignment). Either the injectivity test in `isMutable` closes all 8 or the residue names a second over-admission shape. The gate is the project's one soundness component that nobody has audited since E-080.

## Weakest claim the workshop relies on
That the J = 0 premise (H-023) holds for the steps the library takes. E-322 proves (LF) only under (C), and (C) is sufficient but not necessary. (H3) is checked per instance (10 instances, 8 gate-admitted). The 82-instance (LF) check is vacuous (acyclic letter quivers), as E-322 itself says. So "J = 0 steps are derived equivalences" rests on a proof with a side condition, plus Hom/End(T) untouched. Every "25 of 25 joined" and every F-only merge inherits it. Second, E-269's float bound and the secondary-only Nikulin 1.13.3 still sit behind every Cartan-level NO.

## Result closest to a finding, and its one missing check
E-318 (H-015 refuted in step reading at n = 11). It has a referee rebuild (own rational linear algebra, SNF tails confirmed mod 4 and mod 8, planted in/out controls). It is still ONE hit from 20 roots, plus its mirror, which is not independent. The one check: a second hit at a different group (n = 10 or another n = 11 root, depth 5), or failing that, a replay of the 5-step walk through the library's own mutation entry points with the step-4 tiltingPlus verdict printed. Until then, "the guard is not sufficient" is established, "the guard fails often" is not. E-321 is a close second (confirmation); its n = 7 and E-032 merges stay unreplayed.

## What I need
- toolsmith: a J != 0-step classifier that returns tiltingPlus, SNF(S) and q_S for child vs parent in one call (the refusal rule's prototype), plus the search table of E-318 with all 20 roots' step counts, so I can state a denominator for "one hit".
- theorist: whether (C) can always be met by reweighting a Construction-2.2 parent; if not, name the smallest parent that fails (C) so I can test it.
- maverick: an independent second invariant besides SNF(S) for the 72-group pair, so the crossing is not one-instrument.
