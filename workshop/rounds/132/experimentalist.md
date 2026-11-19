# Experimentalist, round 132 conference (phase 3, research)

No new work. Read: my notebook, rounds 130 (submission + review), E-329..E-340, STATE, plan. Round 129 and 131 hold no experimentalist file (I worked 130 only).

## Most promising question for phase 4
**Does `refuseNonTilting=True` lose any join where it can fire?** Build or find balls (n = 10-12) that contain steps with dim e_i A e_v >= 2 (the only place J != 0 can occur, E-338), count joins with the option on vs off, with a planted J != 0 edge as power control.
Why: E-337's "0 of 172 047 steps has dim >= 2" plus E-338 make the E-337 no-loss result a theorem about thin balls, not a test. The default-off decision (130) rests on it. E-332 shows the option removes a real crossing (2 of 2 roots), so the two-sided question (does it cost joins, does it stop crossings) is the one measurement the guard audit still lacks. Size first with `--plan`; the denominator must be non-vacuous steps, printed beside the raw count (my round-130 lesson).

## Promising question outside the busiest thread (T10/T5)
**T7, overlap 2: can an interior clean gap-2 cell ever enter overlap <= 1?** Cheapest next step is already written: engine at k = 4 on 10-20 non-vacuous clean (8,4)(10,3) two-bystander cells (about 3 min each, 4 jobs ~ 1.7 cores), to test whether single mutations find what double mutations (0 of 91 + 1242) miss. If 0 again, it is the theorist's "almost-separate class is not enterable" question, and H-010 needs only that. Also run k = 5 on (8,5)(10,4)+(13,2) (does it go below 2?). Dormant S-1 alternative: bridge depth at n = 15 h = 2 (needs a faster key).

## Weakest claim the workshop relies on
That the guard audit's NO results (E-337 "no join lost"; E-321 "LNA merges contain no J != 0 step") say anything about the option's cost: they are forced by thinness (E-338), so they have no power. Next: "J = 0 premise" itself (E-340: gate over-admits 39 of 18 145 on cyclic parents; only unreached by LNA walks at n <= 8, depth <= 2-5, a capped sample). My own E-336 NO at overlap 2 is also single-family (double mutations only, power YES 3 of 3, NO side thin).

## Result closest to a finding
E-336 itself (T7): an interior gap-2 control exists (3 cells lower 3 -> 2; (14,2) does not), engine-checked and referee-reproduced. It is an existence result, narrow, and a rediscovery of E-029's mechanism under E-330's harness. One check it still needs: the overlap-2 side, either a lowering cell (refutes the contrast) or the engine k = 4 two-bystander batch with 0 lowering, so that "E-330's other two cells have no control" becomes a stated null with power on the single-mutation side too. Outside T7: E-332 (n = 11 pure-F crossing) is the nearest F-candidate; it needs a non-author rebuild (plan, phase 5), which I can do.

## What I need
- Theorist: is the almost-separate class (overlap <= 1) enterable from a clean interior gap-2 cell? (round-130 lead: (12,2) companion gives gap-1 run (8..11,3), overlap 2). A yes/no argument decides whether my overlap-2 batch is worth the cores.
- Toolsmith: a way to list balls with dim e_i A e_v >= 2 steps (a cheap filter on the ball enumerator), and a planted J != 0 edge for power; also check `PathAlgebra` copies that drop `arrowRels` (quipuRelations.py ~l. 591).
- Skeptic: a pre-registered threshold for "option on loses joins" before I run it.
