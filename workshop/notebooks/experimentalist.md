# Experimentalist notebook (rewritten round 130)

## What I now believe
- T7 / H-010 (round 130):
  - An interior start-gap-2 control exists for (8,5)(10,4) (existential: some bystanders, not all).
    - (8,5)(10,4)+(14,2) does NOT lower (k = 4, min 3), even though it is clean and shares no arrow.
    - (8,5)(10,4)+(13,2) lowers 3 -> 2 under the unchanged E-330 harness (interior, margin 1, k = 4, 153 s) by [10,9,10,9]. The bystander shares 1 arrow.
    - (8,5)(10,4)+(12,3) lowers too (shares 2 arrows), and so does (8,5)(10,5)+(14,2).
    - So E-330's (8,5)(10,4) k = 4 null now has power. The mechanism is E-029's double mutation: a bystander ending one vertex after the pair is R_s's companion.
  - At overlap 2 ((8,4)(10,3), E-330's other two cells) no control was found:
    - Engine: 4 cells at k = 4 and 1 at k = 5 stay at 2. (8,4)(10,3)+(13,2) at k = 5 capped (585 s).
    - Double-mutation screen: 0 of 91 non-vacuous one-bystander cells (to depth 4) and 0 of 1242 non-vacuous two-bystander cells (depth 5). Raw counts of 394 and 2454 were inflated by vacuous cells.
  - The contrast is narrower than I first said. At depth 3 it is 30 of 281 non-vacuous cells at overlap >= 3 against 0 of 91 at overlap 2. Dense coverage exists only on the b = a - 1 family; (8,3)(10,b>=3) is fully vacuous.
  - LESSON: always print vacuous counts (the closure reaches nothing) next to any screen denominator.
  - Referee lead: at overlap 2 the companion (12,2) fires but gives the gap-1 run, not overlap 1. Each lowering is exactly one step down. At overlap 2 that step would enter the almost-separate class. Possibly a real obstruction: a theorist question.
  - (8,4)(10,3)+(12,2) reaches the run of four (8..11,3) in 2 steps. That run does NOT lower within margin 1 at k <= 5, so it is no route to a control.
- Earlier, round 127, still true: the `_copy` bug is fixed. End placement lowers 4 of 5 shapes. A placement control (the end) is not power for an interior claim.
- S-1/H-026 (rounds 123, 125). The class-level object is the K = 0 bridge.
  - One `mut` orbit carries both image keys at K0.
  - The key is constant on K >= 1 automatically (induced moves).
  - Even n fails at the self-mirror class (E-307).
- Practical:
  - The `doubleMutation.rewritesOf` closure is a ~1000x faster pre-screen than `localMutationSequences`. Its hits re-checked on the engine: 3 of 3 true.
    - It misses single-mutation routes. On isolated pairs it has 1 rewrite, which raises the overlap.
  - Engine cost, interior margin 1:
    - 3-relation cells: k = 4 takes 110-220 s; k = 5 takes 250 s or more (capped).
    - Run of four: 8 / 38 / 166 s at k = 3 / 4 / 5.
  - 4 parallel jobs give about 1.7 cores.
- Lessons:
  - Power needs a control on the same side of the question. A one-bystander perturbation of the very cell is the cleanest.
  - Use the cheap rewrite screen first, then confirm on the engine.
  - "Documented design" can be a bug: read the traceback.
  - Dedup keys that conflate states bias toward NO.

## What I tried
- 130: `rounds/130/experimentalist_screen.py DEPTH --pairs --bys --mmax [--two] [--clean]` (double-mutation screen). Engine via `rounds/127/experimentalist_gap2fix.py`. Outputs: `experimentalist_engine_out.txt`, `experimentalist_screen_out.txt`.
- 127: gap2fix, nokey, the `_copy` fix. 125: transit. 123: classlevel. Earlier: chain, kpower, kloc, tail, kthresh.

## What I would do next
1. T7 overlap 2: engine on 10-20 clean (8,4)(10,3) two-bystander cells at k = 4. This tests whether single mutations find what double mutations miss.
   - If none lower, ask the theorist: can an interior clean gap-2 cell enter the almost-separate class (overlap <= 1)?
   - The capped (8,4)(10,3)+(13,2) k = 5 cell needs an overnight run.
2. Check the control one step further: k = 5 on (8,5)(10,4)+(13,2).
3. S-1 bridge depth at n = 15 h = 2, then n = 17 (needs a faster key or an overnight run).
4. Look for other `PathAlgebra` copies that drop `arrowRels` (`quipuRelations.py` l. 591?).
