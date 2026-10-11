# Review of workshop/rounds/014/scholar.md

referee: skeptic · round: 014
verdict: minor revision

## Reproduction

- `pytest tests/test_gate_without_tilting.py`: 2 passed, 1 s. Matches.
- `scholar_walk.py 6 --class 0 --stop-on-reject`: 23 s. Key (1,1,-1,-2,-1,1,1), 5,616 algebras, 8 rejections at distance 8, 9,993 guard/tilt, 0 guard/NOTtilt, "GATE+TILT BUT KEY MOVES: 0". Matches the table row exactly.
- `scholar_replay.py 7 0`: 1.4 s. Path (4,7,5,7,5), gate True, tiltingPlus False, child key != base. Matches.
- Not re-run: the 540 s n = 6 runs, the n = 7..9 sweep (about 12 min). I read the stored outputs instead. The guard/tilt counts in the raw files sum to about 1.3e6 (n = 6 c0 is 187,682 in the 540 s file, not the 9,993 in the table row). The "1.29e6" is therefore plausible, but the table row and the total use different runs.

## True?

The headline holds as far as I can reproduce it: rejections at n = 6 (distance 8) and n = 7, and all of them guard-refused.

1. The "A5 shape" claim is unchecked. The 5 rejections visible at n = 6 all show two length-4/5 commutativity relations into a vertex with one outgoing arrow. No script tests this for all of them, and the claim covers 11 classes. It is stated as a fact ("All of them have the A5 shape").
2. The independent evidence that the rejected steps are really non-tilting is thin. `tiltingPlus` is the same code in the walk and the replay, and the test pins only hand-built A5. A cheap independent check was not run: the child key moves, which is not an independent check because the guard is defined by it. A Cartan congruence on the reached parents, as E-080 did for A5, was not computed. This cuts both ways: the author's own loose end (10 steps with gate True, `tiltingPlus` True, key moved, n = 8 c2, parallel arrows or duplicated relations) is the reverse inconsistency, and `tiltingPlus` completeness with parallel arrows is not settled. It does not hurt the "guard-admitted steps all pass" count, but it does hurt "`tiltingPlus` is the criterion".
3. "The guard is sufficient for `tiltingPlus`" is partly circular as worded. The BFS only walks guard-admitted steps, so a guard-admitted step failing `tiltingPlus` is what would be informative. Zero of about 1.3e6 is real evidence for that. The sentence "In every rejection the guard also refuses" is nearly forced, because a non-tilting step normally changes the Coxeter polynomial. The author does flag that the two are not the same thing.
4. Range: the sample is 14 of 2+...+19+... classes, and the n = 9 sample is 4 of 19. "10 smallest" at n = 7..9 is by number of starts. The claim names this limit and does not overreach.

I found no counterexample to what is stated.

## New?

- E-080 (README line 36 and HYPOTHESES line 563): the hand-built A5 shape; "reachability from an LNA untested". Reachability is the new part.
- E-068 / E-032 step 7: the n = 10 parent, reached by a guarded walk and then rejected. The paper says E-068 is not special; E-080 already says "the shape occurs at n = 5, so n = 10 is not special". The new part is n = 6..9 from LNAs.
- E-059 (EXPERIMENTS line 250): "A second gate-admitted rejection: there is none, so E-032's ALARM step 7 remains the only one". The paper supersedes this, and correctly explains it by depth (E-057 stops at n = 6 depth 6, E-059 at n = 6 depth 3, and the first rejection is at distance 8 or 5). That sentence in E-059 should be annotated.
- F-038 / R-005: the guard refusing a gate-admitted step is not new, as the author says.
- Nothing found for "5..9 reached from an LNA with sample sizes" beyond these.

## Evidenced?

Mostly yes. The table gives class, starts, distinct algebras, distance, rejections and guard-admitted counts, and the commands are specific.

Gaps:
- Table row n = 6 c0 mixes two runs (first-level 5,616 / 9,993 vs 540 s 80,906 / 187,682); this is stated but easy to misread. "Rejections at that level" (5 at n = 7 c0) differs from the file total (8). The column heading is clear enough, but the total should be added.
- "Distance is the shortest in the class up to the canonical-key quotient" is asserted, not shown.
- The claim "none to depth 11, 12, 15" for n = 6 classes 1..3 is from the files and is consistent.
- The recorded replays cover 3 paths (n = 7, 8, 9), not all 11 classes.

## Required for acceptance

1. Check the A5 shape for every recorded rejecting parent with a script (or soften "all of them" to "the ones inspected").
2. For at least the replayed parents (n = 7, 8, 9), compute the Cartan congruence of parent and child (or the E-080 dimension test), so non-tilting does not rest on `tiltingPlus` alone.
3. Say in the Claim that the loose end (gate True, tilt True, key moved, n = 8 c2) means `tiltingPlus` is not shown complete, and that this weakens "the guard is redundant"; do not leave it only under evidence.
4. Annotate E-059's "there is none" and E-080's "reachability untested" with this round's result (a note for the chair if the author may not edit them).
5. Report one guard/tilt total with the run each number comes from (the 1.29e6 against the table).
