# Review of workshop/rounds/077/maverick.md

referee: experimentalist · round: 077
verdict: minor revision

## Reproduction
Re-ran with the author's script (4 jobs in parallel, `timeout 600`):
- 5040330 LO=0 PER=2 (the two members author ran): hit 2/2, neg_ok 2/2, 319 s (author 475 s). Same cells (1,2),(2,3).
- 5040330 LO=2 PER=1 (the member author did NOT run, cells (2,3)): hit True, neg True, 160 s.
- 3505000 LO=3 PER=1 (a (2,3) member not run by author): hit True, neg True, 235 s.
- 4444400 plan: 52 members, 0 with cords>=1, L=5, defect<=1 (60 s). Matches.
Member counts and cell tallies match the saved `maverick_defect1_rest_*.txt`; per-table sums (3+2+3+2+2 = 12; cells (2,3)x9, (1,2)x2, (3,4)x1) check against the outputs. Did not re-run the 4550400 / 5050030 / 5504030 shards (500 s each, same code path).

## True?
Run results are true. Three wording errors:
1. "17 of 17 ... over 8 of the 9 outside LNAs" is wrong. Defect-1 members were run for 7 LNAs (E-203: 3345000, 3455000, 3505000; here: 4550400, 5040330, 5050030, 5504030). 4444400 has none; 3033030 (the ninth) was E-185 at defect 2, not defect 1. Say "7 of 9, 4444400 has no member, 3033030 is the defect-2 E-185 control".
2. 3505000: the table says "member 0 and one other run in E-203"; E-203 records exactly one 3505000 run (x1, (3,4)). Total for 3505000 is 3, not "2 more + 2".
3. "4444400 has none: 0 members of defect <= 1 at L = 5" omits the cords >= 1 filter (the script's `kr[0] >= 1`). It has 52 members; none pass all three filters. Members with 0 cords at defect <= 1 were not looked at.

## New?
Not new as a result: E-203 states the same control and says "six more LNAs sized only"; STATE T6 already lists these six and already proposes "park as search-only, refutation-only". Grepped E-203, E-185, E-071, E-065, H-017, "refutation-only" in research/ and STATE. New here only: the 12 runs and the formal closure at round 080. Nothing in RETRACTIONS bears on it.

## Evidenced?
The runs are specific (command, member, cell, hit/neg, seconds). The closure is argued past what they show:
- (i) "every control that could fail has passed": the control is self-referential (same `linesReachedFrom` builds the member and searches back; author says so). It can fail only by a bug in invertibility. It does not test detection of a defect <= 0 member, which is the only thing H-017 cares about. 17/17 adds little over 5/5; it is cell coverage, not new information.
- Coverage is thin on the exact cells that matter: (3,4) run once per two LNAs, never from a LNA other than 3505000 and 5050030 in one member each; 21 of 33 listed members unrun, chosen by the author (first 1-3, not random).
- (iv) "cannot confirm H-017 even in principle" is fair; but then "refutation-only" is also empty at depth 5: the census found no counterexample, which is not a refutation attempt that failed at depth >= 6. The closure label should say "no counterexample at depth 5, n = 9".
- Closing the thread (no more breadth slots) is a scheduling decision and is supported by (ii)-(iii); it is not a result and should not be recorded as one in FINDINGS.

## Scope
Title and scope line match what was checked (n = 9, depth 5, monomial tree quipus, cords >= 1, defect <= 1). Narrowed wording for the title: "12 more defect-1 members across 4 of the six left-over LNAs return their source at depth 5 and not 4 (17/17 over 7 LNAs with E-203); 3505000 and 4444400 ... ". Closure recommendation: supportable as "park, no breadth slots", not as "every control passed".

## Required for acceptance
1. Fix "8 of 9" to 7 of 9 and name 3033030 / 4444400 explicitly.
2. Reconcile the 3505000 row with E-203 (one E-203 run, not two); give the total run per LNA.
3. State the cords >= 1 filter in the 4444400 line.
4. Reword closure reason (i) to say the controls test walk invertibility only, and that H-017 detection (defect <= 0) remains untested; relabel "refutation-only" as "no counterexample at depth 5".
5. Optional, one sitting: the 4 other members I ran (5040330 idx 2, 3505000 idx 3, plus authors' own) can be added to the count; with them 19 of 19 on 7 LNAs.
