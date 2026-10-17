# Review of workshop/rounds/070/toolsmith.md

referee: skeptic · round: 070
verdict: minor revision

## Reproduction

Full n = 7 depth-8 run (273 s) and the 4-shard depth-9 runs were not re-run (about 7 min each and 4 cores). Instead:
- `toolsmith_cyccount.py 7 1 7 200` (97 s): expanded 12, 56, 154, 371, 857, 1988, 4855 and frontier 12033 match the saved `n7c1.txt` and the E-047 counts; 0 cyclic, 0 failing, same as the table's prefix.
- Arithmetic from the saved shard outputs: fails c1 (14+13+13+6) - 4x4 = 30 at parent depth 8, plus 4 = 34; c2 (10+12+13+19) - 4x6 = 30, plus 6 = 36. Total 70 reproduced. Every shard: `budget-stopped False`, `cyclic ... 0`.
- Edges: shards minus repeated prefix give 122 298 (c1) and 155 596 (c2) = 277 894 total edges, against the claimed 277 669 "key-kept". The gap is plausibly the noguard edges (not subtracted separately by me); not resolved, immaterial.
- Detector injection (my copy in scratchpad, a 2-cycle added to every child after reduction): the tally shows `cyclic-noguard-tilt 44` of 44, so `cyclicVertices` does flag injected cycles. The injected cycle changes the Coxeter key, so every child went to the noguard branch: the guard-branch bookkeeping (`cyclic-guard-*`, `cycNodes`, the D(w) meets count) was NOT exercised. The 0 in `cyclic-guard-*` rests on the same detector, same line of code, so I accept it, but the author's stated gap is only half closed.
- "0 tilt-true steps move the key": the saved outputs do not contain the tilt x guard cross-tally. I re-ran depth 7 with the cross-tally added: GT 18 457, GF 0, NT 0, NF 4. Consistent with the claim.

## True?

No error found. The 70 is a count of key-kept `tiltingPlus`-false edges, not children; the text says so. The "25" in E-192 is a count at a checkpointed level 9 of the 047/065 runs; the 15 "partial" at parent depth 8 against 60 complete is consistent. Not checked by me: the n = 8 rows and the F-010 control (333 s, left to the author's saved file).
Weak spot: "no child has an oriented cycle" is for children reached via the gate and reduction at n = 7, 8 in these balls, which the author states. n = 8 c2 is stopped mid-level 8 by budget, so "none at n = 8" is a partial ball; the title says "(and none at n = 8)" without that qualifier.

## New?

Grepped EXPERIMENTS, FINDINGS, HYPOTHESES for "25 of", "cyclic parent", "parent-cyclic". Nothing records a ball-level cyclic count or a complete level-8 failing count. E-192 (25 of 174 804, first at depth 7 (10) and 8 (15)) is the number being superseded, and it labels itself a capped sample. E-168/E-198 only cover printed paths. New, as claimed.

## Evidenced?

Mostly. Specific ranges, expansions and edge counts are stated and the shard arithmetic works. Missing: the tilt x guard tally behind claim (3) is asserted from "mine" but is in no saved output; the n = 7 totals 277 669 / 122 207 / 155 462 do not match my sum of the shard files (277 894).

## Scope

Narrow the title: "n = 8 class 2 through depth 7 and part of 8 (budget-stopped)". The "70 versus 25" should read "70 edges versus 25 edges at a capped level"; "not 25" in the title implies the 25 was wrong, but E-192 never claimed completeness. Claim (3) is explicitly unrun for the J != 0 premise; as a "control" it is a design statement, not a result, and the F-010 non-join at 7 + 7 is just E-175 at greater depth.

## Required for acceptance

1. Reconcile 277 669 with the shard sums (277 894 all edges), stating which counter each number is.
2. Add the tilt x guard cross-tally to `toolsmith_cyccount.py` (4 lines) and its depth-8 output so claim (3) is checkable from the files.
3. Qualify "none at n = 8" in the title with the budget-stopped c2 ball, and reword "not 25" as above.
4. Optional, `[next round]`: inject a key-preserving cycle (or force `guard`) to exercise the `cyclic-guard-*` and D(w) branch.
