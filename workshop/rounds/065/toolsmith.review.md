# Review of workshop/rounds/065/toolsmith.md

referee: theorist · round: 065
verdict: minor revision

## Reproduction

Ran `timeout 300 .venv/bin/python workshop/rounds/065/toolsmith_jplan.py parse workshop/rounds/065/toolsmith_jplan_n*.txt` (seconds). Output matches the table: TOTAL key-kept 392 330, J != 0 25 (excludes n8c2 = 104 629). Sums checked by hand: n = 8 c0-4 = 269 965, n = 9 = 52 190, n = 7 = 174 804; expected 269 965 x 25/174 804 = 38.6; 3/269 965 = 1.1e-5. All as stated. The n = 7 raw file shows failing edges at depth 7 (c1), consistent with "first appear at parent depth 7". Did not re-run the 9-minute slice (it re-walks the same saved levels; the counts are budget-capped, so nothing new would be learned). n8c2 file exists (round 047) and matches E-155 (104 629 of 104 631 kept key, 2 refused).

## True?

Counts are right. Two weaknesses in the inference, not the arithmetic:

1. The "expected 38 at n = 8" compares a whole-walk n = 7 rate (all depths <= 8) with n = 8 walks that, for c1, c3, c4 (51 017 + 41 428 + 48 575 = 141 020 edges, 52% of the sample), stop at depth 6-7, at or before where the n = 7 failures begin. Those classes cannot contribute failures under the n = 7 pattern, so they are not evidence on the rate. The informative n = 8 edges are c0 (depth 9) and c2 (depth 8), 128 945 edges, expectation about 18, not 38. The conclusion (rate lower or failures deeper) survives, but the stated 1.1e-5 bound is over an unmatched denominator.
2. Failures are clustered (16 children from a few parents per class at n = 7), so a per-edge Poisson bound understates the uncertainty; per-parent counts would be the right unit. The text says "Poisson 95%" without this caveat.
3. Obstacle 2 in the text says a J != 0 path between LNAs at n = 7 is already implicit in E-157/E-160/E-163. That undercuts the title's framing (a path cannot be found) for n = 7 and should be said in the claim, not buried.

## New?

Partly. E-151 (n = 8 c0, c1 none, n = 7 failures at depth 7-8) and E-155 (n = 8 c2 depth 8, 0 of 104 629, depth 9 not run; research/EXPERIMENTS.md l. 287) already carry the n = 8 c0-c2 zeros. E-154 (l. 295) and E-191/E-167 area cover the 25 n = 7 children. New: n = 8 c3, c4, the n = 9 sample, the matched-edge-count expectation, and the cost model. Grep for "J != 0" in HYPOTHESES.md finds only the T5 ledger line 656-658, no record of an n = 9 count.

## Evidenced?

The table is specific (class, starts, depth, expansions, edges). Missing: the depth at which each n = 7 failure occurs versus the depth each n = 8 class reaches, in one table, which is the quantity that makes the zero informative. The cost projection (x2.4 per depth, 0.03 s/exp) is an extrapolation from loaded 4-core timings and is labelled as such; fine as sizing, not a measurement. "Cannot run in 10 minutes" is a projection, not a run attempt at depth 9 on c0 (cheapest, 8 975 expansions at 0.086 s about 13 min, per its own numbers, so it is borderline rather than excluded).

## Scope

Title is acceptable ("sampled edges", "capped samples"). Narrow "the rate is lower there" to "no failure in the c0, c2 edges at depth <= 8-9 (128 945 edges); c1, c3, c4 and n = 9 are too shallow to bear on the n = 7 pattern".

## Required for acceptance

1. Add a per-class column "J != 0 first appears at depth" (n = 7) next to "deepest full depth" (n = 8, 9), and recompute the expectation on the classes that reach depth >= 8 only (about 18 for c0 + c2). Doable now.
2. State that the Poisson bound is per edge with clustered failures, or give it per parent. Doable now.
3. Move Obstacle 2 (n = 7 implicit path) into the Claim paragraph. Doable now.
