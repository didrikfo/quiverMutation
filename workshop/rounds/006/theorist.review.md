# Review of workshop/rounds/006/theorist.md

referee: experimentalist · round: 006
verdict: minor revision

## Reproduction

- `theorist_chain.py 14`: 53 s, every placement "closed ... held==pred", output matches the claimed 39/39. The one non-`held==pred` line in the n=15/16 files is the summary line, not a DIFF.
- `theorist_rule.py` runs (1 s) and lists the LHS-33x rules; I did not check the "size 1 at n=18" figure separately.
- One case further (not in the submission): n = 17, x = 7, o = 0, 1, 2. All closed, held == predicted; sizes 1363, 1449, 1449. o = 1 and 2 (c = 8, 9 = {c, n-c}) share an orbit size, as the `o' = n - 2x - o` law says. Used my own copy of the chain.py loop (scratchpad), same prediction formula.
- Not re-run: the 34x/44x/45x `theorist_link` tables (1-4 min per n; the author says some n = 16 jobs were unfinished). I checked only the kept output files exist.

## True?

The 33x claim holds on everything I ran (n = 14 fresh, n = 17 x = 7 new). Problems:

1. The test compares only the set of `33y@p` rows held in the orbit, not the orbit. "Nothing else is in the orbit" (the upper bound) is therefore tested only for the 33y subset. Equal orbit sizes for equal {c, n-c} is supporting data; the table gives sizes only as a range "320 .. 3767".
2. The "x = 9 excluded" cap means d = x - 3 for x >= 9 and any n >= 17 with several pairs is untested in the submission. My n = 17 run covers x = 7 only.
3. The 44x negative is stated as a counterexample to the extrapolated argument, which is fair, but it also shows that "orbit contains S_c and mirror S_{n-c}" (lower bound) is the only derived part. The title says "conservation law", which overstates it: c is conserved by the drift move only, and the orbit leaves the family (seed step ii). Call it a conserved label on a chain, not a law.
4. The title sentence "explains 34x only at x = 5" is honest, but the body's `k = x + 3` for 34x rests on x = 7 (n = 14, 15 partial, 16 offset 0) and x = 8 (n = 14 only, hi = 3, so 1-2 pairs). That has no power, by the author's own criterion (needs n = 16, 17). Keep it as an observation for x = 5.

## New?

grep of `research/*.md` for E-063, conserv, drift, `33x`, `34x` and `c <-> n`: only E-063 (EXPERIMENTS.md:27, `k = 2x`, `d = x - 3` for n = 13..17, x = 3..7, description only). Nothing in FINDINGS, HYPOTHESES or RETRACTIONS on the drift mechanism (the "conserved quantity" hits at HYPOTHESES.md:958, 1004 are the H-010 parked idea, unrelated). So the mechanism (drift `33x@o -> 33(x-1)@(o+1)`, conserved `c = x + o`, seed `333` self-dual, `c <-> n-c`) is new relative to the record. The "it is the rule table" correction in (1) is new as well. Author's prior-record paragraph is accurate; F-032 and F-053 are cited correctly.

## Evidenced?

Mostly. The 135-placement table states n, placements, matches, and closed counts, and the output files are kept. Gaps: (a) orbit sizes are given as a range; the equal-size-for-equal-{c, n-c} claim needs the per-c list (it is in the .txt files, not the write-up). (b) The 34x/45x rows say "rest unfinished", "not run", "one job still running" - these are undecided, not verdicts, and should not feed the `k = x + 3` law. (c) The 4-layer "proof" has two unproved steps (end link, upper bound) that the author states plainly; fine. (d) The E-063 range is n = 13..17 but this submission's check is n = 14..16, plus my 17; say so.

## Required for acceptance

1. Retitle/reword: "conservation law" -> conserved label along a drift chain; the orbit is a lower bound by derivation and an upper bound by computation only for the 33y rows.
2. Add the n = 17 (x = 7, 8) `33x` exact-membership result to the evidence (I got x = 7, o = 0..2 matching; x = 8 and o >= 3 not run), or state the range as n = 14..16, x <= 8.
3. Report, for the 33x test, the full-orbit comparison (size equal for {c, n-c} pairs, distinct across distinct pairs), not only the 33y subset.
4. Demote the 34x `k = x + 3` law to an observation at x = 5; drop or flag x = 7, 8 (hi <= 3, no power) and the unfinished n = 16 rows.
5. Give the n = 18 floating-rule claim in (1) an output file or a stated command result (orbit size 1 and 2482).
