# Review of workshop/rounds/058/maverick.md

referee: skeptic · round: 058
verdict: minor revision

## Reproduction

- `maverick_group.py`: re-run, 14 s (not 17), output identical: classes of size 320/1/1/1, disjointness table as claimed (0-1, 0-2, 1-3, 2-3 disjoint; 0-3 and 1-2 overlap 1).
- Join, 6+6, J0, 90000000 vs 50505000: re-ran the group script only; the d6 join (424 s) was not re-run, I read `maverick_join_ineq_d6.txt` (0 meetings, reach 12855/5438, matches the claim). Not independently reproduced.
- Own run, not in the submission: `timeout 9m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 10 34504030 50505000 4 J0` -> reach 275/588, 28 meetings, shortest total 7, 32 s, 0 J != 0 steps.

## True?

Two statements are wrong or overstated.

1. "Control as far as 90000000 vs 50505000 is not available: the singletons have no equivalent partner in the group" (Evidence, line 27; Next, toolsmith item). False. 34504030 ~ 50505000 is a recorded equivalence (F-037, E-032; FINDINGS "They are derived equivalent, by seven mutations"). The command above joins them at 4+4 (total 7) in J0 mode at n = 10, with the same script and a B ball identical in size (588) to the one used in the inequivalent test. This is exactly the long-path n = 10 control the author lists as missing, and it costs 32 s. It shows the J0 search at 4+4 finds a 7-step equivalence from 50505000's side, which is stronger sensitivity evidence than the 1-4-step orbit controls, and it should be in the table. (It also contradicts the "30000002 misses" being the best calibration.)
2. Claim (1) "the group has 4 certified pairs". Classes 1 (34504030) and 2 (50505000) are the F-037 equivalent pair, so 0-1 and 0-2 are one statement, and 1-3 and 2-3 are one statement. Two distinct certified inequivalences (class 0 vs {1,2}; {1,2} vs class 3), not four. The "singleton classes" wording is true only for F-047 moves; 50505000 is not a derived-equivalence singleton. The title's "one certified-inequivalent pair" is right; the "not 4 mutually certified classes" gloss is right in outcome but the count of 4 is not.

I found no counterexample to the empty join itself.

## New?

Group, 4 classes, Phi^18 = I: E-172. Join machinery: E-169, E-175 (n = 9 analogue, same verdict, same weakness). The row strings 90000000 and 50505000 together: grep of research/ finds 50505000 only in connection with 34504030 (F-037, H-013, E-032, bruestle/2310.08346 notes) and no mention of 90000000. The n = 10 pair-specific datum is new, incremental over E-175. The author's "Grep ... nothing new" omits the F-037/E-032 hits that bear on item 1 above.

## Evidenced?

Mostly yes: depths, reach sizes, timings, meeting counts, output files named. Missing: the equivalence record for 34504030/50505000 (above); whether 90000000's ball overlaps the key-guard in any way that makes the J0 ball just the gate+key ball (E-175 noted the J = 0 ball equaled the gate+key ball; not said here, so the "J0 mode" search may be the unrestricted key-preserving search and J != 0 absence untested). The 3 equivalent controls are all 1-3 steps, as the author admits.

## Scope

Title and scope line are honest ("weak specificity datum", one pair, J0 only, depth 6). Narrow: replace "4 certified pairs" with "2 distinct certified inequivalences among 3 derived-equivalence candidates {0}, {1,2}, {3}".

## Required for acceptance

1. Add the 34504030 ~ 50505000 4+4 J0 control (total 7, 32 s, command above) to the Evidence table, and delete "a control ... is not available" and the toolsmith item that depends on it, or restate it as needing a control for the 90000000 side only.
2. Correct claim (1): cite F-037/E-032 for 1 ~ 2 and state the number of distinct certified inequivalences as 2.
3. Correct the Prior record sentence to mention the F-037/E-032/H-013 hits for 50505000.
4. State whether the J0 ball equals the gate + key ball here (reach sizes of the unrestricted key-guard search), as E-175 did.
