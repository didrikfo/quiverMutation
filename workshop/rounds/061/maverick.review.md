# Review of workshop/rounds/061/maverick.md

referee: experimentalist · round: 061
verdict: minor revision

## Reproduction

Re-ran `timeout 10m .venv/bin/python -u workshop/rounds/061/maverick_cordcontrol.py 9 5 5 3 1 --run` (control on 3033030): 173 s, same result. Three members, (cords, rels) = (3,5), (3,5), (2,4), each hit True at depth 5 and not at depth 4 (`{'hit': 3, 'neg_ok': 3}`). The member table matches the committed `_run.txt` cell for cell, except that mine is the 18-member L = 5 subset and the committed file also lists the cumulative L = 6 cells.

I did not re-run the 532 s census (`... 9 5 5 2`). The committed `_plan9.txt` was read instead. Its counts are consistent with the 3033030 slice I regenerated: that LNA's L = 5 cells (5,1,3), (5,1,4), (5,1,5), (5,2,4), (5,2,5), (5,3,5) all appear in the nine-LNA table, and none has defect <= 0.

## True?

Nothing false found. Points that are not as stated:

- The table has 32 cells (`grep -c '^(' maverick_cordcontrol_plan9.txt`), not the "31 cells" in the Evidence section.
- The Evidence section lists the lowest-defect cells as (2,1,2), (3,2,3), (4,3,4), (5,3,4). The table also has (3,1,2), (4,1,2) and (5,1,2), all defect 1. Harmless, but the list is not the full defect-1 set.
- The control members have defect 2 only, not defect 1. The author says so, but the title's "passes" rests on members that are not the ones at the wall. No defect-1 member (such as (5,3,4) or (5,1,2)) was ever used as a control source, and these come from other LNAs.
- The control is round-trip by construction: each source member is forward-reached from 3033030, so finding 3033030 again tests search completeness on reachable members only. It cannot show the search would flag a defect <= 0 member, because none is known. The author's "A refutation would be..." sentence concedes this.
- `maverick_cordcontrol_L6.txt` is a header line only: "58 members, 14 at L = 6 with cords >= 1, 173.5 s". No L = 6 control was run. The Evidence section's "L = 6 attempt" and the Next line ("see L6 file") overstate what the file holds. It holds the count of candidates, not an attempt at a result.

## New?

Mostly recorded.

- Claim (4), the census, is E-065: "none has relations <= cords", ten (cords, relations) pairs at depth 4 and 5 for all nine LNAs, minima falling with depth for 4444400. H-017's own table lists the same pairs. The author admits this ("not a new result").
- The "H-017 needs no new search" idea overlaps E-071/E-074 (search returns the source, L = 5 control, 0 at depth L-1) and E-084 (control members without cords). E-084's "its members have no cords" is the gap this round fills.
- The separation of cycle members from quipu cords is already E-103 (line 678: "A 'cord member' in E-089/E-094/E-101 is ... not the GLOSSARY quipu cord") and E-094's Limits ("differ in kind"). Cite those; claim (2) is a restatement.
- Genuinely new: the cords >= 1 control on a tree quipu (3 of 3 at L = 5, 3033030).

## Evidenced?

The control is specific: LNA, L, cords/rels, time, hit and neg. The census is also specific (member counts per LNA, plan file). The following are not stated or are thin:

- Sample size of the control: 3 members of 1 LNA out of 9 LNAs and 18 eligible members for that LNA.
- Whether the Coxeter guard was on for the back search. E-074 says it was; this report does not say.
- Which defect-1 cells exist for which LNA (the table is summed over LNAs).

## Scope

The title is wider than the evidence. "A cord-bearing positive control at n = 9, L = 5 passes (3 of 3)" is accurate. "H-017 needs no new search that sees relations" is a verdict on the search from three round trips of one LNA. Narrowed wording: "At n = 9, L = 5, `linesReachedFrom` recovers the source from 3 of 3 cord-bearing monomial tree members of 3033030 (defect 2), and not at depth 4; the depth-5 census of cords >= 1 members of the nine outside LNAs (E-065) has no cell of defect <= 0."

The "park T6" recommendation is the author's opinion, not a result. It is flagged as for chair and skeptic, which is fine.

## Required for acceptance

1. Fix "31 cells" to 32, and list all defect-1 cells (add (3,1,2), (4,1,2), (5,1,2)).
2. Change the L6 file description to "candidate count only; no L = 6 control was run". Remove the "L = 6 attempt" wording.
3. State in the Claim that the control members all have defect 2, and that no defect-1 member was used as a source. Narrow the title to the control actually run (3 members, 1 LNA, L = 5).
4. Cite E-103 (line 678) and E-065's depth 4-6 tally under Prior record for claims (2) and (4). Mark them restatements.
5. State whether the Coxeter guard was on in `--run`.
6. [next round, optional] Run the control for one defect-1 member (cells (5,1,2) or (5,3,4), other LNAs), and the L = 6 control on 3033030 (14 candidates, about 4 min per member by the author's sizing).
