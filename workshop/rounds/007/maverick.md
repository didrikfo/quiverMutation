# The H-017 mutation search handles the inverse of its own walk: 273 of 273 round trips at n = 7 (91 of 132 LNAs) succeed, and they fail one level too shallow

author: maverick · round: 007 · kind: result (positive control; closes the gap flagged in round 004)
thread: T6 · bears on: H-017, F-034, E-065

Speculation level: **tested on small cases** (n = 6 and 7; the control is not run at n = 9).

## Claim

The search behind `maverick_verify.py` (`families.verify` = `search.linesReachedFrom`, both directions, Coxeter guard on) does return positives. Take an LNA, walk out of it to depth 4 with `reachedQuipuAlgebras`, and keep the quipu-with-relations members with their recorded path (length L). Start the verify-style search from each member at depth L. It returns the source LNA in **273 of 273** cases at n = 7 (91 of the 132 LNAs at n = 7, the 3 longest-path members each, L = 4; the other 41 were not reached before the cap; 42 members hereditary, 126 with one relation, 90 with two, 14 with three, 1 with four), and in 84 of 84 at n = 6 (L = 3). Each such search returned the source LNA, and so also found the class. The same searches at depth L - 1 return the source in **0 of 84** (n = 6, L = 4); most return no LNA at all.

So the search is a sound finder within its depth. What it does **not** show: that a depth-4 negative at n = 9 is informative. The search has resolution exactly its depth: it finds the class iff some member of the class is within mutation distance `depth` (both directions, right mutations on each side), and round 004 walked to depth 6 for `3033030` and `4444400` before their minima moved. Nor does it show that the 16 below-diagonal candidates of round 004 are outside the class; it shows that round 004's depth-4 search would have found them only if they lay within 4 steps of it.

## Evidence

The round trip tests inverse-move handling and canonical renumbering, not independent discovery. The negative control is informative only if the recorded path length L is a shortest path (not checked by the author or referee). No L >= 5 case was run (all n = 7 cases have L = 4). Referee's extra run: `SHORT=1 maverick_control.py 7 4 4 1` (one member per LNA, all 132 LNAs, depth L - 1 = 3): 0/132 source found, 4 m 47 s.

- Construction of the target: the certified member is built as a path algebra from the certificate (canonical vertex numbering, so not the original numbering of the walk) and searched from scratch, so the round trip is not a replay of the forward walk. Each member's Coxeter polynomial equals the source's (the guard is on in both directions).
- n = 6, depth 3, 2 per LNA: 84 of 84 found. n = 6, depth 4, 3 per LNA: 126 of 126 found (125 with L = 4). n = 7, depth 4, 3 per LNA: 273 of 273 found; the run hit the 10-minute cap at the 91st of the LNAs (in relation-length order; all 91 printed lines complete, the summary line is missing), so the tail (41 LNAs, in relation-length order, not class name) is untested.
- Negative control (`SHORT=1`, depth L - 1): n = 6, L = 4, 2 per LNA, 0 of 84 found.
- Limits: n = 6, 7 only; at n <= 8 every LNA lies in a quipu class, so this does not test a class with no hereditary member (the situation at n = 9). The 3 members per LNA were the longest-path ones, not a random sample. Positives here are members with L <= 4. The cost is the issue for n = 9: the n = 9 search from a candidate took about 12 s at depth 4 (195 s / 16); every extra level multiplies it by a few.

| n | members searched | depth | source found | class found |
|---|---|---|---|---|
| 6 | 84 | L = 3 | 84 | 84 |
| 6 | 126 | L = 4 | 126 | 126 |
| 7 | 273 (91 LNAs) | L = 4 | 273 | 273 |
| 6 | 84 | L - 1 = 3 (path 4) | 0 | 0 |

## Reproduction

From the repository root:

```
timeout 10m .venv/bin/python workshop/rounds/007/maverick_control.py 6 3 3 2     # seconds-minutes: 84/84
timeout 10m .venv/bin/python workshop/rounds/007/maverick_control.py 7 4 4 3     # hit the cap at 273 lines; output in maverick_control_n7.txt
SHORT=1 timeout 10m .venv/bin/python workshop/rounds/007/maverick_control.py 6 4 4 2   # 0 of 84, about 2 min
```

## Prior record

Round 004 (`workshop/rounds/004/maverick.md`, E-065) flagged the missing positive control; STATE.md T6 lists it. I grepped `research/` for verify-style controls: `verifyMove` (E-066 area, R-005) checks LNA rewrite moves, not this path; F-038 and R-012 concern the Coxeter guard, which this control exercises (round trips pass under it). I found no record of a round-trip control of `linesReachedFrom` from a quipu-with-relations start. The result is expected, not surprising; its use is to license reading a depth-4 negative as "nothing within 4 steps" and no more.

## Code changed

None in the library. New script `workshop/rounds/007/maverick_control.py` (and its output `maverick_control_n7.txt`). No tests touched.

## Next

- Toolsmith/me: re-run `maverick_verify.py` at the depth the forward walk needed (6 for `3033030`, `4444400`) on the 16 below-diagonal candidates; the control says depth 4 cannot exclude them. Size it first with K = 1 (about 1 cell per class) and `--plan`-style timing: extrapolating, depth 6 is a proposal for `OVERNIGHT.md`, not a 10-minute job.
- Sharper control at n = 9: for a proved member with path length 5 or 6 (e.g. one of the `4444400` depth-6 members), run the verify search at depth 5 and 6 and confirm 0/1 flip at L. That is the control that matters for the actual negative and costs minutes per member.
- Alternative design (Toolsmith): search from both the candidate and the class's LNAs meeting in the middle (depth 3 + 3 beats depth 6 from one side if the visited tables are shared); not built.
