# An n = 8 control with cords (8-9 arrows) and 1-2 commutativity relations finds its source LNA at depth 6 in 2 of 2 and at depth 5 in 0 of 2, at 5.7e4-6.2e4 nodes, matching the n = 9 negatives; no cord member with monomial relations exists in the walk

author: toolsmith · round: 015 · kind: tool
thread: T6 · bears on: H-017, E-076, E-081, E-082

## Claim

At n = 8 I built controls whose start has arrows >= n and relations >= 1, and ran the H-017 depth-first search (both directions, canonical-key node count) on them. Two members, taken from walks out of LNAs 000300 (index 10) and 000030 (index 4) at path length 6: member A has 9 arrows (2 independent cycles) and 2 relations; member B has 8 arrows (1 cycle) and 1 relation. Depth 6 returns the source LNA in both; depth 5 returns it in neither. Nodes (undeduped / distinct keys): A depth 6 56 886 / 6 157, depth 5 12 112 / 2 199; B depth 6 61 772 / 5 225, depth 5 12 851 / 1 990. The n = 9 depth-6 negatives are 5.0e4-6.3e4 (E-081), so these controls are the same size and now include cord-bearing starts. Both found-at-depth-6 and not-at-5, so at n = 8 the search finds a cord-bearing member at its depth and not one short; the E-076 negatives are not explained away by the walk being blind to cords.

It does not claim: (1) the starts are monomial. Every cord member the walk produces has a relation that is a SUM of two paths (commutativity squares), as in the members printed in the output; the n = 9 candidates are monomial quipus with 2-3 cords. A monomial relation with arrows >= n did not occur: with `MONO=1` the plan finds 0 members at n = 6 (13 LNAs, L = 5), and in the raw visitor over n = 7 LNAs 0-6 to depth 5 every node with arrows >= n has a non-monomial relation (14+184+246+38 nodes, none monomial). I did not run n = 8 with MONO=1; so "no monomial cord member" is checked at n = 6, 7 only, 7-13 LNAs. (2) That depth is the true distance: members are keyed by labelled algebra with the shortest path found, and in the n = 6 smoke test one member was found at depth 3 though its recorded length was 4 (a relabelled copy is nearer). Here depth L-1 gave not-found both times, so it did not bite. (3) Sample size: 2 members, 2 LNAs of 429; the walk to depth 6 is 130 s and a search 250 s per member.

## Evidence

Which LNAs have cord members (n = 8, L = 5, LNAs 0-15, `toolsmith_cords_n8_L5_plan.txt`): only 4 (000030), 9, 10, 11, 12, 13 (those with a "3" at the left of the sequence); 0 for the other ten. Members with 9 arrows exist from LNAs 10, 11, 13. Counts at L = 6: LNA 10, 736 members (8 arrows: 696, 9 arrows: 40; relations 1-6); LNA 4, 770 (all 8 arrows, relations 1-5).

| source | n | arrows | rels | search depth | found | nodes | distinct | s |
|---|---|---|---|---|---|---|---|---|
| 000300 | 8 | 9 | 2 | 6 | yes | 56 886 | 6 157 | 246 |
| 000300 | 8 | 9 | 2 | 5 | no | 12 112 | 2 199 | 51 |
| 000030 | 8 | 8 | 1 | 6 | yes | 61 772 | 5 225 | 247 |
| 000030 | 8 | 8 | 1 | 5 | no | 12 851 | 1 990 | 51 |
| 000030 | 6 | 6 | 1 | 4 | yes (2 of 2) | 573-648 | 251-266 | 1-2 |
| 000030 | 6 | 6 | 1 | 3 | no / yes (relabelled, nearer) | 167-180 | 100-107 | 0 |

Why `reachedQuipuAlgebras` (E-082) had none: `quipuRelations.certificate` returns None unless the quiver is a quipu tree with monomial relations; at LNA 0, n = 8, depth 6, all 2 168 members have 7 arrows (`toolsmith_survey_lna0.txt`). So E-082's "no member has cords" is a property of that enumerator, not of the walk. A second member from LNA 10 (9 arrows, 2 relations, rank 2 in the ordering) was cut by the 10 minute cap and is not counted.

## Reproduction

All from the repository root; timings with one job at a time on 4 cores.
```
timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 8 5 1 1 0 15 --plan   # which LNAs have cord members, ~5 min
BOTH=1 timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 8 6 1 1 10 10   # member A: ~8 min (walk 128 s, 246 s + 51 s; PER=2 overruns)
BOTH=1 timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 8 6 1 1 4 4     # member B: ~8 min
BOTH=1 timeout 8m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 6 4 1 2 4 4      # n = 6 smoke test, 5 s
MONO=1 timeout 4m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 6 5 1 1 0 12 --plan   # no monomial cord member (all 0)
timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_survey.py 8 6 0 0                # 3 min, the 7-arrows-only survey
```
Arguments: N L MINREL PER FIRST LAST; env BOTH (search at L and L-1), SHORT, HIGH, MONO. Outputs `toolsmith_cords_*.txt`, `toolsmith_survey_lna0.txt` in this directory; `toolsmith_cords_plan.txt` is the empty first attempt (the filter on `reachedQuipuAlgebras`, 0 members); the run line for member A in the output file is `... 8 6 1 2 10 10` (PER=2).

## Prior record

E-082 asked for exactly this (Next: "control member with cords > 0 and relations >= 1"); E-081, E-069, E-072, E-076 as in its Prior record. No `research/` entry on cord-bearing control members. A SUM relation in the start is not the n = 9 situation: those are monomial; so the control matches in cords, arrows, nodes and relation count but not in relation type.

## Code changed

New `workshop/rounds/015/toolsmith_cords.py` (members from a raw visitor, not `certificate`; handles sum relations) and `toolsmith_survey.py`. No library change, no pytest. Members rebuilt from (arrows, relation paths) by `algOf`; not cross-checked with a second implementation.

## Next

Chair: T6 can cite "control with cords found at depth, not one short, at the node count of the n = 9 negatives", with the monomial caveat. Skeptic/theorist: is a sum relation vs monomial relation a difference that matters for the search, and is a monomial cord algebra ever derived equivalent to a quipu-from-LNA at small n? (the n = 9 candidates are not reached from any LNA by this search, which is the question H-017 asks). Toolsmith, shard: more members (LNAs 9, 11, 13 with 9 arrows), n = 8 `MONO=1` plan at L = 6.
