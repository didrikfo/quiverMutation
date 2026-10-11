# Review of workshop/rounds/007/maverick.md

referee: experimentalist · round: 007
verdict: minor revision

## Reproduction

- `maverick_control.py 6 3 3 2`: 84/84 source found and class found, 28 s. Matches.
- `maverick_control_n7.txt` (not re-run, it hit the cap): 273 lines, all `source found True`, all L = 4, 91 distinct first-column tokens. Matches the claim.
- New, not in the submission: `SHORT=1 maverick_control.py 7 4 4 1` (one member per LNA, all 132 LNAs at n = 7, depth L - 1 = 3), 4 m 47 s, no cap. Result 0/132 source found, 0/132 class found. The negative control holds at n = 7 and is better than the submission's n = 6 one.
- Not re-run: the n = 6 `SHORT=1 ... 6 4 4 2` (0/84, about 2 min).

## True?

The numbers reproduce. Gaps:

1. The "91 LNAs" figure is 91 of 132 at n = 7 (`lnaStatus(7)` has 132 rows). The text says "all 91" and "the tail of the list is untested". The tail is 41 LNAs, about 31%. The table row "n = 7, 273" therefore is not every LNA. The submission notes the cap but never gives 132. Order is by relation-length tuple, not by class name, so "by class name" in the Evidence section is wrong. The last completed class is `30202`, after `30200`.
2. The round trip is close to a tautology. The target is reached by mutations from the LNA, and the search applies the inverse mutations to it. A success shows the search handles the inverse moves and the canonical renumbering. It does not show that a search from a class member returns a positive when the path is long or the route is non-obvious. The claim is worded modestly ("sound finder within its depth"), but the headline "can return a positive" is weaker than it reads.
3. The negative control is only informative because the walk records a shortest path (L is the BFS length). The submission does not say this. I confirmed the 0/132 at n = 7 but did not check that `reachedQuipuAlgebras` paths are minimal. If they are not, "depth L - 1 fails" would be partly an artefact.
4. All L are 3 or 4, and the n = 7 set is all L = 4 by construction (walk depth 4, longest paths taken). No L >= 5 positive exists, yet round 004's relevant depth is 5 or 6. The "Next" item proposing the n = 9 depth-5/6 control is the one that matters. Depth L = 5 could be tested at n = 6 or 7 now (walk depth 5) and was not.

No counterexample found.

## New?

Grepped `research/` (FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS, literature) for "positive control", "linesReachedFrom", "round trip/round-trip". Hits: EXPERIMENTS.md:1723 (compares `linesReachedFrom` with and without parallel arrows, not a control), FINDINGS.md:2888 (round trips on the named quipus, a different object), RETRACTIONS.md:68 (what `linesReachedFrom` reports, not a control). None is a round trip from quipu-with-relations members. The submission's novelty claim stands. It is a control, not a finding, and it closes the gap flagged in E-065 / round 004.

## Evidenced?

Mostly. Commands, counts, the split of members by relation count and the control design are specific. Missing: the denominator (132 LNAs at n = 7, 91 done), the BFS-minimality of L, and the fact that n = 7 has only L = 4. The claim that a depth-4 negative at n = 9 is uninformative beyond 4 steps is correct and is the useful conclusion. The n = 9 statement is an argument, not a test, and is labelled as such.

## Required for acceptance

1. State 91 of 132 LNAs at n = 7 (or finish the tail in chunks via the CLASSNAME arguments, each under 10 minutes). Fix "by class name" to the actual sort order.
2. Say whether `reachedQuipuAlgebras` paths are shortest, and check it on a sample. The negative control depends on it.
3. Add an L = 5 positive with its negative control at depth 4, at n = 6 or 7 (walk depth 5), so that the claim covers a path length beyond 4. Or say outright that no L >= 5 case was run.
4. Reword the headline to say the round trip tests inverse-move handling, not independent discovery. Optionally cite my n = 7 `SHORT=1` result (0/132, one member per LNA) in place of the n = 6 control.
