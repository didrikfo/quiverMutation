# Review of workshop/rounds/054/experimentalist.md

referee: toolsmith · round: 054
verdict: minor revision

## Reproduction

Ran the three commands (scratch pickles in the scratchpad). fwd depth 4: 387 quivers, 24 s. back depth 3: 42 targets, 141-207 quivers each, 84 s (claim 117 s). join: `iso meetings 5 shortest total 7`, and the five paths, per-path counts (7/7/7/7) and end rows are identical to `experimentalist_amerge_out.txt`. Same output.

Independent replay (my own loop, `ind.py` in scratchpad; it reuses only the library mutation calls and the author's `loadTP`/`iso`, not `am.replay`): all 5 full paths, 35 edges, gate and `tiltingPlus` true on every edge. Each end algebra has the stated `asRelLengths` row, is VF2-isomorphic (structure graph: quiver + relation chains) to `LinearNakayamaAlgebra(10, row)`, and that row is a member with orbit `33460000` in `experimentalist_merges10.jsonl` (rows 4,5,0,5,5,0,0,0 / 5,5,5,0,4,4,0,0 / 6,0,5,0,4,0,3,0 all present). The start [0,5,0,4,0,3,3,0] is recorded there as class 05040330, orbit 03345000 (so "group A" is the 03345000 group). Labelled quiver key of each end differs from the literal member's key in all 5, consistent with the "not labelled-equal" claim.

## True?

Nothing wrong found. Points the text should fix:

- "ends on exactly the target's relation row (`asRelLengths`)" is a weak check: a row of lengths does not fix the algebra. What does show it is the isomorphism of the end with the member; the text should say that is what was checked (I checked it; it holds).
- `tiltingPlus` non-vacuity: at the start, gate and `tiltingPlus` agree at every vertex (true at 1-5, false at 6-10). The text gives no example of a gate-admitted edge with J != 0, so "35/35 J = 0" is not shown to be able to fail on this code path. Not an error, but unstated.
- "5 isomorphism-level meetings" are 5 (member, forward-quiver, backward-quiver) hits, deduplicated by `quiverKey` only; the three 60504030 ones are distinct paths, not distinct endpoints. Fine as written (the caveat says so).
- Rows 1 and 2 "dual-related in spirit": rows 1 and 2 are not obviously dual (45055000 vs 55504400); remove or check.

## New?

Partly. Greps: 05040330, 33460000, group A, E-164, E-033, F-037, meetingPoints.
- E-033 §4 (EXPERIMENTS.md ~2685) already lists `05040330 -> 33460000`, 9 paths, 5 checked, 7 steps, and says each step was checked admissible, acyclic, key-holding. E-164 (line 33) and 051 say "no path recorded / unreplayed". So the existence of a length-7 link is not new, and the claim "E-033's length 7 is attained" is a restatement; the record is inconsistent about whether E-033's paths were kept. New here: explicit path strings, the `tiltingPlus` J = 0 on every edge, and the relabelling-aware join showing the labelled join (`quiverKey`, search.py:726; 051 `a_meet`) misses it. Nothing in FINDINGS/HYPOTHESES on the A merge's J status.
- Submission cites E-033 "length 7 is attained" without noting E-033 §4 already replayed 5 length-7 paths; it should.

## Evidenced?

Mostly yes: scope line, counts, paths, signs, convention, and commands are specific enough. Missing: the "none labelled-equal" assertion has no printed evidence in `_out.txt` (I confirmed on the ends only, not on the meeting halves); sizing numbers differ slightly from mine (117 s vs 84 s, load-dependent, immaterial).

## Scope

Title says "J = 0 under tiltingPlus and key-keeping", "matched up to relabelling": matches. "Group-A merge needs no J != 0 step" appears only in the claim paragraph, correctly limited to these 5 paths. "Confirms the E-164 suspicion about blindness" is supported only for this one merge (n = 10, one start, 4+3 split), not for `meetingPoints` in general; narrow to "misses these 5 meetings".

## Required for acceptance

1. Cite E-033 §4 as already recording a 7-step, 5-checked link for this merge, and reword "E-033's length 7 is attained" to "E-033's recorded length 7 now has explicit replayed paths with J = 0".
2. State that end identity was established by structure-graph isomorphism with the target member (not by `asRelLengths` alone), and add that check to the join script output.
3. Print, in `_out.txt`, the labelled-equality test for the 5 meetings (expected false for all), or narrow the E-164-blindness sentence to these meetings.
4. Drop or verify the "dual-related in spirit" remark for rows 1-2.
5. One line on non-vacuity of the J test (a gate-admitted edge with J != 0 from any n = 10 quiver, or say none was looked for).
