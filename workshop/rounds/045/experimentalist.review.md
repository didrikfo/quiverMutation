# Review of workshop/rounds/045/experimentalist.md

referee: skeptic · round: 045
verdict: minor revision

## Reproduction

- `experimentalist_t10b.py --depth 3` (3 modes x 10 pairs): 45 s. Output identical to the saved `_d3.txt` except two per-pair wall-time fields (1 s vs 2 s). Every pair: 1 meeting, total 6, in guard, tilt and both; edge counts as in the table.
- `experimentalist_t10b_sweep.py 6 4 100`: 8 s, 3270 gate edges, all J0, 42/42 LNAs, as claimed. n = 7 and n = 8 sweeps not re-run (83 s and 250 s cap); the three sums do add to 56 692.
- Reimplementation check (the requested one): for pairs 0 and 5 I ran the library `search.meetingPoints(A, B, 3, alsoDual=True)` against the author's `reach(..., 'guard')`. Library: one meeting, paths (3, 3). Author: one meeting (3, 3). The set of meeting keys is identical (`samekeys True`). I also read `mutationSearchDepthFirst` (search.py 420-471). Its step is the gate, then the illegal-relation check on arrow relations, then `reducePathAlgebra`, then the key guard (drop if the key is not None and differs from the base key), and no descent from a cyclic quiver. The author's `children` does the same thing in the same order. The only difference is the author's per-`quiverKey` dedup. The library walk is undeduped, so the dedup could only matter if two different algebras shared a key. That cannot happen: the key is the full edge set plus the relation set, and `None` keys (parallel arrows) are not deduped. The match is on 2 of 10 pairs; the other 8 were not run through the library.

## True?

I found no error. The tilt-mode result follows from the edge table: gate = J=0 = key kept on all 2396 edges, so the three walks are the same graph, and the "meets in every mode" statement is not a separate finding. The tilt mode contributes no information that the edge identity does not. The control at depth 2 (no meeting) is consistent with `test_meeting_in_the_middle_reaches_twice_as_far`; I did not re-run it. The relation-dual side is handled the way the library handles it: dualise, walk, then key the dual of the node.

Gap: the "J" test is applied to the algebra actually mutated, which on the dual side is the opposite algebra. The author says so. So "all J = 0" on the dual side means the opposite step passes `tiltingPlus`, not that the forward step is a tilting step. It matters only if someone reads "tilting-only walk" as covering both sides.

## New?

Nothing in `research/` records a tilting-only or J = 0 replay of the F-041 n = 8 merges (grepped FINDINGS and HYPOTHESES for "tilting-only"; checked E-037, E-038, E-086, E-144, E-147). Closest records: E-086 (about 1.3e6 guard-admitted steps, none failing `tiltingPlus`, first rejection at distance 5-8) and E-144/E-139 (tilting-only meet at n = 6, a different question). The result is a direct consequence of E-086 at depth <= 4, as the author says. It is new as a recorded check, not as information.

## Evidenced?

Yes for the n = 8 table: the 10 pairs, modes, edge counts, totals and the control are all stated. Two weaker points:

1. The sweep tally is a count of edges, and "distinct nodes per BFS summed over LNAs" counts the same algebra many times across LNAs. So "56 692 edges" overstates independent evidence. The claim "no J != 0 step occurs" is unaffected.
2. The n = 8 sweep covers 244 of 429 LNAs, chosen by a time cap in enumeration order. It is a sample, and the title says so ("244 of 429").

## Scope

The title is mostly right. It says "guarded depth-4 edges at n = 6, 7, 8 are all J = 0", which is the right narrowing: the claim stops at depth 4 and the author explicitly excludes depth 5-6 merges and E-147-type steps at distance >= 5. Two wording fixes:

- "n = 8 ... all of them" in the scope line holds only for the 10 pairs left open by `lm.ALL_MOVES` with `free=False, edges=True, doubles=True`, which is the F-041 definition. State that dependence.
- "n = 6 (42/42), n = 7 (132/132)" is not "all LNAs" for the walk shape of `mergeReport`, which also runs on the dual and from non-LNA starts. The sweep is forward-only from LNAs.

## Required for acceptance

1. Add to the report the library cross-check (`meetingPoints(..., 3, alsoDual=True)` returns the same single meeting key and (3, 3) split). I checked pairs 0 and 5 only. Run all 10 in one script (about 40 s) and record it in the file.
2. Reword the "tilting-only meets all 10" headline to say it is implied by the edge identity (gate = J=0 = key kept on all 2396 edges), not an independent result.
3. State in the Claim that on the dual side J is tested on the opposite algebra, and that the forward-dual correspondence was not checked (move it from the evidence section into the claim's limits).
4. Say that the 56 692 edge count is summed over LNAs, so the same algebra is counted more than once, and that the n = 8 subset is the first 244 of 429 in enumeration order.
