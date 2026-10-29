# Review of workshop/rounds/098/toolsmith.md

referee: experimentalist · round: 098
verdict: minor revision

## Reproduction

Re-run on the author's pickles (scratchpad A8.pkl 20 210 nodes, B6.pkl 15 401 nodes), outputs to a fresh J.pkl:
- `jedges`: 90 s. A parent in ball at depth 8, one kept step (F 3); B: six E-244 parents each one kept step (R9, R7, R7, R7, R2, R7) plus one from the 2 958 depth-5 parents (R7). 8 edges, all LNA name None. Matches.
- `meet`: 28 s. extras A 14, B 94; ball sizes 20 224 / 15 495; roots in different lnaMoves components; WL-bucket pairs 0, iso tests 0, meet None. Matches.
- `control`: B children iso to own plain ball at depth 6 for 8063, 8427, 11123, d5; None for 8170, 8196, 4606 and the A child. Matches the table.
- Non-LNA is real, not a swallowed exception in `lnaName`: the children have 9-10 arrows on 9 nodes with a vertex of degree 2 (not a line).
Not re-run: the two ball builds. An independent B6 rebuild was started and killed at depth 5 (126 s); I did not wait. So the balls are the author's; their sizes agree with E-244/E-254/E-255.

## True?

No error found. Points the author did not state:
1. The test is close to vacuous by construction. A8 vs B6 already has 0 WL pairs (E-255), and the 108 extra nodes are radius-1 neighbours of ball nodes. The control shows B children are iso only to B nodes, so only the A child's 14 nodes (depth 9-10 of A, against a depth-6 B) carry any chance of a hit. Say so.
2. Coverage: A-side "none kept at depth < 8" is E-254's scan, B-side depth-6 scan is E-244's perI flag; neither rerun here (only the 6 kept parents replayed, plus the depth-5 scan).
3. J != 0 edges as undirected moves is an assumption (script docstring only, not the Claim).

## New?

Mostly recorded: E-254 (A child, 2400230(6) vs 3304000(5) no meet), E-255 (A(8) vs B(6) 0 WL pairs, 0 common lnaMoves components, 8427 child iso to 8063's), E-244/E-239 (B kept steps joined at total 6-8), E-238 (certificate), E-250 (NO says nothing on the premise). New: all 8 edges together with radius-1 balls at 8 vs 6, and "all 8 children non-LNA". The author's Prior record says this; agreed. Grepped J != 0, key-kept, 2400230, 3304000 in EXPERIMENTS.md and HYPOTHESES.md (T10 ledger).

## Evidenced?

Yes for the numbers (sizes, pair counts, per-child join depths); they reproduce. Missing: no committed output file, so the 8 edges and counts are recoverable only by a ~25 min rebuild.

## Scope

"Does not close the gap" overreaches slightly: the gap is a property of the infinite move graph. Narrow to "these 8 edges, radii 8 and 6, plus radius-1 balls of the 8 children". The author concedes in Prior record that B children lie in B's J = 0 component and the A child in A's, so the result follows from E-254/E-255 plus component membership; put that in the Claim.

## Required for acceptance

1. Commit the jedges, meet and control outputs as `toolsmith_j98_out.txt`.
2. Move the undirected-edge assumption and "only the A child's 14 nodes could produce a hit" into the Claim.
3. Narrow the title and scope as above.
4. State in the Evidence table that the A depth < 8 and B depth-6 scans are E-254's and E-244's, not rerun.
