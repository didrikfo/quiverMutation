# At n = 7 classes 1 and 2, in the explored key-guarded graph (undirected shadow, 20 000 expansions), 11 of E-149's 25 failing key-keepers are pendant leaves and the other 14 lie only on seed-to-seed walks of length 15 or more (as expected from parent depth 7-8; shortest merge is 4)

author: experimentalist · round: 049 · kind: result
thread: T10 (ii) · bears on: E-149, E-153, E-152, H-015
scope: n = 7, key classes 1 and 2 (12 and 14 seeds = LNAs + duals), key-guarded BFS capped at 20 000 expansions each (depth 8 closed, depth 9 partly: pos 7 967 of 17 482 in c1, 3 216 of 21 644 in c2); frontier not closed (unexplored at depth 9: 54 percent in c1, 85 percent in c2); one walk per class. Explored graph only; the graph is treated as undirected, so "path" and "walk" are statements about the undirected shadow, not about tilting-reachable merge paths.

## Response to referee

1. Done. "45 percent" was wrong; the unexplored part of depth 9 is 54 percent (c1: 7 967 of 17 482 done) and 85 percent (c2: 3 216 of 21 644 done). Fixed in Claim.
2. Accepted. A walk through a depth 7-8 edge costs about d(u)+1+d(w) >= 14 by BFS construction, so the walk-length result (15/16/18 against a shortest merge of 4) follows from E-149's parent depths and is not independent evidence. The substantive content is the pendant/in-block split; title and Claim reworded accordingly.
3. Done. The contracted graph is `nx.Graph`, undirected. Tilting steps inside the gate need not be reversible, so "lies in the block of s" and "walk length" describe the undirected shadow; they bound what a merge path could use, they do not exhibit one.
4. Done. E-149's 80 978 / 79 143 are not this graph's 80 922 / 79 097: the differences, 56 and 46, equal the "tiltingPlus failures that move the key" row. So E-149's denominators count those key-moving failures as well; its label "key-keeping" is loose for the denominator (the 16 / 9 numerators are key-keeping and match). Here "key-keeping" means the child key equals the base key, and the edge counts exclude the key-moving ones. This is a reconciliation by arithmetic, not a match.
5. Done (`experimentalist_pendant.py`, outputs `experimentalist_pendant_c{1,2}.txt`, on the same checkpoints). c2: all 9 pendant children are unexpanded (still in the queue, depth 9 beyond pos or in the next level); none is expanded-with-no-child. c1: of 16 failing children 13 are unexpanded (2 of them in-degree 1 = the pendant ones, 11 with a second parent), 3 are expanded and have 4 admitted children each (2 at depth 8, 1 at an earlier level). So "dead ends" is dropped: pendant means only "not yet expanded", and the pendant/in-block split is largely a cap artifact (children at depth 8 or 9 that the walk did not reach). The earlier claim that c2's failures are "without an admitted continuation" is withdrawn.
6. `[next round]` Not done this sitting: committing the edge lists (tens of MB) is outside my file limits; cheaper reproduction (smaller cap) would not reach the depth 7-8 failures. The checkpoints remain in /tmp/e49 (will vanish); rerun is the two slices below.
Also: the 56 / 46 key-moving failures are not followed, unchanged. The `merged` flag in the replay is still not recomputed by the blocks script (point (a) of Evidenced) and not done here.

## Claim

The 16 (c1) and 9 (c2) key-keeping `tiltingPlus`-false edges of E-149 reproduce exactly on a walk of the same size. Contract the 12 / 14 seeds to one node s and take the explored key-kept graph (80 922 and 79 097 edges).
- c2: all 9 failing edges are pendant: the child has in-degree 1 and no other edge (in all 9 the child is unexpanded); no failing edge is in the biconnected block of s, so none lies on any path between two LNAs. All 9 are tree edges.
- c1: 7 of 16 are merge edges (child already seen), 9 tree. 14 of 16 lie in the block of s, so on some simple seed-to-seed path, but the shortest seed-to-seed walk through any of them has length 15 to 18 (parent depth 7-8 plus the child's other route). The 2 others (depth 8 to 9) are pendant, as in c2.
- The shortest merge walk anywhere in the same graph has length 4 (88 edges in c1, 136 in c2; next 6, 8, 10). So no failing edge lies on a shortest merge path, and none on a merge walk shorter than 15 in this graph; this follows from parent depth 7-8 and is not independent of E-149.
Not claimed: that failing edges lie on no merge path at larger depth or in the unexplored depth 9 (54 percent in c1, 85 percent in c2); that pendant children stay pendant once expanded (13 of 25 are unexpanded); that they are harmless for the derived class (E-152 left that open); that the walk length 15 is a minimum over all of the class (explored graph gives an upper bound on distances, so the true shortest path through a failing edge could be shorter if the unexplored graph helps). Refuted by: a failing edge on a seed-to-seed walk of length 4 to 10 in a longer walk.

## Evidence

| | c1 | c2 |
|---|---|---|
| seeds / nodes / edges | 12 / 48 394 / 80 922 | 14 / 45 427 / 79 097 |
| failing key-keepers | 16 (7 merge, 9 tree) | 9 (0 merge, 9 tree) |
| parent depth | 7 (5), 8 (11) | 7 (6), 8 (3) |
| in block of s | 14 | 0 |
| pendant (in-degree 1, degree 1) | 2 | 9 |
| walk length through e (in block) | 15 (1), 16 (4), 18 (9) | none |
| shortest merge walk, any edge | 4 | 4 |
| `tiltingPlus` failures that move the key (not followed) | 56 | 46 |

The failing counts match E-149 (16 and 9); its denominators 80 978 / 79 143 exceed these edge counts by the key-moving failures 56 / 46. And the `tiltingPlus`-failing edges inside the key class come with a clean parent in every case (0 tainted nodes before the first failure; c1 has 2 tainted-parent tilt edges and 2 tainted tree edges after the failures, all in the depth 9 slice). Walk length means d_{G-e}(s,u) + 1 + d_{G-e}(s,w) with e removed, in the contracted explored graph (a closed walk through s, so an upper bound on the simple path length). Block test: networkx biconnected components of the simple contracted graph; edge multiplicity was 1 for all 25 failing edges.

Reading: in c2 the failing children are unexpanded leaves (cap artifact, not dead ends); in c1 the failing edge into w with a second parent (e.g. depth 8 to 9, in-degree 2) joins the block only because w has been reached by another route at depth 9, so the cycle through e is as long as 2 x 8. A shortest merge between LNAs is of depth 2 + 2 at n = 7 here and uses no failing edge. This is the n = 7 positive control for E-153: the tally can see failures (it finds all 25), and it finds them off the merge paths. So at n = 7 and n = 8 c2 the claim "no recorded merge needs a J != 0 step" is not contradicted; it is not proved either.

## Reproduction

Each of the four slices below used `--budget-sec 400/450`; about 1 400 s of wall clock in all, two classes in parallel.

```
for c in 1 2; do timeout 10m .venv/bin/python -u workshop/rounds/049/experimentalist_deepreplay.py 7 --class $c --max-exp 20000 --budget-sec 400 --ckpt /tmp/e49/c$c.ckpt; done   # slice 1, exit 2 = checkpoint
for c in 1 2; do timeout 10m .venv/bin/python -u workshop/rounds/049/experimentalist_deepreplay.py 7 --class $c --max-exp 20000 --budget-sec 450 --ckpt /tmp/e49/c$c.ckpt; done   # slice 2, resumes
for c in 1 2; do timeout 5m .venv/bin/python workshop/rounds/049/experimentalist_blocks.py /tmp/e49/c$c.ckpt; done   # a few seconds
```
Outputs: `workshop/rounds/049/experimentalist_deepreplay_c{1,2}_out.txt` (slice 1 then slice 2), `experimentalist_blocks_c{1,2}.txt`. The checkpoints (with the edge lists, tens of MB) are under /tmp and not committed.

## Prior record

E-149 (counts, parent depths), E-153 (n = 8 c2 tally: 0 failing key-keepers, and it named this control as undone), E-148 / E-154 (n = 8 merges are 3 + 3, all J = 0), E-152 (Cartan data cannot separate the 25 children). The edge-level location of E-149's failures relative to merge paths is new; nothing in RETRACTIONS.md bears on it.

## Code changed

New `workshop/rounds/049/experimentalist_deepreplay.py` (round 047 script plus an edge list (parent key, child key, fail) in the checkpoint; queue items carry the parent key) and `experimentalist_blocks.py` (analysis). No library change, no tests needed. Caveat: the walk uses canonical-key dedup, and a node with no canonical key gets a fresh id, so such nodes (a few at depth 9) can never merge.

## Next

- skeptic: freeze the list of 25 failing edges (the F lines in the _out files give path and vertex) and re-run slice 2 for c1 only; check one walk-length-15 case by hand (the depth 7 to 6 edge, its child seen at lower depth).
- toolsmith (T10 i): the 14 c1 children in the block with a second parent are a ready test for tilting-only paths back to an LNA; the c2 pendant children are not on the merge graph, so a path back to an LNA from them would be the E-152 question in its pure form.
- experimentalist: finish depth 9 (c1 and c2 about 10 and 18 minutes each more) and n = 8 depth 9 overnight if wanted; the tally is now cheap to rerun.
