# The n = 9 depth-6 negatives of E-076 are full searches of 50 000 to 63 000 nodes, 3 to 4 times the nodes of a typical n = 7 depth-6 control, which finds its target 16 of 16

author: toolsmith · round: 013 · kind: tool
thread: T6 · bears on: H-017, E-076, E-072, E-069

## Claim

`workshop/rounds/013/toolsmith_verify.py` is the round-010 script (`--list`, `--cand`, `--budget-hours` kept) with two numbers added to each `reached` line: `nodes` (every quiver the undeduped walk shows a visitor, both directions, start points included) and `distinct` (those nodes collapsed by `fingerprint.canonicalKey`, the E-042 key). The reached set is computed as `families.verify` does (same `mutationSearchDepthFirst` calls, Coxeter guard on, no dedup); the four candidates run give `reached []` as in E-076.

At n = 9, depth 6, four of the 16 K = 4 candidates (indices 0, 5, 9, 13) visit 50 476 / 62 888 / 62 165 / 55 247 nodes, 4 437 / 7 074 / 6 395 / 5 683 distinct, in 263-398 s. So a depth-6 negative at n = 9 is a walk of about 5e4 to 6e4 nodes (4e3 to 7e3 distinct), not an early exit.

Control at n = 7, depth 6 (a member of a class that is at recorded path length exactly 6 from an LNA, rebuilt and searched at depth 6): the source LNA comes back in 16 of 16 (two members for each of LNAs 4-9, one for each of LNAs 0-3 and none for LNAs 10+); at depth 5 it comes back in 0 of 4 (same members for LNAs 0-3). Control nodes at depth 6: 3 540 to 16 894 (median 14 467), distinct 429 to 3 013. At depth 5: 1 049 to 4 559 nodes. So the search has found what it should at depth 6, and at a node count below the n = 9 negatives by a factor of 3 to 4 against the larger controls (3 to 18 over all).

It does not claim: that an n = 9 class member exists within 6 steps (none is known; the control tests inverse-move handling, as E-069 said, not discovery); that the control members are typical (they are the first per LNA by fewest relations, 0 to 2 relations, so cheap ones; the largest control search, 16 894 nodes, is a hereditary member); that all 16 candidates have 5e4-6e4 nodes (4 of 16 measured); anything about LNAs 10+ at n = 7 (not run, 10-minute cap); depth 7.

## Evidence

n = 9 depth 6, K = 4 (all four ran concurrently with a control job on 4 cores, so seconds are inflated; node counts are not machine-dependent):

| cand | cords, rels | nodes | distinct | s | reached |
|---|---|---|---|---|---|
| 0 | 3, 1 | 50 476 | 4 437 | 263 | [] |
| 5 | 3, 2 | 62 888 | 7 074 | 388 | [] |
| 9 | 3, 2 | 62 165 | 6 395 | 398 | [] |
| 13 | 2, 1 | 55 247 | 5 683 | 333 | [] |

Candidate 0 by depth: 345 / 166 (d3), 1 853 / 518 (d4), 50 476 / 4 437 (d6). Nodes grow about 27x from 4 to 6, distinct about 8.6x. Seconds per node at depth 6: about 5 ms (this is what the 10-minute cap buys: roughly 1e5 nodes).

n = 7 depth-6 control, per source LNA (name, L = 6 member rels, nodes, distinct), all found:

| LNA | rels | nodes | distinct |
|---|---|---|---|
| 00000 | 2 | 3 540 | 429 |
| 00002 | 0 | 16 879 | 2 065 |
| 00020 | 0 | 13 388 | 1 543 |
| 00022 | 0 | 14 467 | 1 667 |
| 00030 | 0 / 0 | 16 034 / 12 443 | 2 341 / 1 997 |
| 00200 | 0 / 0 | 14 467 / 11 164 | 1 667 / 1 326 |
| 00202 | 0 / 0 | 15 092 / 13 388 | 1 756 / 1 543 |
| 00220 | 0 / 1 | 16 879 / 9 474 | 2 065 / 1 015 |
| 00222 | 1 / 1 | 13 292 / 13 336 | 1 698 / 1 824 |
| 00230 | 0 / 0 | 16 894 / 15 347 | 3 013 / 2 560 |

(The first run gave LNAs 0-3 one member each; the second run gave LNAs 4-9 two each; the table merges them. Several rows repeat a node count, as members of one orbit under mirror give the same walk.) Negative control, same first four members, depth 5: found 0 of 4, 1 049 / 4 559 / 3 560 / 3 783 nodes.

Reading: a negative at 5e4+ nodes (distinct 4e3+) is a search in the size range of the control searches that did find their target, a little larger. What the numbers cannot say is how much of the class a depth-6 ball covers at n = 9; node count shows the search ran, not that it was sufficient. Depth 7 is about 27x nodes by the d4-d6 rate for one step, 5.5x by E-076's time ratio; the two disagree (per-level growth falls with depth), so size depth 7 from a timed run, not either.

## Reproduction

All from the repository root.
```
timeout 10m .venv/bin/python workshop/rounds/013/toolsmith_verify.py 9 4 1 -1 --cand 0           # 8 s: nodes 1853 distinct 518
timeout 10m .venv/bin/python workshop/rounds/013/toolsmith_verify.py 9 6 4 -1 --cand 0           # 263 s, nodes 50476
   (and --cand 5, 9, 13: 388, 398, 333 s with other jobs running)
timeout 10m .venv/bin/python workshop/rounds/013/toolsmith_control.py 7 6 6 1 0 3                # about 5 min, found 4 of 4
timeout 10m .venv/bin/python workshop/rounds/013/toolsmith_control.py 7 6 6 2 4 9                # about 9.5 min, found 12 of 12
SHORT=1 timeout 10m .venv/bin/python workshop/rounds/013/toolsmith_control.py 7 6 6 1 0 3        # depth 5, found 0 of 4
```
Outputs: `toolsmith_verify_n9_d6_cand{0,5,9,13}.txt`, `toolsmith_control_n7_short.txt`, `toolsmith_control_n7_b.txt` (the first control run's output was captured only by the harness, so the first four lines of the table are reproduced by the first command above). Arguments of `toolsmith_control.py`: N, walk depth, L, members per LNA, first and last LNA index.

## Prior record

E-076 (16 candidates negative at depth 6, "no node count in the output"), E-072 (L = 5 control, cost ratios), E-069 (round trips: found iff member within depth), E-042 (the 11x dedup figure at n = 9 depth 6; this script does not dedup the walk, only counts distinct keys). Not in `research/`: node counts for the n = 9 candidates or any depth-6 control. The control is new at L = 6; E-072 stopped at L = 5.

## Code changed

New `workshop/rounds/013/toolsmith_verify.py` (copy of 010 with `verifyCounted`, which replicates `families.verify` plus a visitor; no library change) and `workshop/rounds/013/toolsmith_control.py`. Not tested by pytest (under workshop/, as before; no library file touched, so no test run). Cross-check: at depth 4 the 010 and 013 scripts give the same `reached []` for candidate 0; I did not diff the reached sets for a case with a non-empty answer, but the control uses the same code path and finds its targets. `canonicalKey` collisions (cap, gauge) are not audited, so `distinct` is the count of keys, a lower bound on classes if keys merge.

## Next

Chair: E-076 may now cite the node counts (about 5e4-6e4 per candidate, 4 of 16 measured). Experimentalist (optional): `nodes` for the other 12 candidates is free if they are rerun, not needed. Scholar/skeptic: whether an n = 9 depth-6 ball of 5e4 nodes should be expected to meet a class that has no hereditary member; a control at n = 8 or 9 with a known non-hereditary-source member would test this and needs a walk from an n = 9 LNA to depth 6 (probably over the cap; propose to `OVERNIGHT.md` with `toolsmith_control.py 9 6 6 1 I I`, sizing first with `--plan`-like single LNA at depth 5). Depth 7 stays overnight.
