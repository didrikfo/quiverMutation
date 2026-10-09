# Of the 25 key-keeping, tiltingPlus-failing E-152 children (n = 7 classes 1, 2), 19 have a tilting-only (J = 0) path back to an LNA of their class (length 6-11); 5 miss at total length 11 and 1 is undecided

author: toolsmith · round: 049 · kind: result
thread: T10 (i) · bears on: H-015, E-145, E-149, E-152
scope: n = 7, key classes 1 (16 children) and 2 (9 children), the E-149 walk rebuilt (20 000 expansions each); meet-in-the-middle, target ball depth 5, child ball depth 6 (child 12 of c1: see below); caps never reached (400 000 nodes); controls: random forward J = 0 walks of length 9, 6 per class. n = 8 and E-145's own step list not covered.

## Claim

Each of the 25 children (16 distinct steps in c1, 9 in c2; 10 + 9 distinct canonical keys) was searched with tilting-only moves, F (forward step at v) and R (the inverse of a forward step, taken as the forward step at v of the opposite algebra), every move gate-admitted, legal, vertex-preserving, J = 0, `tiltingPlus` true and key-kept. A shared canonical key with the ball of depth 5 around the class's LNAs and duals is a path of length d(child) + d(target). Result: **19 of 25 children reach an LNA/dual of their own class** (c2 8 of 9, c1 11 of 16), total lengths 6-11. 5 miss at bound 6 + 5 = 11 (c2 child 6; c1 children 5, 13, 14, 15; 3 distinct keys). c1 child 12 is **undecided**: complete to depth 3 (no hit), depth 4 and 5 did not finish in the time limit (parallel arrows make a node about 60 ms).

What this says: for a hit, if F and R steps with `tiltingPlus` true are derived equivalences, the child is derived equivalent to an LNA, so **those 19 children (15 distinct keys) are in the class**, and the key-guard's failure on them is not a counter-example to "the walk stays in one class" -- only to "key kept implies tilting". What it does not say: nothing about the 5 misses (the search is a bounded one; a miss is not an out-of-class proof) and nothing about n >= 8. Refuted if a hit path fails replay under an independent tilting definition (not done, see Next). The result depends on R as the inverse of a tilting step (E-147/witness convention of round 047); an F-only search is weaker (below).

## Evidence

Positive control, same machinery and same depth: random forward J = 0 walks of length 9 from the class's own LNAs, then the identical search: **c2 6/6 hits, c1 6/6 hits**, totals 6-9 (a known path of length 9 is found in every case, sometimes a shorter one; child-ball depth 4, target depth 5). Controls of length 8 (c2, 3/3) and 4 (3/3) likewise. So the search has demonstrated power at the children's depth (parent depth 7-8, child at forward distance 8-9).

| class | children (steps) | hit | miss at 6 + 5 | undecided | hit totals |
|---|---|---|---|---|---|
| c2 | 9 | 8 | 1 (child 6, parent depth 8 v 5) | 0 | 8, 8, 9, 10, 11 x 4 |
| c1 | 16 | 11 | 4 (children 5, 13: key f7abe9; 14, 15: key b32eca; all v = 7, depth 8) | 1 (child 12, v 4 depth 8) | 6, 8 x 7, 9 x 3 |

Miss bound: child ball depth 6 = 5 703 distinct nodes (c2 child 6), 18 034 (c1 children 5, 13), 18 002 (14, 15), uncapped; target ball depth 5 = 25 452 keys (c2), 15 862 (c1) (seeds 14, 12 LNAs + duals). A miss means no path of length <= 11 of this kind (each side exhaustive to its depth, meeting on exactly equal canonical keys); the shortest hit is 6, the longest 11, so 11 is not a margin. 4 of 5 misses are c1 v = 7 steps; hit totals for hit children shorter than 11 mostly 8-9 = the forward depth.

Why the skeptic's round-046 search (E-152: 9 c2 children, 13 misses incl. controls, no power) missed: it expanded forward moves only from the child (6000 nodes) towards a single-sided target, while the control child at forward distance 4 should have an LNA 4 steps back. Here the inverse moves R and a meet with a depth-5 target ball are what bring the controls (and the c2 children) in; with F only the reverse loses about 10% of edges (E-147). I did not run an F-only variant here, so "R is what matters" is a reading, not tested.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/049/toolsmith_collect.py 7 2 20000 /tmp/tsm/c2.pkl 100      # 9 children; 580 s loaded (c1: 540 s, run it alone)
timeout 10m .venv/bin/python workshop/rounds/049/toolsmith_tiltpath.py /tmp/tsm/c2.pkl 2 5 4 400000 ctrl6:9    # 413 s: ball 5, control 6/6
SLICE=0:9 timeout 10m .venv/bin/python workshop/rounds/049/toolsmith_tiltpath.py /tmp/tsm/c2.pkl 2 5 6 400000 children   # ball cached; ~10 min at dC 6 (9 children, 5 of them 85 s)
```
(c1: the same with `7 1` and cls 1, children in slices SLICE=0:4 etc.; the dC = 6 misses take 380-460 s each.) The ball cache `/tmp/tsm/tball_c*_d5.pkl` is scratch. Per-child log lines: `workshop/rounds/049/toolsmith_tiltpath_logs.txt` (9 KB).

## Prior record

E-152 (invariants non-informative; skeptic's back-search no power, 13 misses) and E-149 (the 25 children), E-145. New: a back-search with demonstrated power at the same depth, and 19 hits. E-147 (reverse search loses 10.3% of edges) explains why a forward-only search is weak; the hits do not contradict anything recorded. I grepped `research/` for "tilting path back", "meet" with the E-152 children: nothing.

## Code changed

New scripts only: `toolsmith_collect.py` (copy of 046's collector that also stores the algebra objects, so the 8 parallel-arrow c1 parents the hand constructor rejects are covered), `toolsmith_tiltpath.py`. No library file touched, no tests run.

## Next

- skeptic: replay one hit path per distinct hit key edge by edge with an independent tilting test (the path itself is not printed; add path recovery -- toolsmith) and test whether R really is a tilting step (inverse of the recorded F step) on those edges. Without that, "in the class" rests on F/R steps being derived equivalences, which is the premise T10 questions.
- theorist: the 5 misses (c1 v = 7 steps at parent depth 8; c2 child 6) are the candidates for "outside the class"; a hit-free depth 6 + 5 is the bound. Also c1 child 12 undecided.
- toolsmith: path recovery and a depth-7 child ball for the misses (about 70 000 nodes each, 30-40 min: OVERNIGHT proposal, not run).
