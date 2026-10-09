# Of the 25 key-keeping, tiltingPlus-failing E-154 children (n = 7 classes 1, 2), 19 are connected to an LNA of their class by a printed path of J = 0 `tiltingPlus` steps (length 6-11, F and R moves, premise: these steps are derived equivalences); 5 miss at total length 11 and 1 is undecided; their parents are connected too

author: toolsmith · round: 049 · kind: result
thread: T10 (i) · bears on: H-015, E-147, E-151, E-154
scope: n = 7, key classes 1 (16 children) and 2 (9 children), the E-151 walk rebuilt (20 000 expansions each); meet-in-the-middle, target ball depth 5, child ball depth 6 (child 12 of c1: see below); caps never reached (400 000 nodes); controls: random forward J = 0 walks of length 9 (6 per class) and 11 (3, c2). Statement is conditional: derived equivalence is the `tiltingPlus` + J = 0 criterion, not tested independently. n = 8 and E-147's own step list not covered.

## Response to referee

1. **R reworded; premise stated.** R is not the inverse of a recorded F step. It is the forward step at v of the OPPOSITE algebra, dualized back (`moves`, kind R); it is a tilting step whenever F is, because D(A) ~ D(B) iff D(A^op) ~ D(B^op) (op-duality; independent of E-149's inverse convention). The only premise is: J = 0 plus `tiltingPlus` true implies derived equivalence. No independent derived-equivalence test exists in the repo and none was run; the referee's replay re-tests the same predicates plus the necessary Cartan congruence. "In the class" below means "derived equivalent to an LNA, given that premise". Claim and Next are corrected accordingly.
2. **Parents: checked directly, and the sentence cut.** I do not claim anything from the walk's ancestry (the E-151 walk keeps J != 0 children in the frontier, so a depth-8 parent may descend from a failing step). Instead each parent was searched itself (`toolsmith_paths.py ... parents`, dC 4 against the depth-5 target ball): **all 16 distinct c1 parents and all 9 c2 parents** (including the parents of the 5 + 1 children that miss) have a J = 0 `tiltingPlus` path to an LNA/dual of the class, total 7 or 8, every replay ok (`toolsmith_paths_logs.txt`). So the parents are in the class by their own paths, and the 19 hits say the children are in the class too; the failing step is a step between two class members. The sentence "not a counter-example to the walk stays in one class" is dropped: what the hits show is that the children are in the class, nothing about the walk beyond that.
3. **Filters.** In every search `stats` has tp = 0 and key = 0 (also in the new runs): with J = 0 the `tiltingPlus` and key filters never fired. The search graph is gate + legal + vertex-preserving + J = 0 (the J = 0 filter does fire, 56-569 times per run). The list "J = 0, tiltingPlus, key kept" in the Claim is a description of what is checked, not three active filters.
4. **L = 11 control and F-only variant (both run).** Control: `toolsmith_tiltpath.py /tmp/tsm/c2.pkl 2 5 6 400000 ctrl3:11`: 3/3 hits (totals 6, 11, 10), as the referee found. F-only (env `FONLY=1`, F moves on both sides, own target ball 3 119 keys instead of 25 452): the 9 c2 children, dC 6: **0/9 hits**; but the F-only control at L = 9 (dC 4) is 1/6, against 6/6 with R. So the F-only search has no power at this depth and its 0/9 says nothing about the children; the explanation "R is what separates this from round 046" is now tested in the sense that the same machinery without R fails its own positive control (the two searches differ in other ways too: one-sided target, 6000 nodes).
5. **Paths printed** (cached balls were still in /tmp/tsm, so no 10-minute rebuild was needed beyond a ~5 min ball with parent pointers per class): one hit per distinct hit key, **15 paths (7 in c1, 8 in c2)**, in `toolsmith_paths_logs.txt`, produced by `toolsmith_paths.py`. Each is a child-side move list (child -> meeting algebra) and an LNA-side move list (LNA/dual -> meeting algebra, traversed backwards, i.e. inverse tilting steps); both lists were replayed with fresh move generation and the canonical keys compared (all "replay ok"). The same replay is available for the 25 parents. Example (c2 child 2, total 8): child F2 F7 F4 meets LNA/dual #10 after R6 R5 R4 F5 R7. Not done: a depth-7 child ball for the 5 misses and the undecided c1 child 12 (overnight, as proposed); the referee's Cartan-congruence check on my printed paths.

## Claim

Each of the 25 children (16 distinct steps in c1, 9 in c2; 10 + 9 distinct canonical keys) was searched with tilting-only moves, F (forward step at v) and R (the forward step at v of the opposite algebra, dualized back; tilting by op-duality), every move gate-admitted, legal, vertex-preserving, J = 0 (the only filter that fires), `tiltingPlus` true and key-kept (these two never fired). A shared canonical key with the ball of depth 5 around the class's LNAs and duals is a path of length d(child) + d(target). Result: **19 of 25 children reach an LNA/dual of their own class** (c2 8 of 9, c1 11 of 16), total lengths 6-11. 5 miss at bound 6 + 5 = 11 (c2 child 6; c1 children 5, 13, 14, 15; 3 distinct keys). c1 child 12 is **undecided**: complete to depth 3 (no hit), depth 4 and 5 did not finish in the time limit (parallel arrows make a node about 60 ms).

What this says: for a hit, if J = 0 plus `tiltingPlus` steps (F and R) are derived equivalences (the premise, not independently tested), the child is derived equivalent to an LNA, so **those 19 children (15 distinct keys) are in the class**, and the key-guard's failure on them is only a failure of "key kept implies tilting". The parents are in the class by their own paths (Response 2). What it does not say: nothing about the 5 misses (a bounded search; a miss is not an out-of-class proof) and nothing about n >= 8. Refuted if a printed hit path fails an independent tilting test (the referee's replay: 49 edges pass the predicates and the Cartan congruence, not an independent derived-equivalence test).

## Evidence

Positive control, same machinery and same depth: random forward J = 0 walks of length 9 from the class's own LNAs, then the identical search: **c2 6/6 hits, c1 6/6 hits**, totals 6-9 (a known path of length 9 is found in every case, sometimes a shorter one; child-ball depth 4, target depth 5). Controls of length 8 (c2, 3/3) and 4 (3/3) likewise. So the search has demonstrated power at the children's depth (parent depth 7-8, child at forward distance 8-9).

| class | children (steps) | hit | miss at 6 + 5 | undecided | hit totals |
|---|---|---|---|---|---|
| c2 | 9 | 8 | 1 (child 6, parent depth 8 v 5) | 0 | 8, 8, 9, 10, 11 x 4 |
| c1 | 16 | 11 | 4 (children 5, 13: key f7abe9; 14, 15: key b32eca; all v = 7, depth 8) | 1 (child 12, v 4 depth 8) | 6, 8 x 7, 9 x 3 |

Miss bound: child ball depth 6 = 5 703 distinct nodes (c2 child 6), 18 034 (c1 children 5, 13), 18 002 (14, 15), uncapped; target ball depth 5 = 25 452 keys (c2), 15 862 (c1) (seeds 14, 12 LNAs + duals). A miss means no path of length <= 11 of this kind (each side exhaustive to its depth, meeting on exactly equal canonical keys); the shortest hit is 6, the longest 11, so 11 is not a margin. 4 of 5 misses are c1 v = 7 steps; hit totals for hit children shorter than 11 mostly 8-9 = the forward depth.

Why the skeptic's round-046 search (E-154: 9 c2 children, 13 misses incl. controls, no power) missed: it expanded forward moves only from the child (6000 nodes) towards a single-sided target. The F-only variant of this machinery (Response 4) fails its own control (1/6 at L = 9, children 0/9), while with R the controls are 6/6 (L = 9) and 3/3 (L = 11). Consistent with "R plus the meet is what gives power"; not a proof that this was the only difference.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/049/toolsmith_collect.py 7 2 20000 /tmp/tsm/c2.pkl 100      # 9 children; 580 s loaded (c1: 540 s, run it alone)
timeout 10m .venv/bin/python workshop/rounds/049/toolsmith_tiltpath.py /tmp/tsm/c2.pkl 2 5 4 400000 ctrl6:9    # 413 s: ball 5, control 6/6
SLICE=0:9 timeout 10m .venv/bin/python workshop/rounds/049/toolsmith_tiltpath.py /tmp/tsm/c2.pkl 2 5 6 400000 children   # ball cached; ~10 min at dC 6 (9 children, 5 of them 85 s)
```
(c1: the same with `7 1` and cls 1, children in slices SLICE=0:4 etc.; the dC = 6 misses take 380-460 s each.) The ball cache `/tmp/tsm/tball_c*_d5.pkl` is scratch. Per-child log lines: `workshop/rounds/049/toolsmith_tiltpath_logs.txt` (9 KB).

## Prior record

E-154 (invariants non-informative; skeptic's back-search no power, 13 misses) and E-151 (the 25 children), E-147. New: a back-search with demonstrated power at the same depth, and 19 hits. E-149 (reverse search loses 10.3% of edges) explains why a forward-only search is weak; the hits do not contradict anything recorded. I grepped `research/` for "tilting path back", "meet" with the E-154 children: nothing.

## Code changed

New scripts only: `toolsmith_paths.py` (round 049 response: parent-pointer ball, paths, parents) and `toolsmith_paths_logs.txt`; `toolsmith_tiltpath.py` got the `FONLY=1` switch; `toolsmith_collect.py` (copy of 046's collector that also stores the algebra objects, so the 8 parallel-arrow c1 parents the hand constructor rejects are covered), `toolsmith_tiltpath.py`. No library file touched, no tests run.

## Next

- skeptic: the 15 printed paths (`toolsmith_paths_logs.txt`) can be replayed without a ball build; the open point is an independent derived-equivalence test of a J = 0 `tiltingPlus` step (F or R, R valid by op-duality). Without that, "in the class" rests on the premise T10 questions.
- theorist: the 5 misses (c1 v = 7 steps at parent depth 8; c2 child 6) are the candidates for "outside the class"; a hit-free depth 6 + 5 is the bound. Also c1 child 12 undecided.
- toolsmith: (path recovery done) a depth-7 child ball for the misses (about 70 000 nodes each, 30-40 min: OVERNIGHT proposal, not run).
