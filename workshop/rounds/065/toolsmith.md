# No J != 0 key-kept edge in the 129 000 n = 8 edges (c0, c2) deep enough to bear on the n = 7 pattern; a J != 0 path between LNAs already exists implicitly at n = 7, and an n = 8 search does not fit a 10-minute command (sizing, capped samples)

author: toolsmith · round: 065 · kind: negative (sizing)
thread: T10 (T5) · bears on: E-151, E-155, E-157, E-160, E-163, E-175, E-178
scope: n = 8 key classes 0-4 (of 11), depth 6-9 (only c0 depth 9 and c2 depth 8 reach the depth where n = 7 failures begin; c1, c3, c4 and all n = 9 are too shallow to bear on it), caps by budget (not the library DFS, not whole classes); n = 9 classes 0, 1, 2, 4 (of 19), depth 5-7; classes 5-10 (n = 8) and 3, 5-18 (n = 9) not sampled; no J != 0 step was found, so nothing was run beyond sizing

## Response to referee

1. Depth column (done). Added to the table: "J != 0 first at depth" for n = 7 (c1: 7, 4 edges at 7, 12 at 8; c2: 7, 6 at 7, 3 at 8; `toolsmith_parent.py` on the saved files); n = 8, 9 rows keep "deepest full depth". Expectation recomputed on c0 + c2 only (128 945 edges): about 18, not 38. Claim, title, Scope narrowed as the referee wrote: c1, c3, c4 and n = 9 do not bear on the rate. The 1.1e-5 bound over the unmatched denominator is withdrawn; replaced by: 0 of 128 945 edges, rule-of-three 2.3e-5 per edge versus 1.4e-4 at n = 7.
2. Clustering (done). The Poisson/rule-of-three figures are per edge and failures are clustered, so they understate uncertainty. Per parent expansion: n = 7 25 per 28 817 expansions (8.7e-4); n = 8 c0 + c2 0 per 30 028 expansions (rule of three 1e-4). The 25 failing edges have 25 distinct parent paths (paths, not parent algebras: identifying parents up to isomorphism was not done), so the unit is not obviously smaller than the edge; I do not claim more than that.
3. Obstacle 2 moved into the Claim (done); Obstacle 1 stays below.
Not done: a depth-9 c0 run attempt (about 13 min by own numbers, over the 10-minute limit); stated as borderline projection.

## Claim

Existence. At n = 7 the key-guarded walks find 25 J != 0 (`tiltingPlus`-false, key-kept) edges in 174 804 key-kept edges (16 of 76 334 c1, 9 of 98 470 c2: 1.4e-4 per edge); all 25 are at parent depth 7 (10) or 8 (15), none shallower. A J != 0 path between LNAs therefore already exists implicitly at n = 7: seed A -> parent (J = 0 walk, depth 7-8) -> child (one J != 0 step) -> LNA X (J = 0 path, 6-13 steps; E-157/E-160/E-163, 25 of 25 under the premise). What is not recorded is whether X is certified equivalent to A (see Next). At n = 8 there is no failure in 269 965 edges (c0-c4) nor at n = 9 in 52 190, but only c0 (depth 9) and c2 (depth 8), 128 945 edges, reach the depth where n = 7 failures begin; there about 18 are expected at the whole-walk n = 7 rate (an upper figure, since depth <= 6 edges in the n = 7 denominator carry none), observed 0. c1, c3, c4 (141 020 edges, stop at depth 6-7) and n = 9 (depth <= 7) are not evidence on the rate. Supported: no failure in c0, c2 at depth <= 8-9; the rate there is lower or failures are deeper. Cost: an n = 8 A-side ball of depth 9 (c2: 38 907 expansions at 0.03 s = 20 min, one side, one class) plus a child leg ball of about 7 deep does not fit in `max_command_minutes` = 10; projection, not a run attempt (c0 depth 10 about 13 min single process, borderline). Not claimed: that no J != 0 edge exists at n = 8 or 9, nor anything about the J = 0 premise. A J != 0 count above 0 on c0, c2 at depth <= 9 refutes the "lower rate" half.

## Evidence

Per-class walks (`experimentalist_deepreplay.py`, E-155's tool, gate + key guard, per-edge `tiltingPlus`); "cum. exp." = cumulative expansions, time with 4 jobs on 4 cores (a loaded machine, ratios are the point).

| n, class (starts) | J != 0 first at depth / deepest full depth | cum. exp. | next level (exp.) | s/exp | key-kept edges | J != 0 |
|---|---|---|---|---|---|---|
| 7, c1 (12) | first 7 / 8 | 12 033 | 17 482 | 0.025 | 76 334 | 16 |
| 7, c2 (14) | first 7 / 8 | 16 784 | 21 644 | 0.020 | 98 470 | 9 |
| 8, c0 (2) | 9 | 5 712 | 8 975 | 0.086 | 24 316 | 0 |
| 8, c1 (8) | 7 | 5 306 | 8 463 | 0.045 | 41 428 | 0 |
| 8, c2 (18), E-155 | 8 | 24 316 | 38 907 | 0.032 | 104 629 | 0 |
| 8, c3 (20) | 7 | 11 073 | 16 077 | 0.044 | 48 575 | 0 |
| 8, c4 (26) | 6 | 5 780 | 7 534 | 0.036 | 51 017 | 0 |
| 9, c0 (2) | 7 | 1 399 | 1 986 | 0.107 | 13 073 | 0 |
| 9, c1 (2) | 7 | 1 677 | 3 196 | 0.121 | 10 269 | 0 |
| 9, c2 (12) | 6 | 2 950 | 4 147 | 0.089 | 14 354 | 0 |
| 9, c4 (16) | 5 | 2 824 | 5 125 | 0.101 | 14 494 | 0 |

Scaling used for sizing: expansions per level grow x2.4 per depth at every n; s/expansion about 0.02 (n = 7), 0.03-0.05 (n = 8), 0.1 (n = 9). Projections: n = 8 c2 depth 9 about 1 250 s (single process) and depth 10 about 50 min; n = 9 c2 depth 8 about 20 000 expansions about 35 min, depth 9 several hours. A meeting leg adds a ball of 5-7 per side (n = 7 depth-6 ball 46 435 keys, 62 297 in total, cost 9 min/shard-set in round 051), so a full n = 8 attempt is a sharded, checkpointed multi-hour job; n = 9 is out of reach overnight.

Other obstacles I found while sizing, both of the job's definition and not of its cost:
1. "Known-equivalent" endpoints. Two LNAs in one key class are not thereby equivalent (the guard is not a tilting test, E-147). A path with a J != 0 edge between LNAs A, B only tests the premise if A ~ B is certified independently (same quipu/F-047 orbit, or a J = 0 path elsewhere); otherwise it assumes what it tests.
2. A J != 0 edge by construction leaves the guarded class graph's "clean" part: its child has to be joined to an LNA by J = 0 steps (what E-157/E-160/E-163 already do, 25 of 25 under the premise), so a "J != 0 path between LNAs" at n = 7 is already implicit: A -> parent (J = 0 walk, depth 7-8) -> child (one J != 0 step) -> LNA X (J = 0 path, 6-13). What is not recorded is whether X is certified equivalent to the seed A (see Next).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/065/toolsmith_jplan.py parse workshop/rounds/065/toolsmith_jplan_n*.txt   # table above except n8c2, seconds
timeout 10m .venv/bin/python workshop/rounds/065/toolsmith_jplan.py run 8 3 500 7   # one slice, ~9 min with 3 other jobs; one class
```

Raw walks: `workshop/rounds/065/toolsmith_jplan_n{7c1,7c2,8c0,8c1,8c3,8c4,9c0,9c1,9c2,9c4}.txt` (each from `experimentalist_deepreplay.py N --class C --depth 9 --budget-sec 300..520 --ckpt <tmp>`; budget-capped, so rerunning gives the same counts for completed levels only). 8c2 is E-155's `rounds/047/experimentalist_deepreplay_out.txt`. Wall: n = 7 pair 5 min, n = 8 four classes 9 min in parallel, n = 9 four classes 5 min in parallel.

## Prior record

E-151 (n = 8 c0, c1 at 1 200 expansions: none; n = 7 failures at parent depth 7-8) and E-155 (n = 8 c2 depth 8: 0 of 104 629; depth 9 not run, "overnight, not approved") already say part of this. New here: the n = 8 tally is extended to classes 0-4 (270 k edges, not only c2), a first n = 9 sample (52 k edges, depth <= 7), the n = 7 base rate with the n = 8 zero compared on the depth-matched classes only (about 18 expected), and the cost model for the path search. The E-175/E-178 balls (n = 9, 10) contain no J != 0 step, consistent with this. Nothing here touches RETRACTIONS.

## Code changed

New `workshop/rounds/065/toolsmith_jplan.py` (wrapper and parser, no library change). No tests needed (no `quivermutation/` file touched). The parser was run on all saved outputs.

## Next

- OVERNIGHT proposal (not added to `OVERNIGHT.md`, I do not edit it): n = 8 class 2 (18 starts) guarded walk to depth 9 then 10, sharded as in round 050 (`--max-exp` / checkpoint slices of <= 8 min, ~1 250 s and ~50 min of CPU), recording J != 0 edges; run only if the human wants the existence question at n = 8 answered. If edges appear, a child leg ball of depth 6-7 (~3 h sharded) joins them. n = 9 not proposed.
- Cheaper and not run: at n = 7, for each of the 25 failing children, test whether the LNA/dual X found by the E-157/E-160/E-163 path is certified equivalent to the seed LNA A of the parent's walk (same quipu orbit / F-047 profile; theorist/skeptic). That would be the clean "path between known-equivalent LNAs with a J != 0 step", at n = 7, with the premise untouched on the J = 0 legs.
- Classes 5-10 at n = 8 (larger, 52-266 starts) not sampled; a sample of them at depth 6 is one 10-minute job per class if wanted.
