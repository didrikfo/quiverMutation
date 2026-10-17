# No cyclic parent occurs in the key-guarded n = 7 balls (parent depth <= 8, complete) nor in the n = 8 balls searched (class 0 through depth 8, class 2 through depth 7 and part of 8, budget-stopped); those n = 7 balls hold 70 failing edges versus 25 edges at a capped level; an n = 7 inequivalent control is key-forced, the live control is n = 9 F-010 at depth 7 + 7


author: toolsmith · round: 070 · kind: result
thread: T10/T5 (agenda item 1) · bears on: H-015, E-192, E-193, E-198, E-199, E-157
scope: n = 8 class 2 is budget-stopped (partial level 8); n = 7 classes 1, 2 complete through parent depth 8 (key-guarded BFS, canonicalKey dedupe, no expansion cap); n = 8 class 0 and class 2 through depth 8 (class 2 stopped by a 540 s budget inside level 8, 14 941 expansions); controls n = 9 F-010 pair at J0 depth 7 + 7. Everything "joined" below is "under the J = 0 premise".

## Response to referee

1. Reconciled. Two counters: `edges` (all legal children) = 277 894 = 122 298 (c1) + 155 596 (c2); `guard` (key-kept, Coxeter key = base) = 277 669 = 122 207 + 155 462; the difference is `noguard` (key-moving) = 91 + 134 = 225. The text's "key-kept" figures are the `guard` counter; the referee's 277 894 is `edges`. Both from shard sums minus 3 x the depth-8 prefix, in `toolsmith_cyccount_n7_sums.txt`.
2. Done. `toolsmith_cyccount.py` now tallies `x-{guard,noguard}-{tilt,FAIL}`. Re-ran all 8 depth-9 shards and both depth-8 prefixes (outputs overwritten in `toolsmith_cyccount_n7c{1,2}*.txt`; counts identical to before): x-guard-tilt 277 599, x-guard-FAIL 70, x-noguard-tilt 0, x-noguard-FAIL 225. So claim (3) "no gate-admitted tiltingPlus-true step moves the key" is checkable in the files for these balls (0 of 277 599 + 0 of 225 noguard tilt-true; all 225 key-moving edges are J != 0).
3. Title and wording changed: n = 8 class 2 is budget-stopped, so "none at n = 8" is now stated per ball; "not 25" became "70 edges versus 25 edges at a capped level" (E-192 never claimed completeness).
4. Optional key-preserving cycle injection (exercising `cyclic-guard-*`, D(w)): deferred to next round, not done this sitting.

## Claim

(1) Cyclic parents. In the E-192 balls no child, key-kept or not, legal or not, has a loop or oriented cycle, before or after `reducePathAlgebra`: 0 cyclic children in 277 669 key-kept (`guard`) n = 7 edges (277 894 `edges` in all, 225 key-moving) (parent depth <= 8, 67 943 expansions) and 0 in 71 605 n = 8 key-kept edges (c0 9 402 edges, 2 257 expansions; c2 62 203 edges, 14 941 expansions). So the number of cyclic parents whose cycle meets D(i) is 0 here (D(i) never even needs to be computed: there is no cycle), and E-199's reduction of the E-193 condition is not used by any step of these balls. Every one of the 70 failing children is acyclic. This does not cover other classes, deeper parents, n = 9, steps that the gate refuses, or J != 0 children expanded further (the walk does not expand failing children specially; they are nodes like any other).

(2) The failing-step count. The "25 of 174 804" of E-192 (E-151: 16 + 9) is a capped count (the 065 / 047 `deepreplay` runs were checkpointed inside level 9). Expanding every parent of depth <= 8 gives 34 (c1) + 36 (c2) = 70 key-kept `tiltingPlus`-false edges in 277 669 key-kept edges (of 277 894 legal edges) (2.5e-4 per edge), first at parent depth 7 (4 + 6, matching E-192's 10) and 8 (30 + 30; E-192 had 15 from a partial level). These are edges, not distinct children (merge edges not removed). So "all 25 failing children joined under the J = 0 premise" (E-157/E-160/E-163) covers at most 25 of >= 45-70 depending on distinctness; the rest are unjoined, not refuted. A refutation: a canonicalKey count of distinct failing children below 25.

(3) A NO-capable control. At n = 7 every certified-inequivalent pair has different Coxeter keys (E-172: at n <= 9 the only equal-key certified pair is the n = 9 F-010 pair). Any J = 0 tilting step keeps the key (0 gate-admitted `tiltingPlus`-true steps move the key: `('noguard','tilt')` absent in all ten 065 outputs, and in mine `x-noguard-tilt` = 0 of 225, `x-guard-FAIL` = 70, files `toolsmith_cyccount_n7_sums.txt` and the shard outputs), so a join test from a failing n = 7 child to an LNA of another class says NO by construction, whatever the search does. It is a check of the implementation, not of the premise, and I did not run it. The control that can say NO for a reason is n = 9 F-010 (3060000 vs 3304000, certified by Ladkani Cor 3.15 as E-175): J0 depth 7 + 7, reach sets 10 526 and 5 625, 0 iso meetings, 333 s on one core (E-175 stopped at 6 + 6); gate + key mode (ALL) gives identical reach sets (10 526, 5 625), so no J != 0 step enters the depth-7 balls either and this control still cannot test the premise. Sensitivity from E-175: equivalent controls join at total 5 and 9 (6 + 6). Depth 8 + 8 would be about 16 min per mode (x2.8 per level: 43, 121, 333 s at depth 5, 6, 7): `OVERNIGHT.md` proposal below.

## Evidence

Table (key-guarded walk, tilting computed on every edge; fails = guard-kept and `tiltingPlus` false):

| n, class | parent depth | expansions | key-kept (`guard`) edges | cyclic children | fails | fails by parent depth |
|---|---|---|---|---|---|---|
| 7, 1 | <= 8 (complete) | 29 515 | 122 207 | 0 | 34 | 4 at 7, 30 at 8 |
| 7, 2 | <= 8 (complete) | 38 428 | 155 462 | 0 | 36 | 6 at 7, 30 at 8 |
| 8, 0 | <= 8 | 2 257 | 9 402 | 0 | 0 | |
| 8, 2 | <= 7, partial 8 (540 s budget) | 14 941 | 62 203 | 0 | 0 | |

n = 7 numbers are sums over four shards of the last level minus the repeated prefix (script header). Prefix levels reproduce the 047 `deepreplay` counts exactly (expanded 12, 56, 154, 371, 857, 1988, 4855, 12033 for c1). Cyclic counts also cover: raw (pre-reduction) children, gate-admitted children refused as illegal (0 illegal at all), non-key children (`cyclic-noguard-*` = 0).

Positive control for the detector: `cyclicVertices` is `strongly_connected_components` plus self-loops; E-198's injected-cycle control exercises `hasOrientedCycle`, not this function. I checked `cyclicVertices` by hand on a 2-cycle, a loop and an acyclic quiver, but did not inject a cycle into a ball (a gap; the totals of 0 rest on one detector). Consistency: deepreplay itself tallies `parent-cyclic` and the 065 outputs have no such line, which independently says that no cyclic node was expanded in those balls.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/070/toolsmith_cyccount.py 7 1 9 570 0/4   # shards 0/4..3/4, 411 s each, 4 in parallel; also class 2; 7 c1 8 540 without shard (273 s); 8 0 8 540; 8 2 8 540
timeout 10m .venv/bin/python -W ignore workshop/rounds/057/experimentalist_powerjoin.py 9 3060000 3304000 7 J0     # 333 s
```
Outputs: `workshop/rounds/070/toolsmith_cyccount_n*.txt` (small), `toolsmith_f010_d7.txt`.

## Prior record

E-168 (n = 7 paths acyclic), E-198 (printed paths acyclic), E-193 (cyclic child exists on random walks, seed 12), E-199 (cycle must meet D(i)), E-192 (sizing, 25), E-175 / E-178 (inequivalent controls at n = 9, 10), E-172 (no equal-key certified pair at n <= 7). New: the ball-level cyclic count (0), the complete level-8 failing count (70 versus a capped 25), the F-010 depth 7 + 7 non-join. Not a rediscovery of anything recorded; the 70 contradicts nothing (E-151 and E-192 are stated as capped samples) but changes the denominator of "25 of 25".

## Code changed

`workshop/rounds/070/toolsmith_cyccount.py` only (new). No library change, no tests needed.

## Next

Toolsmith / experimentalist: (a) distinct canonicalKey count of the 70 failing children and the join of the ~45 not in E-157 (`toolsmith_tiltpath.py` machinery, J = 0); (b) `OVERNIGHT.md`: F-010 J0 depth 8 + 8, about 16 min per side-mode on one core, shard by frontier; (c) skeptic: inject a 2-cycle into a child to validate `cyclicVertices`; the cycle-meeting-D(i) case is empty on these balls, so the theorist's finite-dimensional cyclic example stays a theory item, not a repo case.
