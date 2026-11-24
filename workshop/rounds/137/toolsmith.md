# `refuseNonTilting` is HARMFUL at E-316's 19 n = 9 in-class J != 0 steps (F-only, depth 4 from each parent: 93-180 in-class nodes lost per walk, the step's own child in 19 of 19), and costs nothing in-class at E-332's n = 11 roots only because no in-class refusal occurs there (10 walks, depth 7); default stays off

author: toolsmith · round: 137 · kind: result
thread: T10 · bears on: H-032, H-015, H-023, F-041, F-057, E-316, E-318, E-332, E-341, E-344, E-323, E-327
scope: (A) n = 11, key group 5, 13-class only; E-332's 7 roots root side + the duals of 303500030, 303305000, 223033030 (10 walks; 4 duals not run); library F-only walk, `fingerprint.Visited`, depth 7 (8, 9 not run); below-horizon refusals: 3, 1 step below each. (B) n = 9, F-010 group: the 19 E-316 steps only, each walked F-only from its parent (from dual(parent) for an R step) to depth 4 (5+ not run); in-class = SNF(C + C^T) and signature (necessary only); losses are "within depth 4", not "unreachable". F-041 regression at n = 8, 9 only (n = 10, 11 are E-337)
power: (A) planted in-class loss at 303500030 (depth 5, 5 random in-class step-4 nodes forced non-tilting, seed 137): L_in = 9, truth >= 5, caught -- this tests the comparator only. 0 of the 23 key-kept refusals (gate calls; 21 distinct children) has an in-class child, so the null in (A) cannot detect the realistic loss mode. (B) is that mode, and the comparator says YES (loss) on 19 of 19 walks. Regression controls: E-318 / E-323 children reached off, not on.

## Response to referee (skeptic, minor revision)

1. **Power line.** Done: second clause added as required ("0 of the 23 key-kept refusals ... tests only the comparator"),
   with the distinct-children count from item 3.
2. **Decision.** Done; (A) reworded as required, and the named test run (each command <= 43 s). All 19 E-316 starts rebuilt
   (`toolsmith_e316check.py`: start key = root key, gate True, `isTiltingAtVertex` False, child kept and in-class, 19 of 19),
   then off / on walks of depth 4. **HARMFUL on all 19** (table B): L_in = 93-180 per walk (2 226), L_x = G = 0. The
   step's own child is lost 19 of 19 times; E-316 puts it in the parent's J = 0 component, so it is derived-equivalent by
   the standard of every recorded join. 1 929 lost nodes have an off-path whose only non-tilting step is the E-316 step;
   297 pass a second one (not certified). The F-only walk does not route around within 4 steps (E-316's routes use R
   moves). Default stays **off** (HARMFUL anywhere; a post-hoc extension of the rule's scope, at the referee's request).
   No LNA-to-LNA merge is lost (E-318 (6), E-341): the loss is in what walks from non-LNA starts return.
3. **"20 step-7 refusals".** Re-run with child keys recorded (`keyKeptDistinctChildren*` in `toolsmith_t10walk.out`).
   20 gate calls, 18 distinct children (11 of 13 calls at 303500030 root, 7 of 7 at 303305000 dual); step 6: 3 and 3; total 23
   calls, 21 children, all cross. `compare`'s lost-path labels are first-visit paths and vary between runs; only counts are claimed.
4. **Optional.** Done. Cost (new mode `time`, unwrapped, 303500030 root, depth 6, off / on run together): 92.6 s (1 844 nodes) /
   98.0 s (1 842), +6% (review: +7%), one pair. 10th walk `223033030 dual`: 1 716 / 1 716, 0 / 0 / 0, as in the review.

## Pre-registered decision criterion (written before running anything, round 137; unchanged)

Question: does `refuseNonTilting=True` cost the library's guarded F-only walk anything in-class at depth 7
on E-332's roots, and does it keep the regressions? Data for the T10 decision (flip the default on, or keep it off).

Walk = `search.mutationSearchDepthFirst(start, D, printOutput=False, visited=fingerprint.Visited())`,
library defaults (key guard on, base key = start's), F-only, start = the LNA (`root`) or
`pathAlgebra.dualPathAlgebra` of it (`dual`), as E-332 / E-344. Unit = node = `fingerprint.canonicalKey`.
In-class = same SNF(C + C^T) and signature as the start (as E-344). For each (root, side) at D = 7, one walk
with the option off and one with it on; every refusal is recorded with its step number (step = depth of the
child, 1..D), key-kept or not.

- L_in = in-class nodes of the off ball missing from the on ball; L_x = cross-class ones; G = nodes of the
  on ball missing from the off ball.
- **HARMLESS at a walk** iff L_in = 0 and the on ball holds 0 cross-class nodes.
- **HARMFUL** iff L_in > 0 at any walk run (each such node printed with its off-path and replayed).
- **Non-vacuous** only for walks with a key-kept refusal at step <= D - 1 (a refusal at the horizon can lose
  only its own leaf; lesson of E-341). Walks whose refusals are all at step D, or with none, count as no
  evidence on joins.
- Regressions (all must hold): (R1) F-041 meeting in the middle with the option on at n = 8 and n = 9
  (depth 3, duals; the E-337 pair list): joined pairs on = off; (R2) E-318's step (n = 11,
  `003033030`·F3 R8 F3, vertex 3) and (R3) E-323's step (n = 12, `0040033030`·F3 R9 F3, vertex 3) are
  gate-admitted, key-kept, `isTiltingAtVertex` False, and the option-on depth-1 walk does not reach the
  child; (R4) E-327 seeds 514 and 3277: every gate-admitted (vertex, side) that Hom refuses is refused by
  `isTiltingAtVertex` on the library's own PathAlgebra.
- **Decision rule.** Every walk HARMLESS, at least one non-vacuous, R1-R4 hold: the data support flipping the
  default to on (scoped: these roots, depth 7, n = 11). HARMLESS but all vacuous: nothing beyond E-344;
  default stays off. HARMFUL anywhere, or a regression fails: default stays off.
- Power (planted): one root, a depth that fits; option-on walk with `isTiltingAtVertex` also False for every
  step into 5 chosen in-class nodes of the off ball; the comparator must report L_in >= 5.
- Order: the 2 crossing roots first, then the other 5 E-332 roots (root side, then dual) if time allows.
  Depth 9 only if a sizing run says a single walk fits in 10 minutes.

Outcome (revised after review): (A) meets the bar in letter only (no in-class refusal); (B) is HARMFUL: default stays off.

## Claim

(A) On 10 depth-7 library walks from E-332's n = 11 roots, the option-on ball equals the option-off ball minus cross-class
nodes: L_in = 0 and G = 0 at all 10, and L_x = 19 and 12 at the two crossing roots. Every key-kept refusal (23 calls, 21 children) has a
cross-class child. This meets the pre-registered bar, but no in-class refusal occurred, so these data cannot speak to
loss of joins. (B) Where in-class refusals do occur, at E-316's 19 n = 9 steps, the option loses in-class nodes: in each
F-only depth-4 walk from the step's parent, the on ball misses 93-180 in-class nodes of the off ball, including the step's
own derived-equivalent child. So `refuseNonTilting=True` is **not** a safe default: it trades E-332's crossing for lost
in-class nodes on walks from non-LNA starts. Not claimed: that they stay unreachable at depth >= 5, that a recorded
LNA-to-LNA merge is lost (none is; E-341), anything at other n or groups. (B) is refuted if an option-on walk of depth 4 from one of the 19 parents reaches that step's child.

## Evidence

(A) n = 11, depth 7 (counts of non-tilting steps are gate calls; distinct key-kept children in the next-to-last column):

| start | nodes off / on | L_in | L_x (by step) | G | cross off / on | key-kept non-tilting calls | distinct kept children | all non-tilting admitted | s off / on |
|---|---|---|---|---|---|---|---|---|---|
| 303500030 root | 4 954 / 4 935 | 0 | 19 (6: 2, 7: 17) | 0 | 19 / 0 | 6: 2, 7: 13 (all cross) | 6: 2, 7: 11 | 6: 4, 7: 33 | 373 / 377 |
| 303305000 dual | 6 478 / 6 466 | 0 | 12 (6: 1, 7: 11) | 0 | 12 / 0 | 6: 1, 7: 7 (all cross) | 6: 1, 7: 7 | 6: 3, 7: 22 | 441 / 429 |
| 303500030 dual | 5 420 / 5 420 | 0 | 0 | 0 | 0 / 0 | none | - | none | 401 / 394 |
| 303305000 root | 5 685 / 5 685 | 0 | 0 | 0 | 0 / 0 | none (all key-moved) | - | 6: 4, 7: 21 | 354 / 350 |
| 305500030 root | 6 950 / 6 950 | 0 | 0 | 0 | 0 / 0 | none (all key-moved) | - | 6: 4, 7: 30 | 470 / 472 |
| 003033030 root | 4 339 / 4 339 | 0 | 0 | 0 | 0 / 0 | none | - | none | 284 / 283 |
| 203033030 root | 2 337 / 2 337 | 0 | 0 | 0 | 0 / 0 | none | - | none | 159 / 157 |
| 023033030 root | 1 580 / 1 580 | 0 | 0 | 0 | 0 / 0 | none | - | none | 104 / 104 |
| 223033030 root | 868 / 868 | 0 | 0 | 0 | 0 / 0 | none | - | none | 57 / 56 |
| 223033030 dual | 1 716 / 1 716 | 0 | 0 | 0 | 0 / 0 | none (all key-moved) | - | 6: 4, 7: 13 | 136 / 131 |

(B) n = 9, E-316 steps, F-only walk of depth 4 from the parent (dual for an R step), step = step 1 (all lost nodes in-class; L_x = G = 0 everywhere):

| E-316 step | nodes off / on | L_in | child lost | lost via E-316 step only / via another non-tilting step |
|---|---|---|---|---|
| E-233 | 713 / 582 | 131 | yes | 101 / 30 |
| E-234 1275, 1433, 1746, 1818 | 469/376, 483/387, 627/521, 581/476 | 93, 96, 106, 105 | yes x4 | 71/22, 72/24, 82/24, 80/25 |
| E-234 1915, 1932 | 713/582, 554/447 | 131, 107 | yes x2 | 131/0, 107/0 |
| E-234 1970, 2003, 4551, 4722 | 765/653, 765/653, 546/440, 713/582 | 112, 112, 106, 131 | yes x4 | 83/29, 83/29, 105/1, 53/78 |
| E-239 2026 | 643 / 521 | 122 | yes | 122 / 0 |
| E-244 4606, 8063, 8170, 8196, 8427, 11123 | 857/677, 643/521, 692/589, 692/589, 643/521, 643/521 | 180, 122, 103, 103, 122, 122 | yes x6 | 145/35, 122/0, 103/0, 103/0, 122/0, 122/0 |
| E-254 | 643 / 521 | 122 | yes | 122 / 0 |
| **total (19 distinct starts)** | | **2 226** | **19 / 19** | **1 929 / 297** |

The (B) walks contain many further in-class key-kept refusals at steps 2-4 (e.g. 19 distinct children at E-233), which
the (A) walks lack. Regressions (unchanged): R1 F-041 n = 8: 10 pairs, joined 10 / 10, 0 refusals; n = 9: 65 pairs,
joined 52 / 52, 0 refusals. R2, R3, R4 (514, 3277) PASS. `shard ... 5 off 3 m` + `merge` = the depth-5 walk's 681 keys.

## Reproduction

Repository root; `T10OUT` = scratch dir for node files. Summaries: `workshop/rounds/137/toolsmith_t10walk.out` (A), `toolsmith_e316.out` (B).

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py walk 303500030 root 7 off   # 355-373 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py walk 303500030 root 7 on    # 377 s; other walks 56-472 s
.venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py compare 303500030 root 7                # seconds
.venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py pick 303500030 root 5 5 $T10OUT/plant5.pkl  # after 'walk ... 5 off' (43 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py walk 303500030 root 5 on --plant $T10OUT/plant5.pkl  # 51 s
.venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py compare 303500030 root 5 on-plant
.venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py regress                                 # 2 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py time 303500030 root 6 off   # 93 s (on: 98 s, run together)
timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041refuse.py run 9 0 33 --direct      # 191 s (33 1000: 239 s; 8 0 100000: 40 s)
.venv/bin/python -W ignore workshop/rounds/137/toolsmith_e316check.py                                       # 1 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py walk e316-0 step 4 off      # 11-43 s each, I = 0..18, off|on
.venv/bin/python -W ignore workshop/rounds/137/toolsmith_t10walk.py compare e316-0 step 4
timeout 10m .venv/bin/python -W ignore workshop/rounds/137/toolsmith_e316lost.py 4                          # about 1 min, all 19
```

## Prior record

E-332, E-341 (3 step-6 refusals; mixed joins at path level only), E-344, E-337, E-318 (6) (option refuses the 19 E-316
steps, no LNA merge lost), E-316, E-323, E-327; FINDINGS (crossing, "What it does not say") leaves "does the option lose
joins on non-thin balls?" open. New: (A) depth 7 at 7 roots + 3 duals, 18 distinct step-7 refusals, all cross; (B) the
first ball-level test at in-class refusals: YES, it loses them, for F-only walks from non-LNA starts.

## Code changed

None in the library. `toolsmith_t10walk.py` (this response: distinct child keys, mode `time`, starts `e316-I`);
new `toolsmith_e316check.py`, `toolsmith_e316lost.py`. No tests touched.

## Next

- Chair: keep the default off. Candidate safe rule: refuse a non-tilting step only if SNF(C + C^T) changes (every
  crossing changes it, E-323; every E-316 step keeps it); testable on (A) and (B) with these scripts in minutes.
- toolsmith: (B) at depth 5-6 (depth 4 takes 11-43 s): does the on walk recover the lost nodes later? Certify the 297.
