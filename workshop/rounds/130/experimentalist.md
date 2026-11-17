# An interior start-gap-2 control exists for (8,5)(10,4): one bystander, (13,2) or (12,3), ending at or past the pair's far end, makes it lower 3 -> 2 under E-330's harness (margin 1, k = 4), while (14,2) does not; for the overlap-2 shape (8,4)(10,3) no control was found (5 engine cells, k <= 5; 0 of 91 + 1242 non-vacuous screened cells), so only one of E-330's three non-vacuous cells acquires power

author: experimentalist · round: 130 · kind: result (revised after review)
thread: T7 · bears on: H-010, E-330, E-313, E-288, E-271, E-029
scope: pairs (8,a)(10,b), a <= 7, b <= 7 (20 shapes), interior (min start >= 7, ends >= 7 vertices away), margin 1. Engine (`localMutationSequences`, unchanged E-330 harness): 8 pair+bystander cells at k = 4, 2 at k = 5 (one capped at 585 s, undecided), run of four at k = 3, 4, 5. Screen: closure under `doubleMutation.rewritesOf` only (2-mutation moves; misses sequences made of single mutations), on a quiver one vertex longer than the engine's (screen L = hi + 9, engine L = hi + 8), depth <= 4 (one bystander, s 2..20, m 2..8) and <= 5 (two bystanders, overlap-2 pairs, s 2..17, m 2..7). Counts are of non-vacuous cells (closure reaches at least one rewrite): 1150 of 1600 one-bystander and 1212 of 2454 two-bystander clean cells are vacuous, and (8,3)(10,b), b >= 3, is vacuous throughout. Not covered: k >= 6 by the engine; three or more bystanders; any proof that an overlap-2 control cannot exist.
power: YES side: the screen's 3 depth-2 hits on E-330 shapes were re-run on the engine and all 3 lower (3 of 3). NO side (overlap 2 has no control): at depth 3 the same screen lowers 30 of 281 non-vacuous clean cells of pair overlap >= 3 (46 of 281 at depth 4), against 0 of 91 at overlap 2. This is power for the double-mutation family only; single-mutation sequences at overlap 2 are tested only by the 4 engine cells at k = 4 and 1 decided cell at k = 5.

## Response to referee

1. *Narrow the first Claim sentence; record (14,2).* Done. The Claim is now existential. I re-ran (8,5)(10,4)+(14,2) at k = 4 on the engine: base 3, 11 LNAs, min 3, lowered 0, in 3 min 30 s (the referee's 3 min 54 s run agrees). It is in the engine table. The screen agrees: of 86 non-vacuous clean bystanders of (8,5)(10,4), 2 lower at depth 2 and 3 at depth 3.
2. *Vacuous / non-vacuous counts in table and title.* Done, with `experimentalist_vacuity.py` (imports `experimentalist_screen.screen` unchanged). My counts match the referee's table, with one exception. On the b >= a, overlap >= 3 row I get 504 clean and 22 non-vacuous, where the referee has 510 and 28 (the 482 vacuous and 15 lowering agree, so 22 is the difference). I use 15 of 22. Title: 0 of 91 + 1242. The fully vacuous (8,3)(10,b), b >= 3, is stated, and "none of the 11 shapes" is gone.
3. *Power NO side with non-vacuous denominators.* Done (30 of 281 against 0 of 91). I also went one depth further: at depth 4 it is 46 of 281 against still 0 of 91.
4. *Depth-3 run and per-overlap counts in `experimentalist_screen_out.txt`.* Done. Depths 2, 3, 4 (one bystander), depth 5 (two bystanders) and the isolated-pair counts are appended, about 1 s each (5 s for two bystanders).
5. *Depth of "exactly one rewrite"; L mismatch.* At depth 1 the isolated pair has exactly one rewrite (R at the first relation, [8,8], which raises the overlap). That holds only for the b = a - 1 shapes. They reach 4 rows at depth 2 and 6 at depth 3. The b >= a shapes reach 0 at every depth. The screen uses L = hi + 9 and the engine uses L = hi + 8, so they do not run on the same quiver. This is stated in Scope. The screen hits were confirmed on the engine's own L.
- Not done this sitting: the k = 5 batch with the capped cell (9 min 46 s, at the cap; an overnight run is needed to decide it).
- The referee's lead, that the overlap-2 companion (12,2) fires but yields the gap-1 run (8..11,3) rather than overlap 1, is added to Next.

## Claim

Placing E-330's shape (8,5)(10,4) in the interior with one extra relation can give cells that lower the maximum overlap from 3 to 2 under exactly E-330's harness (interior, margin 1, k = 4). The ones found end at or just past the pair's far end. (13,2) shares one arrow and lowers by [10, 9, 10, 9]. (12,3) shares two arrows and lowers by [10, 10, 10, 10]. Not every such bystander works: (14,2) shares no arrow and does not start next to a pair start, but it does not lower (k = 4, min 3, lowered 0). The same pattern holds for (8,5)(10,5) + (14,2). The base overlap 3 is the pair's own in each cell, and the lowered LNA has lost the pair's overlap ((8,5)(10,4)(13,2) -> (8,3)(9,4)(11,4), overlaps 2, 2). So E-330's null for (8,5)(10,4) at k = 4 now has a control on the same side of the question: the harness says YES on an interior gap-2 cell that differs from it only by one bystander. This is E-029's interior double mutation (the bystander ending one vertex after the pair is the companion R_s needs).

For the overlap-2 shape (8,4)(10,3), which is E-330's other two non-vacuous cells (k = 4, 5), I found no control. On the engine, 4 single-bystander cells at k = 4 and 1 at k = 5 do not lower, and 1 is capped at k = 5. In the double-mutation screen, 0 of 91 non-vacuous clean one-bystander cells lower (to depth 4), and so do 0 of 1242 non-vacuous clean two-bystander cells (to depth 5), for (8,4)(10,3) and (8,4)(10,4). The screen's contrast is narrower than a 20-shape split:
- On the (8,a)(10,a-1) family, which is the only one with dense coverage, it lowers 0, 0, 3, 5, 7 of 78, 83, 86, 87, 86 non-vacuous cells for pair overlap 1 to 5 (depth 3).
- On (8,a)(10,b), b >= a, it lowers 0 of 8 at overlap 2 and 15 of 22 at overlap >= 3.
- On (8,3)(10,b), b >= 3, the screen is vacuous.

I do not claim that no control exists at overlap 2. A cell that lowers (8,4)(10,3) or (8,4)(10,4) in the interior, below base 2, with clean bystanders, would refute the contrast. Single-mutation sequences or k >= 6 could hold one.

## Evidence

Engine, `experimentalist_gap2fix.py` (round 127, unchanged), interior, margin 1 (full output `experimentalist_engine_out.txt`):

| cell | k | L | base | LNAs | min | lowered | best (sequence) | s |
|---|---|---|---|---|---|---|---|---|
| (8,5)(10,4)+(13,2) | 4 | 22 | 3 | 11 | **2** | 1 | (8,3)(9,4)(11,4) [10,9,10,9] | 153 |
| (8,5)(10,4)+(12,3) | 4 | 22 | 3 | 10 | **2** | 2 | [10,10,10,10] | 139 |
| (8,5)(10,5)+(14,2) | 4 | 23 | 3 | 5 | **2** | 1 | (8,3)(9,5)(12,4) [10,9,10,9] | 219 |
| (8,5)(10,4)+(14,2) | 4 | 23 | 3 | 11 | 3 | 0 | (non-lowering, same kind) | 210 |
| (8,4)(10,3)+(12,2) | 4 | 21 | 2 | 16 | 2 | 0 | run of four (8..11,3) [10,10] | 114 |
| (8,4)(10,3)+(12,3) | 4 | 22 | 2 | 4 | 2 | 0 | | 169 |
| (8,4)(10,3)+(13,2) | 4 | 22 | 2 | 18 | 2 | 0 | | 191 |
| (6,3)+(8,4)(10,3) | 4 | 21 | 2 | 2 | 2 | 0 | | 161 |
| (8,4)(10,3)+(12,2) | 5 | 21 | 2 | 25 | 2 | 0 | | 252 |
| (8,4)(10,3)+(13,2) | 5 | | | | | | CAP 585 s, undecided | |
| run (8,3)(9,3)(10,3)(11,3) | 3/4/5 | 21 | 2 | 8/16/26 | 2 | 0 | | 8/38/166 |

(For comparison, E-330's isolated (8,5)(10,4) at k = 4: 4 LNAs, min 3.) The (12,2) cell reaches a run of four length-3 relations in 2 mutations. That run does not lower within margin 1 at k <= 5, so this route gives no control at k <= 7 inside the window.

Screen, `experimentalist_screen.py` / `experimentalist_vacuity.py` (`--clean`: base equals the pair's own overlap, and no bystander start is at a pair start +-1). Clean / vacuous / non-vacuous cells, and lowering cells at depth 2 / 3 / 4:

| pair overlap | clean | vacuous | non-vacuous | lower d2 | d3 | d4 |
|---|---|---|---|---|---|---|
| 1 ((8,3)(10,b), 6 shapes; b >= 3 all vacuous) | 428 | 350 | 78 | 0 | 0 | 0 |
| 2 ((8,4)(10,b), 5 shapes) | 394 | 303 | 91 | 0 | 0 | 0 |
| 3 (4 shapes) | 337 | 242 | 95 | 5 | 7 | 10 |
| 4 (3 shapes) | 263 | 168 | 95 | 7 | 11 | 17 |
| 5 (2 shapes) | 178 | 87 | 91 | 7 | 12 | 19 |
| 2, two bystanders (depth 5) | 2454 | 1212 | 1242 | 0 | | |

By family at depth 3, non-vacuous / lower: b = a - 1 gives 78/0, 83/0, 86/3, 87/5, 86/7 for overlap 1..5. b >= a gives 0/0, 8/0, 9/4, 8/6, 5/5. Two bystanders: (8,4)(10,3) has 1053 non-vacuous cells and (8,4)(10,4) has 189, with 0 lowering in each. On the isolated pair, a b = a - 1 shape has exactly one rewrite at depth 1 (R at the first relation, [8,8], which raises the overlap), 4 at depth 2 and 6 at depth 3. A b >= a shape has none. It agrees with E-330.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/experimentalist_gap2fix.py 4 1 540 "8,5:10,4:13,2;8,5:10,5:14,2;8,5:10,4:12,3" --jobs 3   # 3 min 40 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/experimentalist_gap2fix.py 4 1 560 "8,5:10,4:14,2" --jobs 1   # 3 min 30 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/experimentalist_gap2fix.py 4 1 560 "8,4:10,3:12,2;8,4:10,3:12,3;8,4:10,3:13,2;6,3:8,4:10,3" --jobs 4   # 3 min 12 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/experimentalist_gap2fix.py 5 1 585 "8,4:10,3:13,2;8,4:10,3:12,2" --jobs 2   # 9 min 46 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/experimentalist_gap2fix.py 5 1 500 "8,3:9,3:10,3:11,3" --jobs 1   # 2 min 47 s (k = 3, 4: 8 s, 38 s)
timeout 5m .venv/bin/python -W ignore workshop/rounds/130/experimentalist_screen.py 2 --pairs "<20 shapes, see screen_out>" --bys 2..20 --mmax 8 --clean   # 1 s
timeout 5m .venv/bin/python -W ignore workshop/rounds/130/experimentalist_screen.py 5 --pairs "8,4:10,3;8,4:10,4" --bys 2..17 --mmax 7 --two --clean   # 5 s
timeout 5m .venv/bin/python -W ignore workshop/rounds/130/experimentalist_vacuity.py D        # D = 2, 3, 4; 1 s each
timeout 5m .venv/bin/python -W ignore workshop/rounds/130/experimentalist_vacuity.py 5 --two  # 6 s
timeout 5m .venv/bin/python -W ignore workshop/rounds/130/experimentalist_vacuity.py 3 --iso  # isolated pairs
```
The 20 shapes are `8,a:10,b` for 3 <= a <= 7, 2 <= b <= 7, 8 + a < 10 + b. Output: `experimentalist_screen_out.txt` (all screen runs, including the vacuity tables).

## Prior record

- E-330 / E-313: the isolated nulls, with no control.
- E-288: stepped chains (8,10,12) with first length <= 4, no lowering. (8,5)(10,4)(12,3) was not among its patterns. Its (8,5)(10,5)(12,5) is vacuous.
- E-271: 0 of 4 informative single-bystander controls for (8,4)(10,4). This agrees with the 0 of 2 non-vacuous found here.
- E-029 / H-010 (2026-09-17): interior double mutations lower the overlap when a relation crosses r. The mechanism is not new. New here: an interior gap-2 control under the E-330 harness, and the screen's overlap-2 / overlap->=3 contrast (narrow, as above). The referee's grep and mine find no record of either.

## Code changed

None in the library. New scripts `workshop/rounds/130/experimentalist_screen.py` and `experimentalist_vacuity.py`. No tests touched.

## Next

- Theorist: why the contrast at pair overlap 2. Every screen lowering takes the overlap down by exactly one step. Referee's lead: the overlap-2 companion (8,4)(10,3)+(12,2) does fire, but in 2 mutations it gives the gap-1 run (8..11,3), still overlap 2, not overlap 1. Is entering the almost-separate class from a clean interior gap-2 cell impossible? That is what H-010 needs at overlap 2.
- Experimentalist: engine on 10-20 non-vacuous clean (8,4)(10,3) two-bystander cells at k = 4 (about 3 min each), to test whether single mutations find what double mutations miss. The capped k = 5 cell needs an overnight run. Check whether (8,5)(10,4)+(13,2) goes lower still at k = 5.
