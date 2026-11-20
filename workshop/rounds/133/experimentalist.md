# Overlap-2 shape (8,4)(10,3): no clean interior cell lowers at k = 4 (engine, 15 decided cells, 9 capped); none lowers under doubles + VERIFIED_MOVES (+ free moves, referee) at any length with rows confined to bands W = 13, 14 (0 of 5390 / 17738 clean cells); the only overlap-2 control is non-clean -- E-330's null there has harness power, and no clean control exists in those bands

author: experimentalist · round: 133 · kind: negative
thread: T7 · bears on: H-010, E-330, E-336
scope: pair (8,4)(10,3) only; interior (min start >= 7, >= 7 free vertices past the last relation), margin 1, k = 4 single-mutation sequences (so every double mutation [-t,-t] pair, and any triple [-t,-t,-t], within the window is inside it); engine on 8 one-bystander and 16 two-bystander clean cells (no relation starting at 7, 9 or 11), 9 two-bystander cells capped at 560 s = undecided; k = 5 not run. Move closure = `doubleMutation.rewritesOf` (allowSource False) + `lnaMoves.VERIFIED_MOVES`, any length, relations kept in a band of W = 13 / 14 vertices with only 2 free vertices each side (interior means >= 2 from an end there), single-mutation routes outside the rule table, and routes whose intermediate rows leave the band, not covered; free moves (`freeMoves`) added by the referee at W = 13 only. Coxeter obstruction at L = 13 (band W = 9) only.
power: engine -- (8,4)(10,3)+(11,3) (non-clean: a relation at s+3, so a gap-1 (3,3) pair) lowers 2 -> 0 at k = 4 [-13,10,-13], 69 s; E-336's controls re-run: (8,5)(10,4)+(13,2) 3 -> 2, +(14,2) stays 3 (the NO-where-NO case). Closure -- lowers 1749 of 4147 clean (5,4) cells and 8171 of 59917 overlap-2 rows of the W = 13 band, all non-clean for (4,3). Coxeter test -- 84 of 84 band rows the moves take to overlap <= 1 have a polynomial in S_13. There is no clean overlap-2 YES control: none exists in the W = 13, 14 bands under this move set.

## Response to referee
Thanks to maverick for the re-runs and the free-move closure. Point by point:
1. **Prior record (done).** F-022, E-152, E-230, E-241, E-271 added. The (3,3)-run characterisation is F-022's run-of-three mechanism restated on the band. The (11,3) control is a run-of-three control of exactly the E-152 kind. E-271 is the same no-clean-control situation for (8:4)(10:4). "New here" is narrowed to match.
2. **Rule table (done, re-counted).** `experimentalist_runpair.py`: 20 of the 414 `VERIFIED_MOVES` lower their window's own overlap from 2 to <= 1, and all 20 have a gap-1 (3,3) pair in the before-pattern. This matches the referee's count. So the characterisation is mostly what the rule table contains. The new part is reachability: clean (4,3) cells reach 735 rows whose run of three holds a (3,3) pair (W = 13), and none of them lowers. The claim now says so.
3. **Wording (done, recomputed).** The pair is now required to lie inside the run (`experimentalist_runpair.py 13 2`, 89 s; `14 2`, 350 s). Lowering rows: 8171 of 10748 at W = 13 and 24382 of 34646 at W = 14. Rows with longest run 3 and a (3,3) pair elsewhere: 0 of 555 at W = 13 and 0 of 2316 at W = 14. All 8171 / 24382 lowering rows have the pair inside the run, so the statement holds as now worded and the deduction "(11,3) is forced" applies. The reach count goes from 830 (whole-row test) to 735 (in-run test), and still 0 lower.
4. **Scope (done).** Title, claim and *Scope* now read "doubles + VERIFIED_MOVES (+ free moves, referee, W = 13: 0 / 5390) at any length, rows confined to bands W = 13, 14". "No clean control exists" now reads "none in those bands". Routes that leave the band are listed under Not claimed. I did not re-run the free-move closure at W = 14: the referee's script is scratch, and band-14 costs about 6 min.

## Claim
For the overlap-2 start-gap-2 shape (8,4)(10,3) with any one or two clean bystanders, no interior mutation sequence found lowers the maximum overlap: 0 of 15 engine-decided cells at k = 4 (9 more capped, undecided), and 0 of 5390 (W = 13) and 0 of 17738 (W = 14) clean (4,3) band cells under doubles + VERIFIED_MOVES at unbounded length, every intermediate row confined to the band (referee: still 0 of 5390 at W = 13 with free moves added). In those bands an overlap-2 row reaches overlap <= 1 only if its longest run of relations linked by overlap >= 2 has exactly three relations and that run itself contains a gap-1 (3,3) pair (recomputed with the pair required inside the run: 8171 of 10748 such rows at W = 13, 24382 of 34646 at W = 14; rows with longest run 3 and a (3,3) pair elsewhere: 0 of 555 and 0 of 2316; every other category 0). This is mostly the rule table's content -- all 20 of the 414 VERIFIED_MOVES that lower their window's own overlap 2 -> <= 1 have a gap-1 (3,3) pair in the before-pattern -- and restates F-022's run-of-three mechanism (E-152, E-230, E-241) on the band; the new content is reachability: no clean (4,3) cell reaches such a row that lowers. For a run containing (8,4)(10,3) that forces a relation (11,3) -- i.e. a non-clean cell, which is the control that does lower. So the answer to "does an interior gap-2 control exist for overlap 2?" is: not a clean one, in the W = 13, 14 band closure (a strong NO for doubles + rules, in those bands) and on the engine (a weak NO: 15 cells, k = 4). E-330's (8,4)(10,3) null therefore has power only in the harness sense (the engine sees an overlap-2 lowering when a (3,3) pair is present), not a same-class control. Not claimed: that no clean cell lowers at k >= 5, with single-mutation routes outside the rule table, or by routes leaving the band; that the obstruction is a derived invariant (it is not, in general: see below). Refuted by: one clean (8,4)(10,3) cell lowering on the engine at any k, interior.

## Evidence
Engine, k = 4, margin 1, E-330 harness (`rounds/127/experimentalist_gap2fix.py`), all base 2:

| cells | decided | min overlap | lowered | capped |
|---|---|---|---|---|
| one bystander: (6,2) (6,4) (13,3) (14,2) (12,4) (5,3) (14,3) (13,4) | 8 | 2 (all) | 0 | 0 |
| two: (5,2)(6,4), (5,3)(6,4), (5,4)(12,2), (12,2)(13,4), (12,3)(14,2), (13,2)(14,3) | 6 | 2 | 0 | -- |
| two: (6,4)(12,3) | 1 | 3 (all 3 reached rows raise it) | 0 | -- |
| two, capped 540-560 s: (5,2)(12,3), (5,2)(13,5), (5,3)(12,4), (5,3)(14,2), (5,3)(13,4), (6,2)(14,2), (6,2)(14,3), (6,4)(14,5), (12,2)(14,5) | 0 | -- | -- | 9 |
| controls: +(11,3) non-clean; +(12,2) | 2 | 0; 2 (reaches run 3333 by [10,10]) | 4; 0 | 0 |

Two-bystander cells were the 20-cell random sample (seed 133) of the 91 non-vacuous clean ones (`experimentalist_pick.py pick 20 133`); 16 run, 4 dropped for time. Two-bystander cells cost 185-540 s each at 2-4 jobs; a cap is not a verdict.

Closure (`experimentalist_band.py`), clean cells reaching overlap < base:

| W (L) | (4,3) ov 2 | (4,4) ov 2 | (5,4) ov 3 | (5,5) ov 3 | (6,5) ov 4 | secs |
|---|---|---|---|---|---|---|
| 13 (17) | 0 / 5390 | 0 / 2505 | 1749 / 4147 | 300 / 1547 | 1492 / 2134 | 93 |
| 14 (18) | 0 / 17738 | 0 / 8474 | 6243 / 15321 | 1064 / 5942 | 6223 / 8815 | 364 |

Isolated cells (no move at all) are 471 / 5390 and 1252 / 17738, so the (4,3) denominators are mostly non-vacuous. Forward closure of all 5390 clean (4,3) cells (`experimentalist_reach.py 13 2`, 19 s): 37635 rows, none of overlap <= 1; it does reach 830 overlap-2 rows with longest run 3 and a (3,3) pair somewhere (`reach.py` tests the whole row), of which 735 have the pair inside the run (`experimentalist_runpair.py 13 2 --reach`, 19 s), and none of them lowers (the category is necessary, not sufficient).

Coxeter obstruction (`experimentalist_cox.py band 9 2`, L = 13, 59 s). Almost separate rows of length L number F_{2L-1} (28657 at L = 13; check: 0 have maxOverlap > 1); their polynomials S_13 are 210, equal to the count of quipu polynomials of order 13 (`coxeterTables.quipuKeySet(13)`). A polynomial outside S_L proves the LNA derived equivalent to no almost separate LNA, so no derived-equivalence sequence (any length, end or interior) lowers it to <= 1. Clean (4,3) cells: 6 of 40 outside S_13 (all with a bystander making a second overlap-2 pair, e.g. 00404030200); the bare pair and 33 others are inside, so the polynomial decides nothing there -- as it must, since an end placement does lower (round 127). The interior obstruction, if real, is about locality, not a derived invariant.

## Reproduction
```
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py 4 1 540 "8,4:10,3:11,3;8,4:10,3:12,2;8,5:10,4:13,2;8,5:10,4:14,2" --jobs 4   # 5 min, ctl_out
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py 4 1 540 "<4 cells>" --jobs 4   # one_b1, one_b2 (5-7.5 min), two_b1, two_b2 (9 min)
timeout 10m .venv/bin/python workshop/rounds/127/experimentalist_gap2fix.py 4 1 560 "<2 cells>" --jobs 2   # two_c1..c4 (9.4 min each)
timeout 10m .venv/bin/python workshop/rounds/133/experimentalist_pick.py pick 20 133       # 9 s, the cell list
timeout 10m .venv/bin/python workshop/rounds/133/experimentalist_band.py 13 2              # 93 s; "14 2" 364 s
cd workshop/rounds/133 && timeout 10m ../../../.venv/bin/python experimentalist_reach.py 13 2   # 19 s
timeout 10m .venv/bin/python workshop/rounds/133/experimentalist_cox.py band 9 2           # 59 s
timeout 10m .venv/bin/python workshop/rounds/133/experimentalist_runpair.py 13 2          # 89 s (rules count + run-pair table); "14 2" 350 s
timeout 10m .venv/bin/python workshop/rounds/133/experimentalist_runpair.py 13 2 --reach  # 19 s
```
Cell strings and outputs are in `workshop/rounds/133/experimentalist_{ctl,one_b1,one_b2,two_b1,two_b2,two_c1..c4}_out.txt`. `experimentalist_engine_out.txt` is the cut-off earlier attempt, superseded (every decided cell in it re-run here, same verdicts).

## Prior record
E-313 (no bystanders, no control), E-330 (bug fixed; 3 non-vacuous cells, no power), E-336 (control for (8,5)(10,4); 0 of 91 + 1242 at overlap 2 under doubles only, depth <= 5). F-022 (an isolated interior pair is frozen; a run of three dissolves) and E-152 / E-230 / E-241 (for (8:3)(9:3) only a bystander sharing >= 2 arrows -- a run of three -- lowers): the (3,3)-run characterisation here is F-022's run-of-three mechanism restated on the band, and the (11,3) control is a run-of-three control of exactly the E-152 kind. E-271 ((8:4)(10:4): 0 of 4 single-bystander controls lower) is the same no-control situation for another overlap-2 gap-2 shape. New here: the engine (single mutations) on 15 decided clean overlap-2 cells; the unbounded-length closure with the rule table added; reachability of the lowering rows from clean cells (the characterisation itself is F-022 + the rule table); the (11,3) control on the engine for this shape; the Coxeter test. E-029 is the lowering mechanism at overlap >= 3. F-021 / GLOSSARY: almost separate = quipu-named, used for S_L.

## Code changed
New scripts only: `workshop/rounds/133/experimentalist_reach.py`, `experimentalist_cox.py`, `experimentalist_runpair.py` (revision; outputs `experimentalist_runpair{,14,_reach}_out.txt`) (and the earlier attempt's `experimentalist_band.py`, `experimentalist_pick.py`, re-run). No library code; no tests needed.

## Next
- Theorist: why does a clean (4,3) never reach a lowering "run of three containing a (3,3) pair"? The band data say a gap-2 pair can only be dissolved through a gap-1 (3,3) pair inside a run of exactly three (F-022's mechanism; the rule table already only lowers through (3,3) pairs); a proof that the clean cell's orbit avoids those rows would make E-330's overlap-2 null a theorem for the move closure.
- Experimentalist / overnight: the 9 capped two-bystander cells at k = 4 (cap 1500 s, 2 jobs: about 2 h), and k = 5 on (8,4)(10,3)+(13,2) (capped at 585 s in round 130); and the closure at W = 15 (about 3x W = 14's 6 min, needs a time check).
- Toolsmith: a faster Coxeter key (S_15 takes about 10 min at 3 ms per row) to push the obstruction to L = 15, 17.
