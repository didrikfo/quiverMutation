# Overlap-2 shape (8,4)(10,3): no clean interior cell lowers at k = 4 (engine, 15 decided cells, 9 capped), none lowers under the library's derived moves at any length (bands W = 13, 14: 0 of 5390 / 17738 clean cells), and the only overlap-2 control is non-clean -- E-330's null there has harness power but no clean control (none exists in the move closure)

author: experimentalist · round: 133 · kind: negative
thread: T7 · bears on: H-010, E-330, E-336
scope: pair (8,4)(10,3) only; interior (min start >= 7, >= 7 free vertices past the last relation), margin 1, k = 4 single-mutation sequences (so every double mutation [-t,-t] pair, and any triple [-t,-t,-t], within the window is inside it); engine on 8 one-bystander and 16 two-bystander clean cells (no relation starting at 7, 9 or 11), 9 two-bystander cells capped at 560 s = undecided; k = 5 not run. Move closure = `doubleMutation.rewritesOf` (allowSource False) + `lnaMoves.VERIFIED_MOVES`, any length, relations kept in a band of W = 13 / 14 vertices with only 2 free vertices each side (interior means >= 2 from an end there), single-mutation routes outside the rule table not covered. Coxeter obstruction at L = 13 (band W = 9) only.
power: engine -- (8,4)(10,3)+(11,3) (non-clean: a relation at s+3, so a gap-1 (3,3) pair) lowers 2 -> 0 at k = 4 [-13,10,-13], 69 s; E-336's controls re-run: (8,5)(10,4)+(13,2) 3 -> 2, +(14,2) stays 3 (the NO-where-NO case). Closure -- lowers 1749 of 4147 clean (5,4) cells and 8171 of 59917 overlap-2 rows of the W = 13 band, all non-clean for (4,3). Coxeter test -- 84 of 84 band rows the moves take to overlap <= 1 have a polynomial in S_13. There is no clean overlap-2 YES control: the closure says none exists in its move set.

## Claim
For the overlap-2 start-gap-2 shape (8,4)(10,3) with any one or two clean bystanders, no interior mutation sequence found lowers the maximum overlap: 0 of 15 engine-decided cells at k = 4 (9 more capped, undecided), and 0 of 5390 (W = 13) and 0 of 17738 (W = 14) clean (4,3) band cells under the library's derived moves at unbounded length. In those bands an overlap-2 row reaches overlap <= 1 only if its longest run of relations linked by overlap >= 2 has exactly three relations and contains a gap-1 (3,3) pair (8171 of 11303 such rows at W = 13, 24382 of 36962 at W = 14; every other category 0 of 48614 / 153734). For a run containing (8,4)(10,3) that forces a relation (11,3) -- i.e. a non-clean cell, which is the control that does lower. So the answer to "does an interior gap-2 control exist for overlap 2?" is: not a clean one, in the closure (a strong NO for doubles + rules) and on the engine (a weak NO: 15 cells, k = 4). E-330's (8,4)(10,3) null therefore has power only in the harness sense (the engine sees an overlap-2 lowering when a (3,3) pair is present), not a same-class control. Not claimed: that no clean cell lowers at k >= 5 or with single-mutation routes outside the rule table; that the obstruction is a derived invariant (it is not, in general: see below). Refuted by: one clean (8,4)(10,3) cell lowering on the engine at any k, interior.

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

Isolated cells (no move at all) are 471 / 5390 and 1252 / 17738, so the (4,3) denominators are mostly non-vacuous. Forward closure of all 5390 clean (4,3) cells (`experimentalist_reach.py 13 2`, 19 s): 37635 rows, none of overlap <= 1; it does reach 830 rows of the "run of three with a (3,3) pair" kind, but only ones that never lower (the category is necessary, not sufficient).

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
```
Cell strings and outputs are in `workshop/rounds/133/experimentalist_{ctl,one_b1,one_b2,two_b1,two_b2,two_c1..c4}_out.txt`. `experimentalist_engine_out.txt` is the cut-off earlier attempt, superseded (every decided cell in it re-run here, same verdicts).

## Prior record
E-313 (no bystanders, no control), E-330 (bug fixed; 3 non-vacuous cells, no power), E-336 (control for (8,5)(10,4); 0 of 91 + 1242 at overlap 2 under doubles only, depth <= 5). New here: the engine (single mutations) on 15 decided clean overlap-2 cells; the unbounded-length closure with the rule table added; the (3,3)-run characterisation of which overlap-2 band rows lower; the non-clean (11,3) control; the Coxeter test. E-029 is the lowering mechanism at overlap >= 3. F-021 / GLOSSARY: almost separate = quipu-named, used for S_L.

## Code changed
New scripts only: `workshop/rounds/133/experimentalist_reach.py`, `experimentalist_cox.py` (and the earlier attempt's `experimentalist_band.py`, `experimentalist_pick.py`, re-run). No library code; no tests needed.

## Next
- Theorist: why does a clean (4,3) never reach a lowering "run of three with a (3,3) pair"? The band data say a gap-2 pair can only be dissolved through a gap-1 (3,3) pair inside a run of exactly three; a proof that the clean cell's orbit avoids those rows would make E-330's overlap-2 null a theorem for the move closure.
- Experimentalist / overnight: the 9 capped two-bystander cells at k = 4 (cap 1500 s, 2 jobs: about 2 h), and k = 5 on (8,4)(10,3)+(13,2) (capped at 585 s in round 130); and the closure at W = 15 (about 3x W = 14's 6 min, needs a time check).
- Toolsmith: a faster Coxeter key (S_15 takes about 10 min at 3 ms per row) to push the obstruction to L = 15, 17.
