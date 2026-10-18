# Experimentalist notebook (rewritten round 073)

## What I now believe
- Free-stretch rule (E-114/E-146/E-200/E-206): a vertex at distance D >= 3 from every relation support gives one image class; n = 9 -> 8 6/6, n = 10 -> 9 10/10 (key+orbit source labels, 382 cells). D >= 2 fails (real), D >= 1 fails more. The rule is mostly vacuous: 32 of 42 n = 10 labels have no D >= 3 cell.
- 073: the E-190 pair 90000000/00000030 cannot test the rule (90000000 has a relation over all 10 vertices, D = 0 everywhere). Unrestricted, every deletion of 90000000 gives 80000000 and 00000030 gives 00000300 (i <= 7) or 0000002 (i >= 8): same class P^(6)_(1,1) except the three deletions inside the relation end. Agreement is 00000030's own 7-of-10 behaviour, not a transport rule.
- Exact classifyLength(10): the quipu seed + move orbits place 2967/4862 rows (61.0 %, 36 classes) in 311 s; 1895 unplaced rows, timed searches 0.9-5.7 s at depth 3-5 (main pass <= ~3 h worst case); the resolve step is unsized. E-032 recorded 1010 rows/9 h on the old engine. classify.py has no --plan but has --budget-hours and exits 2 on budget (already what the OVERNIGHT proposal needs).
- Coxeter key == class at n = 8, 9 (067); no fixed i or class+i+coverage rule predicts the image (E-194, E-200). T4 floating H-020 widths 6..8 pass; T7 right LNA rule exact n = 5..10.

## What I tried
- 073: `workshop/rounds/073/experimentalist_pair.py` (2 s), `experimentalist_plan10.py` (plan 5.5 min; --run checked at n = 7 only: exit 0, the exit-2 path is classifyLength's own and not exercised).
- 070: `rounds/070/experimentalist_free10.py`, `_free10_d2.py`, `_pool.py`. 067: `_deltable.py`, `_delsplit.py`. 065: `_delcompat.py`, n = 9 class CSV.
- Lessons: grep EXPERIMENTS.md for E-114/E-146 first; background jobs piped through tail show nothing until the end; any command over ~2 min goes to the background.

## What I would do next
1. Run the --run loop overnight (OVERNIGHT Menu 4), then redo free10 on exact labels; first check the --run path on n = 8 (seconds).
2. Is 80000000 ~ 00000300 at n = 9 joined by orbit+mirror moves (E-190 one length down)?
3. S-1 item 2 for H-020 head/tail at n = 12 -> 11 on the lone-3 class (E-146).
4. n = 11 -> 10 with the D >= 3 filter to cover more unfitted labels (size first).
