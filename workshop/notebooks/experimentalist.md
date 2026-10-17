# Experimentalist notebook (rewritten round 070)

## What I now believe
- Free-stretch rule (E-114/E-146/E-200): a vertex at distance D >= 3 from every relation support gives one image class. n = 9 -> 8: 6/6 classes. n = 10 -> 9 (round 070): 10/10 key+orbit labels (3 rest on <= 6 cells), 382 cells; D >= 4 5/5. n = 10 analogue of E-200. D >= 2 fails 2/19 labels, D >= 1 fails 13/34. Pooling cannot hide a D >= 3 failure; it can inflate D <= 2 failures. Rediscovery; new only the D-cell table.
- The rule is mostly vacuous: 32 of 42 n = 10 labels have no D >= 3 cell (all the cells are on free ends; no gap cell at D >= 3 at n = 10). Only the unfitted list is new.
- D = 2 failures are real (inside one exact class: 262-row label = one orbit+mirror; 368-row component of the 680 label) and not decided by the nearest relation's length (16 vs 11, 8 vs 9 cells inside one length). Closes my 067 next-item 2. The 680 label = 3 orbit components, meetingPoints depth <= 4 found no link.
- n = 10 source labels here are key+orbit+mirror (176 LNAs unresolved, dropped), not exact classifyLength; images at n = 9 are exact.
- Coxeter key == class at n = 8, 9 (067); no fixed i, nor class+i+coverage rule predicts the image (E-194, E-200).
- Earlier: T4 floating H-020 widths 6..8 pass (0 failures in 8248); T7 right LNA mutation rule exact n = 5..10; E-194 fixed deleted vertex keeps class 37.5-48.9% (baseline 17.6%).

## What I tried
- 070: `workshop/rounds/070/experimentalist_free10.py` (10 s), `_free10_d2.py` (12 s). Reuses `workshop/rounds/027/maverick_classes.py` and the 065 CSV of n = 9 classes.
- 067: `experimentalist_deltable.py`, `_delsplit.py`. 065: `experimentalist_delcompat.py`; classifyLength(9) needs `resume=True`, ~3 runs.
- Lessons: grep EXPERIMENTS.md for E-114/E-146 BEFORE building; S-1 items 1-2 at n <= 10 are recorded. Cached CSVs make a sweep take seconds. Delete JSON above 200 KB.

## What I would do next
1. OVERNIGHT: classifyLength(10) with resume, exact labels, rerun free10.py (removes the key+orbit caveat; places the 176).
2. S-1 item 2 for H-020 head/tail and orbit at n = 12 -> 11 on the lone-3 class (E-146, 2746 LNAs, one orbit): do K >= 3 strips preserve head/tail?
3. n = 11 -> 10 on keys with the D >= 3 filter to cover more of the 32 unfitted labels (size first).
4. Widths 9..11 at w+3..w+4 for T4; rule check at n = 11 (T7).
