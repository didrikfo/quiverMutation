# Experimentalist notebook (rewritten round 075)

## What I now believe
- Row convention (bit me in 073): an n-vertex LNA is an n - 2 entry row (`allRelationLengths(n)`, `lines.className`); E-210 wrote the n = 9 images of the E-190 pair as 80000000 / 00000300 (8 digits = n = 10). The n = 9 pair is 8000000 / 0000030. Check digit count against n before trusting a name.
- 075: E-190 phenomenon holds at n = 9 and the whole family A_n = rad^(n-1), B_n = 3-relation on the last four vertices, n = 6..12: same key, A_n has no move (orbit 1), B_n orbit n*2^(n-5), never joined (joined only at n = 5). Deletion: A_n -> A_{n-1}; B_n -> B_{n-1} for i <= n-3, hereditary 0..02 otherwise. Join/non-join is preserved, orbit size follows a recurrence. Class equality only exact at n = 9, 10.
- Free-stretch rule (E-114/E-146/E-200/E-206): D >= 3 from every relation support gives one image class (n = 9 -> 8 6/6, n = 10 -> 9 10/10). Mostly vacuous (32 of 42 n = 10 labels have no D >= 3 cell). The E-190 pair is outside its domain.
- Exact classifyLength(10): seed places 61.0 % (311 s); main pass unsized (<= ~3 h arithmetic, resolve step unknown). classify.py has --budget-hours, exits 2 on budget. Proposal sits in 073 submission.
- Coxeter key == class at n = 8, 9; no fixed i or class+i+coverage rule predicts the image (E-194, E-200). T4 widths 6..8 pass; T7 right LNA rule exact n = 5..10.

## What I tried
- 075: `rounds/075/experimentalist_join9.py` (orbitReport, movesJoin, derivedOrbits + mirrorRow), `experimentalist_delpair.py` (family n = 5..12, 10 s).
- 073: `experimentalist_pair.py`, `experimentalist_plan10.py`. 070: `_free10.py`, `_free10_d2.py`, `_pool.py`. 067: `_deltable.py`, `_delsplit.py`. 065: `_delcompat.py`.
- Lessons: grep EXPERIMENTS.md first; mirror for rows of another n must use len(row)+2 (edge moves can change n in an orbit); forward orbit absent != not equivalent; background jobs piped through tail show nothing until the end.

## What I would do next
1. Overnight --run loop for exact classifyLength(10), then redo free10 on exact labels (check the --run path on n = 8 first).
2. S-1 item 3: an LNA with a long core equivalent to a quipu/almost-separate LNA, delete outside the core; the A_n/B_n family has no core, so it says nothing.
3. S-1 item 2 for H-020 head/tail at n = 12 -> 11 on the lone-3 class (E-146).
4. A join search for `8000000` ~ `0000030` (moveless vs mobile row) as the long-path control; ask toolsmith.
5. n = 11 -> 10 with the D >= 3 filter (size first).
