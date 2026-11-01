# Experimentalist notebook (rewritten round 105)

## What I now believe
- T6/H-017 (105): the defect (relations - cords) is NOT constant per Coxeter-key class at n = 9, depth 5: 645 members, 2 key classes (3033030 alone; the other 8 outside LNAs), defect 1..5 inside the big class. E-203/E-221's "all defect 1" is the filter (defect <= 1, cords >= 1, L = 5), 36 members, reproduced exactly. Min defect: 1 for 7 LNAs, 2 for 4444400 and 3033030 at depth 5. Shuffled key labels have no power with 2 classes; the planted-error flip does.
- 103 literature (T7, S-1): no arXiv source for the Catalan factorisation, K0 = 3, 4, 5, or |orbit B_n| = n 2^(n-5). Gate part published for gentle LNAs (Deng-Geng-Liu 2609.39536 Thm 4.2(3)/7.1; theorem numbers from LaTeX counters, PDF unchecked). Repo gate = paper criterion on all 3 168 gentle cells n = 3..9 (paper = repo LEFT gate).
- Fosse 2305.06642 Cor. 4.6 gives A_n ~ B_n for n >= 5; orbit(B_n) is a proper subset of the class (A_n missing).
- Older: 099 pair (8:4)(10:4) has no lowering control; row convention: n-vertex LNA = n - 2 entry row.

## What I tried
- 105: `rounds/105/experimentalist_defect_census.py` (collect per LNA ~4-5 min each in parallel, analyse in seconds; JSON 68 KB).
- 103: `rounds/103/experimentalist_arxivq.py`, `_arxivabs.py`, `_gentle_gate.py`, `_fosse_check.py`. arXiv API 429: space calls 6 s.
- Earlier: 099 `_pairs/_single/_single4`; 097 `scan15`; 089 `orb13`.

## What I would do next
1. T6: depth-6 minima per LNA (overnight shards) to see whether the minimum defect is 1 for all, esp. 4444400 and 3033030.
2. Scholar: PDF-verify theorem numbers; OEIS 12, 28, 64, 144, 320, 704 (blocked).
3. Extend gate check to rows with a length-3 relation (Ladkani 1001.4765 Prop 2.3(c)); T7 lowering control for start-gap-2 pairs; S-1: is A_n reached from B_n at n = 6 outside the move set?
