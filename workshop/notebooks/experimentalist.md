# Experimentalist notebook (rewritten round 115)

## What I now believe
- S-1/H-026: K0 = min(h, n-4-h) for lone-3 classes at key level now fits 10 points: n = 13 h=2,3,4; 15 h=3,4,5; 17 h=5,6; 19 h=5 (K0=5), 6 (6), plus central 19, 21 (7, 8). Round 115 predictions (5, 6, 3) all held. Still key level, head ends, K >= K1 only; no class level, no derivation.
- Two-3 key classes are different: 1, 2 or 3 image keys, no law found (n=13 (2,7): split to K=2, (4,8): 3 keys at K=1, (3,6): none; n=15 (3,9): 3 keys at K=2). Counts 780/265/129/65 do not reappear. Caution: a key class may contain rows not of the two-3 shape.
- E-263 (4,7) "12 keys" = writer bug (R-015), K0 = 4.
- T10/H-015 (110): key guard refuses real steps at n = 6 class 0; open n=7 c1, c2 record 6. T7 (107): no lowering control at k<=4. T6 (105): defect not constant per key class.
- Practical: tail row = shifted row + zeros; image = drop an end entry. Outputs to files; prediction file first. Job size: rows(n-K1) ~ 10 min limit; n=19 K1=3 is ~50 min (overnight proposal if wanted).

## What I tried
- 115: `rounds/115/experimentalist_kloc.py`, `_nonlone.py` (p= positions). 114: `_tail.py`; 113 `_kthresh.py`; 110 `_refusal.py`; 107 `_chain.py`; 105 `_defect_census.py`.

## What I would do next
1. Non-lone: pick p,q small at n = 15/17 with predictions per key count; try to find the rule for depth of split (maybe shorter stretch before the first 3 plus gap).
2. n = 19 h=5 at K1=3 (overnight, ~50 min) to see lower K; n = 23 tail.
3. T10 n=7 class 1 chain; T7 at k=5; T6 depth 6.
