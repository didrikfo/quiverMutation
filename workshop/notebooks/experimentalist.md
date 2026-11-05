# Experimentalist notebook (rewritten round 117)

## What I now believe
- S-1/H-026: all ten (twelve) earlier fits were ODD n. Round 117: at even n the law K0 = min(h, n-4-h) holds for non-central lone-3 classes (n = 10, 12, 14, 16, 10 classes) but the self-mirror class h = n-4-h has NO split at K >= 1 (n = 10, 12, 14; n = 16 K >= 3; n = 18 K >= 5, out of sample, h=6 K0=6 ok). Wording: refutes the even-n EXTENSION, not H-026 as worded; E-146 had n = 12. My pre-registered K0(5,5) = 5 at n = 14 was falsified. H-026 needs rewording; refinement fitted post hoc.
- Rivals min(h, n-5-h), min(h, n-3-h) (short-side reading) equal the law off-centre: no data can separate them there; R1/R4/R5/long-side reading are refuted (power control; planted K0-1 miss caught 8/8).
- Class level (n = 15 h = 3): image key constant on forward-rewrite components (sub-orbits, not full orbit; 1677 components; shuffle null 52 mixed vs 0). Components are sub-orbits only; q_S / derived class untested. My "not separated" prediction was wrong.
- Two-3 classes: no law (round 115). E-263 (4,7) was a writer bug (R-015).
- Practical: kpower.py `table N K1 h,h,..` does a one-pass multi-class table (rows(n-K1) enumeration: 208k rows = 45 s); `orbits` uses lnaMoves.rewritesOf union-find.
- Lesson: when every fit is at one parity (or one side), the rival that differs only at the other parity is the real test; list what the fits cannot separate before choosing rivals.

## What I tried
- 117: `rounds/117/experimentalist_kpower.py` (+ prediction file, outputs *_out.txt). 115: kloc, nonlone. 114 tail; 113 kthresh; 110 refusal; 107 chain; 105 defect census.

## What I would do next
1. Explain the self-mirror exception (head/tail symmetry) and test non-lone self-mirror classes ((p,p) two-3 rows) at even n.
2. Class level properly: q_S or Z[Phi] on rows of the two image-key populations at n = 15 h = 3, K = 3; full orbit with reverse edges (needs rows(15) enumeration, ~12 min, chunk it).
3. Tail end at even n; T10 n=7 class 1 chain; T7 at k=5; T6 depth 6.
