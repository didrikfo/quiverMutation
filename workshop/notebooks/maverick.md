# Maverick's notebook

## What I now believe (after round 050)
- T6/H-017: the E-063 signature (of C+C^T) is congruent to the Euler form's, so it is a Cartan-matrix function and a derived-class invariant. It is constant on a class, hence blind to (cords, relations) of members. It cannot test H-017; it only separates classes. H-017 stays OPEN. The real open item is a positive control for the search (monomial cord member at n = 8, E-092) and n = 9 depth 7 (OVERNIGHT).
- Separation power (maverick_sigpower.py, r050): LNAs with a quipu's polynomial but a different signature: 0 at n <= 9, 2 at n = 10 (34504030, 50505000), 16 at n = 11 (all UNPLACED). F-010 cospectral pair at n = 9: same signature and Smith form, not separated.
- T1/T2 (H-021', s = n - k(c)) should close as description (E-056, E-060..E-071); pairing is a property of the forward reduced walk (E-065).

## S-1 lone 3, key level (round 042)
- Free-end K-threshold law, three lengths: K0 = 3, 4, 5 first fail at n = 11, 13, 15 (15 unrun at class level, E-139, E-144).
- Failure is the lone 3 with h != K; self-mirror (4,4) never fails. Deletion from the shorter free side; untested for other cores.

## What I tried
- r050: maverick_sigpower.py (n = 6..11 signature tables). r047: note on E-059..E-071. r042: maverick_n12.py. r039: n = 13 class.

## Watch for
- Image comparisons by key, not label: "differs" sound, "same" not.
- An invariant constant on a class cannot measure a within-class statistic; check that before designing a power test.
- One class per length is a small sample.

## Next
- Are the 16 n = 11 UNPLACED LNAs extensions of the 2 n = 10 ones (vertex addition)?
- S-1 n = 15 K0 = 5: `--plan` first. Non-lone cores.
- Depth 7 at n = 9 for H-017 only as OVERNIGHT proposal.
