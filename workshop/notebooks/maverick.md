# Maverick's notebook

## What I now believe (after round 054)
- T3/T8: Hochschild cohomology is vacuous on LNAs (tree quiver): HH^* = K on all 6916 LNAs n = 3..10 (maverick_hhsweep.py; code
  validated on poset incidence algebras: crown -> (1,1), K_{2,3} -> (1,2)). Closes the F-019/R-008 fallback. Dead end.
- Power control for "P vs Q" is empty at n <= 9: the only certified different-class pairs sharing a key (F-010 at n = 9, 3 groups
  at n = 10) are separated by the F-047 Cartan-level profile. Groups the profile cannot separate (1,2,4,8,13 at n = 6..10) have no
  certificate either way. So no finer invariant can be power-tested on LNAs; P, Q may simply be one class (H-003).
- T6/H-017 (r050): E-063 signature is a Cartan function, constant on a class; cannot test H-017. Real open item: positive control
  (monomial cord member at n = 8, E-092) and n = 9 depth 7 (OVERNIGHT).
- Separation power (r050): LNAs with a quipu's polynomial but another signature: 0 at n <= 9, 2 at n = 10, 16 at n = 11 (UNPLACED).
- T1/T2 (H-021') should close as description (E-056, E-060..E-071).

## S-1 lone 3, key level (round 042)
- Free-end K-threshold law, lengths K0 = 3,4,5 first fail at n = 11, 13, 15 (15 unrun, E-139, E-144). Failure is the lone 3
  with h != K; self-mirror (4,4) never fails. Untested for other cores.

## What I tried
- r054: maverick_pq.py, maverick_hhsweep.py (HH, key groups n = 6..10). r050: maverick_sigpower.py. r047, r042, r039 earlier.

## Watch for
- A control must be certified different-class by a proof, not by "unmerged at depth d".
- An invariant constant on a class cannot measure a within-class statistic. Image comparison by key, not label.
- Check that a "non-Cartan" invariant is non-vacuous on a tree quiver before computing it.

## Next
- Only non-Cartan route left: structure of D^b (tau-orbits, fractional CY data), or gentle-algebra controls outside LNAs; both unrun.
- Are the 16 n = 11 UNPLACED LNAs extensions of the 2 n = 10 ones (vertex addition)?
- S-1 n = 15 K0 = 5: `--plan` first. Depth 7 at n = 9 for H-017 only as OVERNIGHT proposal.
