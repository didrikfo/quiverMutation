# Maverick's notebook

## What I now believe (after round 055)
- T3/T8: HH^* = k on all 6916 LNAs n <= 10, a THEOREM already on file (2312.14699 note, 0805.1018 Prop 5.1, EXPERIMENTS ~l.2084, R-008).
  My 054 "nobody computed it" was wrong; it is a code check only. Code now has relation-bearing controls (rad^2 = 0 crown+sink gives (1,2),
  the incidence version (1,), the LNA version (1,)); so all-(1,) is not a relation bug. HH_* is also dead on LNAs (argued, not run).
- "Key" = Coxeter polynomial. At n = 10: 40 key groups; F-047 orbits (tables+free+edges+doubles) 113, 25 groups with >= 2 orbits; plus
  mirror join 71 classes, 16 groups; 3 split by profile, 13 unresolved. Same objects, different unit (maverick_recon.py).
- Power control: no certified-inequivalent pair with equal key and equal F-047 profile at n <= 9. Unresolved classes are UNCERTIFIED, so
  "distinct class" never means inequivalent. Certified pairs with different profile remain for agreement testing (F-010; 3 groups at n = 10).
- Candidate C (idea): object-level Serre periodicity S^a(A) ~ A[b]. Vacuous if Phi has infinite order: true for the F-010 pair and 2 of
  3 certified n = 10 groups. Live at one n = 10 certified group (4 classes, Phi^18 = I exactly; matrix-checked, char poly (T+1)^2(T^2-T+1)(T^6-T^3+1),
  -1 eigenspace dim 2 so diagonalisable there). Unimplemented (needs minimal complexes).
  The periodicity obstruction is KNOWN (0911.5137 Cor 1.9, 1310.1557 2.9); only the finite-order-Phi x certified cross-table is mine. Entropy = spectral radius: unchecked folklore, dropped.
- T6/H-017: E-065 signature is Cartan-level; real open items are the monomial cord positive control (E-094) and n = 9 depth 7 (OVERNIGHT).
- Separation power (r050): LNAs with a quipu's polynomial but another signature: 0 at n <= 9, 2 at n = 10, 16 at n = 11 (UNPLACED).

## S-1 lone 3, key level (round 042)
- Free-end K-threshold law, lengths K0 = 3,4,5 first fail at n = 11, 13, 15 (15 unrun, E-141, E-146). Failure is the lone 3 with
  h != K; self-mirror (4,4) never fails. Untested for other cores.

## What I tried
- r055: maverick_control.py, maverick_recon.py, maverick_fcy.py, maverick_phiorder.py. r054: maverick_pq.py, maverick_hhsweep.py.
  r050: maverick_sigpower.py. r047, r042, r039 earlier.

## Watch for
- grep research/ AND research/literature/ (fractional, Calabi, periodic, Serre too; r055 review caught me) for the closing theorem before announcing a "dead end" or "nobody computed".
- A positive control must exercise the feature under test (relations), not just the code path (posets).
- A control must be certified different-class by a proof, not "unmerged at depth d".
- An invariant constant on a class cannot measure a within-class statistic; image comparison by key, not label.
- Denominators: always say what the unit is (orbit, class, polynomial group) when quoting counts against F-047.

## Next
- Candidate C at the one live n = 10 group: implement S on complexes of projectives over an LNA; scholar to check the fractional-CY literature.
- Gentle-algebra controls outside LNAs (equal Cartan, different AG) remain unrun. 16 n = 11 UNPLACED LNAs vs the 2 at n = 10 (vertex addition?).
- S-1 n = 15 K0 = 5: `--plan` first. Depth 7 at n = 9 for H-017 only as OVERNIGHT proposal.
