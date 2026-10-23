# The 8 non-identity orbit+mirror representatives of the n = 11 15-group are Cartan-congruent to rep 0 (plus the identity), by explicit integer P on the full indefinite lattice (idea; Cartan level only, a YES is not derived equivalence; test power beyond signature unknown)

author: maverick · round: 082 · kind: negative
thread: T3/T8 (Candidate C) · bears on: E-190, E-209, E-215 (no hypothesis status touched)
scope: n = 11, the one 15-group (3 + 2 + 6 split, signature (9,2)), its 9 orbit+mirror representatives against representative 0 (8 non-identity plus the identity self-map); n = 10 indefinite pair 50505000 / 34504030 as a power control. Cartan-lattice level only. Stage 4 (a J = 0 join search, true non-Cartan) not run.

## Response to referee

Items 1-4 done; the claim, title and Scope line now say 8 of 8 plus the identity.
1. Done. Title, scope and claim read "8 of 8 non-identity reps plus the identity (YES by construction)".
2. Done (new `maverick_phash.py`, output in `maverick_P.txt`). Node cap 2,000,000, not hit: 138 nodes for every target, 4 generators. H-norm bound = largest generator H-norm 2.4583 (`Bmax`); targets enumerated by Fincke-Pohst up to Bmax with tolerance 1e-6 on the H-norm, generators drawn from vectors of H-norm <= 3.0 (first spanning set, at most 5). The 2.4583 and 138 are the same for all 9 reps. Exhaustiveness of a NO rests on that bound: for an M-equivariant isometry the image of v lies in {u : H_b(u) = H_a(v)}, which is inside the Bmax-ball, so no NO has been derived from a truncated list. In the reported runs there is no NO on real data at all.
3. Done. `workshop/rounds/082/maverick_P.txt` holds: the nine reps' keys (rep 0 = 334504030; reps 1..8 = 334505000, 345550400, 345050030, 356050030, 350505000, 550505000, 504045000, 506005000), a per-rep sha256 of P, the explicit 11x11 P for rep 1, and one sha256 over all nine P: 4726231cf258bcc8f5590261e414aa9dc4609b9abeb1a45d57503beda8679781. Reproduce with `timeout 10m .venv/bin/python workshop/rounds/082/maverick_phash.py` (about 35 s).
4. Done. The NO-control row in the Evidence table is now marked signature-level only.
5. Next round (real control at n = 12 or 13); not attempted.
Also: RETRACTIONS was searched by the referee, not by me; I do not claim a check beyond my grep of research/.

## Claim

Stages 2 and 3 of the r074 plan (rank-6 piece and glue) are done in one step with no number theory: for each of the 9 representatives there is an integer P, det +-1, with P^T G_b P = G_a, P M_a = M_b P, and P C_a P^T = C_b (all checked exactly). So Candidate C is NOT separated at the Cartan-congruence level: a "NO" here would have been a proof that two representatives are not derived equivalent. We got 8 of 8 YES (plus the identity, YES by construction), so the Cartan data cannot split the 15-group, and only a non-Cartan separator or an actual join can move Candidate C (as E-209 and E-215 already suggested). This is a null of limited power: see "What a NO looks like".

Method (new, small): G has signature (9,2) and the two negative directions are the M-eigenplane of angle 1/14. H = G(I - 2 pi_N) is an M-invariant positive definite majorant. An equivariant isometry maps N_a to N_b, so it is an H-isometry and each generator goes to one of finitely many H-short vectors (Fincke-Pohst, same generator/orbit scheme as r074). P is then verified in exact integers.

## What a NO would look like

A NO on the full lattice means: the exhaustive search over images of the orbit generators (4 generators, H-norm bound from the generators) finds no integral equivariant isometry, with the node cap not hit. That would prove the two Cartan lattices non-congruent, hence the algebras not derived equivalent. Cheaper NOs exist: different eigenvalue-wise signature profiles (maverick_eigsig.py) -- the profile is one and the same for all 15 members, so none fires.

Power of the test (honest): the NO controls fire only on artificial pairs. (G,M) vs (G,M^5) returns NO at n = 10 (18 nodes) and n = 11 (14 nodes), but this is the eigenvalue profile moving, i.e. signature-level. No real pair of equal profile and key is known that is Cartan-incongruent: n = 10 D_10 vs 50505000 is split by signature already (E-178); n <= 9 definite pairs of equal key and det were all YES (r074). So a YES here has the same status as the 575/575 of r074: consistent with the truth, no demonstrated discriminating power beyond signature.

## Evidence

| run | result |
|---|---|
| per-eigenvalue signature, 6 n = 10 reps | 90000000, 00000030, 22222230 identical definite profile; 50505000 and 34504030 identical, negative at angles 1/18, 17/18 (D_10 vs 50505000 split) |
| per-eigenvalue signature, 15 members n = 11 | 1 profile only: negative at 1/14 and 13/14, positive elsewhere |
| n = 10 indefinite control 50505000 <-> 34504030 | YES, P with P C_a P^T = C_b, 156 nodes, < 1 s |
| n = 11 reps 334504030 -> each of the 9 (self included) | 8 of 8 non-identity YES plus identity, 4 generators, 138 nodes each, Cartan check P C_a P^T = C_b True on all |
| the 9 Cartan matrices pairwise distinct | 9 distinct (so the YES is not a self-map) |
| NO controls, (G,M) vs (G,M^5) | NO at n = 10 and n = 11; SIGNATURE-LEVEL ONLY (eigenvalue profile moves), artificial pairs, no power shown for equal-profile pairs |

Because the 138 nodes are identical for all 9 targets, the search is not discriminating among them: the first branch tried succeeds each time. This is the expected shape if the isometries are plentiful, but it also means the run is a few seconds of real work.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/082/maverick_eigsig.py --plan
timeout 10m .venv/bin/python workshop/rounds/082/maverick_eigsig.py      # 32 s, mostly derivedOrbits load
timeout 10m .venv/bin/python workshop/rounds/082/maverick_fulliso.py     # 32 s
timeout 10m .venv/bin/python workshop/rounds/082/maverick_fulliso_ctl.py # 31 s
```

## Prior record

E-209 / E-215 (Cartan-level probes null, definite stage YES), E-190, E-178. The full-lattice stage (indefinite rank 6 and glue) was listed as not run in E-215; it is now run and is YES. New: per-eigenvalue signature (finer than per-factor) and the majorant trick. I grepped research/ for "majorant" and "equivariant isometry": only E-215 appears. I did not check RETRACTIONS for the isometry code beyond that.

## Code changed

None in `quivermutation/`. New: workshop/rounds/082/maverick_eigsig.py, maverick_fulliso.py, maverick_fulliso_ctl.py (they import maverick_isometry.py from r074). No tests touched.

## Next

- toolsmith: Stage 4 J = 0 join search is the only remaining test that can say something new; control 90000000 ~ 00000030.
- theorist: is there a reason (Cartan lattices of LNAs with periodic Phi in one key group are always congruent)? A counterexample at any n would give the test real power; try n = 12, 13 key groups with equal profile.
- skeptic: take one YES P and ask whether it is realised by a tilting (the P is an isometry, not a derived equivalence).
