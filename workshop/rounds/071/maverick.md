# n = 11 periodic groups: SNF(xC+C^T), Z-conjugacy of Phi (15-group only, 14 pairs) and (rank, signature, det, theta<=8) of the pieces ker f(Phi^T) are constant in each group (null; low power, n = 10 control run: not blind); the 15-group has 9 orbit+mirror classes (a lower bound on derived classes, not a count); the 16 unplaced LNAs sit in none of the four groups

## Response to referee

Verdict was minor revision; all four required items done.

1. Retitle: done. Title now says "(rank, signature, det, theta<=8) of the pieces ker f(Phi^T)" instead of "per-factor isometry classes", and marks Z-conjugacy as 15-group only (14 pairs against one reference). The indefinite pieces (rank 4/6/8) are compared by (rank, signature, det) only; theta is used where the piece is definite. This is a necessary invariant of the isometry class, not the class. "9 orbit+mirror classes" is a lower bound on the number of derived classes under known moves, not a count of derived classes.
2. The 16 unplaced: done. "Phi not periodic" is cap-limited (order cap 200), not a proof. The placement rests on the Coxeter key (differs from all four group keys) and the signature (9,1,1) against (9,2,0). Claim and Scope reworded.
3. Prior record: done. E-197 already recorded constant per-factor SNF and q_G per group at n = 11 (and listed SNF(xC+C^T) and Z[Phi]-isometry as not run); E-066 already has the orbit+mirror/key refinement (n = 10, 12..16 cores, not n = 11 per group); an n = 7 SNF(xC+C^T) line exists near EXPERIMENTS.md:386. Delta against them: (a) the x-scan at n = 11 over all four groups, which for the 15-group is implied by Z-conjugacy of Phi (SNF(xC+C^T) = SNF(xI - Phi^-1) up to a unit) and so adds content only for the other three groups; (b) 14 explicit integral conjugators in the 15-group; (c) the piece invariant (kernels of cyclotomic factors of Phi^T with restricted G), which goes beyond per-factor SNF/q_G only by signature and theta on definite pieces; (d) the orbit+mirror class counts 9/10/25/213 at n = 11; (e) the placement of the 16.
4. n = 10 control: RUN (`timeout 10m .venv/bin/python workshop/rounds/071/maverick_control10.py 10 8`, 143 s, out `maverick_control10.out`; copy of maverick_pieces.py with the indefinite filter dropped and a guard for degenerate rank-2 pieces at eigenvalue 1 in the smaller groups). In the E-190 group (323 members, factors t+1, t^2-t+1, t^6-t^3+1) there are 2 tuples: 321 members (including 90000000 = D_10 and 00000030, a known equivalent pair that orbit+mirror moves do not join) with pieces (2,(2,0),36,theta (4,2)), (2,(2,0),27,theta (6,6)), (6,(6,0),27,theta 2:18,4:108,6:234,8:234); and 2 members (34504030, 50505000, known equivalent, F-037) with (2,(2,0),36,(6,4)), (2,(2,0),3,...), (6,(4,2),3, indef). So the piece code (i) gives the same tuple on the known equivalent pairs 90000000 ~ 00000030 and 34504030 ~ 50505000, and (ii) separates the known inequivalent pair D_10 vs 50505000 (E-190/E-178), on the definite piece's theta and det as well as on signature. Limit: that pair is already separated by the signature of G ((0,0,10) vs (2,0,8)), so the control shows the code is not blind and does not false-split equivalent pairs, but is no evidence of power beyond signature. The n = 11 null stays low power; it is no longer labelled uncontrolled, but "controlled" means this one known-inequivalent and two known-equivalent pairs at n = 10.


author: maverick · round: 071 · kind: negative
thread: T3/T8 (Candidate C) · bears on: E-190, E-197, E-162, F-047, E-154
scope: n = 11; all four periodic-indefinite Coxeter groups (535, 64, 24, 15; 638 LNAs); Z-conjugacy only on the 15-group (14 pairs against one reference); pieces compared by (rank, signature, det) only where indefinite, with theta<=8 where definite; level: idea. n = 10 control run once (one known-inequivalent pair separated by signature already, two known-equivalent pairs agree). No search, no full isometry decision. 16 unplaced: Phi non-periodicity is order-cap (200) limited; placement rests on signature and Coxeter key.

## Claim

Candidate C (a group is more than one derived class although every Cartan-level invariant agrees) is NOT decided at n = 11. Four further Cartan-level probes were run in a form that could say NO, and all gave the NULL outcome: (1) SNF(xC + C^T) for x = -5..5 and SNF(a Phi + b I) for 7 pairs (a,b): one value per group, all four groups; (2) the 15-group: for each of the 14 other members an explicit P in GL_11(Z) with P Phi_a = Phi_b P, exact check (Z-conjugacy of Phi, F-047 level); (3) for every member and every cyclotomic factor f of the Coxeter polynomial, the invariant (rank, signature, determinant, and theta series up to norm 8 where definite; a necessary invariant of the isometry class of (ker f(Phi^T), G), not the class) of the piece: one value per group. New counts: inside the 15-group the orbit-plus-mirror moves leave 9 classes (6 of two LNAs, 3 singletons); 24-group 10, 64-group 25, 535-group 213. The 16 unplaced LNAs of E-162 lie in none of the four groups (one Coxeter key, a quipu's; signature (9,1,1); Phi not periodic within the order cap 200, which is a cap result, not a proof: the placement rests on signature and key). It does NOT claim each group is one derived class, nor that the groups are Z-isometric as (Z^n, G, Phi): the full isometry question (a P with P^T G_b P = G_a, glueing across pieces, the indefinite rank-6 piece of the 15-group having infinitely many isometries) was not decided.

What would have been a NO: two members of one group with a different SNF(xC+C^T), a different (signature, det, theta) on some ker f(Phi^T), or no integral conjugator of Phi (the last cannot be certified by search, only by the first two). None occurred. The claim is refuted by a single such pair.

## Evidence

Convention checked numerically (repo, `-C^-T C = Phi`): the isometry group of G = C^-1 + C^-T acting on column vectors contains Phi^T (P^T G P = G), not Phi (P G P^T = G). Under C' = P^T C P one has G' = P^-1 G P^-T and Phi' = P^-1 Phi P, so the invariant lattice structure is (Z^n, G, Phi^T). Also (I - Phi) is invertible (no eigenvalue 1) and C^-1 = (I-Phi)^-1 G, so Cartan congruence IS the Phi^T-equivariant isometry of G (no extra data).

| group | LNAs | ord Phi | factors of Coxeter polynomial (ranks) | probes constant | orbit+mirror classes |
|---|---|---|---|---|---|
| 15 | 15 | 42 | (t+1)^3 [rank 3], t^2-t+1 [2], 6th-degree [6] | SNF x=-5..5, 7 (a,b); pieces; Z-conj 14/14 | 9 (sizes 2,2,2,2,2,2,1,1,1) |
| 24 | 24 | 24 | 1, 2, 4, 4 | same, no conj run | 10 |
| 64 | 64 | 72 | 1, 4, 6 | same | 25 |
| 535 | 535 | 30 | 1, 2, 8 | same | 213 |

15-group pieces (one value for all 15): ker(Phi^T+I): rank 3, positive definite, det 168, theta (6,6),(8,4); rank 2: det 3, theta (2,6),(6,6),(8,6); rank 6: signature (4,2), det 7 (indefinite, theta not computed). SNF(2C+C^T) = (1^8, 3, 3, 387) in the 15-group (non-cyclic at 3, unlike the other groups: (1^10, 2979), (1^10, 2907), (1^10, 3315)), x = 0 gives all ones as expected.

Z-conjugacy: Hom lattice {P : P Phi_a = Phi_b P} has rank 17 in all 14 cases; LLL-reduced, a unimodular element with exact check P Phi_a = Phi_b P was found within 8-12 small-coefficient combinations in each case (reference 334504030). Control: a random unimodular conjugate of Phi_a (30 elementary column operations) found at the same rate (9 tries). So F-047 holds in the 15-group by explicit P. The orbit classes are 9, but Phi does not see them.

Mirror and orbit: `freeMoves.derivedOrbits(11, free=True)` plus `freeMoves.mirrorRow` union-find; no orbit leaves its group (class size = members in group for every class printed), a sanity check on the key invariance.

Placement of the 16: Coxeter key shared by all 16 with the quipu P^(1,1)_(1,5,1), key different from all four group keys; C + C^T signature (9,1,1) against (9,2,0) for the groups; Phi not periodic (order cap 200). They are separated from the four groups by two independent Cartan-level facts and need no join test.

Power, honestly: (a) a first run with Phi in place of Phi^T (a mismatch, not a congruence invariant: kernels of f(Phi), Gram of G) took 9 distinct tuples on the 15 members, which shows the piece code is not constant by construction; that the count equals the 9 orbit+mirror classes was NOT examined (it may be a coincidence or a labelling artefact; do not read it as separation). (b) cross-group disjointness of the piece tuples holds, but groups already differ by key. (c) No known inequivalent pair at n = 11 shares key, signature and these invariants, so power against a hard NO is UNESTABLISHED here. The n = 10 control was then run with the filter dropped (Response to referee, item 4, `maverick_control10.py`): equivalent pairs agree, D_10 vs 50505000 differs, but that pair is already split by signature.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/071/maverick_snfx.py 11 5            # 65 s, out: maverick_snfx.out
timeout 10m .venv/bin/python workshop/rounds/071/maverick_orbits.py 11             # 62 s, out: maverick_orbits.out
timeout 10m .venv/bin/python -u workshop/rounds/071/maverick_conj.py 11 15 3000    # 36 s, out: maverick_conj.out
timeout 10m .venv/bin/python -u workshop/rounds/071/maverick_pieces.py 11 8        # 60 s, out: maverick_pieces.out
timeout 10m .venv/bin/python workshop/rounds/071/maverick_place16.py               # 55 s, out: maverick_place16.out
timeout 10m .venv/bin/python workshop/rounds/071/maverick_control10.py 10 8        # 143 s, out: maverick_control10.out (n = 10 control, indefinite filter dropped)
```
(`maverick_conj.py` is imported by `maverick_pieces.py` for `refine`. The 16 names come from `workshop/rounds/050/maverick_sixteen.py 11`, 10 min.)

## Prior record

E-197 (per-factor SNF, q_G constant per group at n = 11; SNF(xC+C^T) and Z[Phi]-isometry listed as not run), E-066 (orbit+mirror vs key refinement, n = 10, 12..16), E-190 (four groups, C open at n = 11), E-162 (the 16), F-047 / E-154 (Z-conjugacy; n = 7 explicit P). Grepped `research/` for "ker f(Phi)", "theta series", "orbit+mirror classes" at n = 11 per group: no hit. New: the x-scan, 14 explicit conjugators, the piece invariant, the 9/10/25/213 class counts, the placement of the 16. Note SNF(xC+C^T) = SNF(x I - Phi^-1) up to a unit (C^T = -C Phi^-1), so it is a function of the Z[t]-module (Z^n, Phi) only: it adds nothing beyond F-047 and, since 14 of 14 Phi are Z-conjugate, its constancy is implied in the 15-group (checked independently for all four groups). The piece invariant is the only probe here that uses G and Phi together.

## Code changed

None under `quivermutation/`. New scripts in `workshop/rounds/071/`: maverick_control10.py, maverick_snfx.py, maverick_orbits.py, maverick_conj.py, maverick_pieces.py, maverick_place16.py. No tests run (no repo file touched).

## Next

Not a push to continue Cartan-level probes on these groups: three levels agree and the remaining Cartan-level question is a class-number statement (genus of Hermitian lattices over the Z[Phi]-order; indefinite and rank above 2, so strong approximation suggests one class per genus, stated from memory, UNVERIFIED). A theorist could decide that, or I can attempt the staged explicit isometry (enumerate isometries on the definite pieces, glue, test integrality) on the 15-group to get a Cartan congruence for pairs of different orbit classes (E-154 style); that would still be Cartan-level. To separate classes one needs a non-Cartan invariant or a join: toolsmith: J = 0 join search between two of the 9 orbit+mirror classes of the 15-group (nine n = 11 classes; 90000000 ~ 00000030 at n = 10 is the long-path control for the join search, E-190).
