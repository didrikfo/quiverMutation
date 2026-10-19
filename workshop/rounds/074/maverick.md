# Staged Cartan-congruence test, definite stage only: control 90000000 ~ 00000030 gets an explicit P (Cartan congruence checked directly); 15 of 15 members of the n = 11 15-group pass on the two definite pieces (idea; NO power unestablished; the indefinite rank-6 piece and the glue are not run)

author: maverick · round: 074 · kind: tool + null
thread: T3/T8 (Candidate C) · bears on: E-190, E-209, E-178 (no hypothesis status touched)
scope: n = 10 definite LNAs with periodic Phi (577 definite, groups 256 and 321; 575 = 320 + 255 members tested against their group rep) and the n = 11 15-group, definite pieces only (ranks 3 and 2 of 3 + 2 + 6); idea level; nothing about derived equivalence is proved. The NO branch has NOT been shown to separate any two Cartan lattices of equal key and signature (none exists in the searched range), so "575 of 575 YES" is a null of unestablished power, not evidence for "one class".

## Response to referee

New script `maverick_controls.py` (parts 1, 2, 3; whole run 63 s). Verdict was minor revision; all four required items answered.

1. P convention and direct Cartan check: DONE. The code returns P with P^T G_b P = G_a and P M_a = M_b P. The Cartan relation in this convention is **P C_a P^T = C_b** (i.e. Q = P^T gives Q^T C_a Q = C_b); the line in "Plan" below, "P^T C_b P = C_a", had the direction reversed and is corrected here. Part 1 of the script checks `P C_a P^T == C_b` on every YES of the n = 10 survey (member vs rep): 575 of 575 hold, and in all 575 the other direction `P^T C_b P == C_a` fails, confirming the transpose convention. So the isometry-to-Cartan-congruence step is verified exactly on all YES outputs of the n = 10 survey (not yet added to pieces15, where only the definite pieces, not C, are in play).
2. (G, -M) NO control: DONE as a runnable command (`maverick_controls.py 2`): NO after 180 nodes. Labelled artificial and non-Cartan: -M is not the transpose Coxeter matrix of any Cartan matrix. It shows the search can fail, nothing more.
3. Scope: DONE (see Scope line and title). Attempted a real NO-power control (part 3): at n = 4..9, all definite LNAs grouped by (Coxeter key, det G), each member vs the first: 3, 12, 38, 112, 259, 271 pairs, all YES, zero NO (and each YES passes the direct Cartan check). So there is no pair of definite Cartan data of equal key and det in n <= 9 for the test to separate; the real-NO control cannot be built from this sample. Power beyond key/signature remains undemonstrated.
4. Stage 2 stays marked unverified (class-number argument from memory; not to be cited as established in FINDINGS-bound text). Added: the generator count and bounds are in the code (kmax = 4, norm bound 4 for n > 4, node cap 300000, never hit); the survey records P and checks it exactly via isometry().

Also narrowed: "a NO at any stage is a proof of derived inequivalence" holds for a NO on the full lattice, and also for a NO on any single Phi-invariant piece (the pieces are canonical sublattices); a YES on pieces says nothing about the glue or the full lattice. Not done: the pieces15 rerun with direct Cartan check (pieces are not Cartan matrices; N/A), and any indefinite stage (out of this sitting).

## Plan (what the stages are, and what each can say)

Convention (E-209): G = C^-1 + C^-T is invariant under M = Phi^T. Cartan congruence is equivalent to an M-equivariant isometry P: (Z^n, G_a, M_a) -> (Z^n, G_b, M_b) (since C^-1 = (I-Phi)^-1 G); code convention: P C_a P^T = C_b (Q = P^T gives Q^T C_a Q = C_b). A derived equivalence forces such a P (Cartan matrices of derived-equivalent algebras are congruent), so a NO on the full lattice (or on a canonical Phi-invariant piece) is a proof of derived inequivalence; a YES is only "not separated".

- Stage 0 (done, E-197, E-209): SNF, q_G, per-factor (rank, sig, det, theta), Z-conjugacy of Phi (15-group). All null.
- Stage 1 (this round, exact): definite lattices. Pick generators v_1..v_k of least norm whose M-orbits span Q^n; every equivariant isometry sends v_i to a vector u_i of equal norm and equal G(v_i, M^c v_j); enumerate all such tuples (Fincke-Pohst), set P = W(u) W(v)^-1, test integrality, det +-1, P^T G_b P = G_a, P M_a = M_b P in exact arithmetic. Exhaustive, so NO is a proof (within the node cap; the script reports INCONCLUSIVE-cap otherwise).
- Stage 2 (NOT run; proposal): the rank-6 piece L_6 = ker(Phi_7-type factor) of the 15-group is indefinite, signature (4,2), det 7, a Z[zeta_7]-lattice of rank 1. Q(zeta_7) has class number 1, so L_6 = (O, a x conj(y)) with a in the real subfield O^+, and isometry of two such pieces is "a_b / a_a is a norm of a unit of O" -- a finite computation in the rank-2 unit group of Q(zeta_7)^+ (totally positive units mod squares). I have not carried it out; the claim that the class-number argument applies as stated is from memory and should be checked by a theorist or the scholar.
- Stage 3 (NOT run): glue. Z^11 contains L_3 + L_2 + L_6 with finite index; two members are isometric iff the piece isometries can be chosen compatibly with the glue (a finite group computation in the discriminant groups, orders 7 x 3 x ... small). This is where a NO could first appear for the 15-group if pieces all pass.
- Stage 4 (true non-Cartan, toolsmith): the J = 0 join search between orbit+mirror classes, power-controlled by 90000000 ~ 00000030 (E-190, E-210: all 10 deletions land in one n = 9 class, so no deletion invariant either). Stages 1-3 are Cartan-level in disguise: a stronger form of the Cartan test, not a new kind of invariant.

## What stage 1 can say NO to, and what it cannot

Can: any pair of definite (or definite-piece) Phi-lattices not related by an equivariant isometry, e.g. (G, M) vs (G, -M) at n = 10: NO after 180 candidate images (artificial control, shows the NO branch fires; maverick_isometry.py is the library, the one-liner is in Reproduction).
Cannot: (i) anything about indefinite lattices (SKIP), so not the full 15-group nor 50505000/34504030 at n = 10; (ii) inequivalence of pairs that ARE Cartan-congruent. The only inequivalent n = 10 pair known (D_10 vs 50505000, E-178) is already split by signature, and the only definite periodic n = 10 groups (256, 321) are single Cartan-congruence classes (below), so the test has NO demonstrated power at the Cartan level beyond signature at n = 10; its power at n = 11 is unestablished because no known inequivalent n = 11 pair shares key and signature (E-209). Galois twists (G, M) vs (G, M^k), k coprime to 18, are all YES: the test does not manufacture NOs.

## Evidence

| run | result |
|---|---|
| 90000000 vs 00000030 (n = 10, power control, equivalent pair) | YES: explicit integer P, det +-1, 360 candidate tuples, 0.2 s (also reversed and self) |
| all 321 members of the D_10 key group vs 22222230; all 256 of the A_10 group vs 00000000 | 575 of 575 YES (maverick_survey10.py, 43 s) |
| (G, M) vs (G, M^k) k in {5,7,11,13,17}, reps 90000000, 22222230 | 12 of 12 YES (maverick_negctl.py) |
| (G, M) vs (G, -M), 90000000 | NO (180 nodes), artificial NO control |
| 15-group, piece t+1 (rank 3) and piece t^2-t+1 (rank 2), all 15 members vs member 0 (hence all 9 orbit+mirror reps) | 15 of 15 YES on both pieces (56 s, mostly the derivedOrbits load) |

Reading: at the definite stage the 15-group is not separated; this is the Cartan-level null of E-209 repeated with explicit isometries instead of theta counts on the definite pieces (YES there is a strengthening: it is an equivariant isometry, not just equal invariants). It says nothing about the rank-6 piece, where the 9 classes could still differ. The n = 10 control confirms that "equivalent but not joined by orbit+mirror moves" passes with an explicit P, i.e. the test is consistent with known truth; it is not a certificate of equivalence.

Limit of the "complete search" claim: the generator choice requires the orbit span to be full (checked; kmax = 4 generators); the node cap 300000 was never hit.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/074/maverick_isometry.py 10        # 90000000 vs 00000030, <1 s
timeout 10m .venv/bin/python workshop/rounds/074/maverick_survey10.py 10        # 43 s
timeout 10m .venv/bin/python workshop/rounds/074/maverick_negctl.py             # 2 s
timeout 10m .venv/bin/python workshop/rounds/074/maverick_pieces15.py           # 56 s
```
```
timeout 10m .venv/bin/python workshop/rounds/074/maverick_controls.py   # parts: direct Cartan check, (G,-M) NO control, n<=9 NO search; 63 s
```

## Prior record

E-209 / E-197 / E-190: Cartan-level probes null; Z-conjugacy of Phi in the 15-group explicit; pieces compared by theta. New here: explicit equivariant isometries (pieces, and the full D_10 key group at n = 10), and the stage plan. Not new: the D_10 key group being one derived class (E-190, E-178 level). I did not grep RETRACTIONS for isometry-search code beyond the E-numbers above.

## Code changed

None in `quivermutation/`. New: workshop/rounds/074/maverick_isometry.py (library + n = 10 control), maverick_survey10.py, maverick_negctl.py, maverick_pieces15.py. No tests touched.

## Next

- theorist/scholar: confirm the stage-2 reduction (rank-1 hermitian lattices over Z[zeta_7], class number 1, unit-norm criterion), then toolsmith or I implement it; stage 3 glue after.
- toolsmith: stage 4 join search with 90000000 ~ 00000030 as the long-path control (unchanged request).
- Honest outlook: if stages 2-3 also pass, Candidate C stays undecided at n = 11 and only a non-Cartan separator or a join can move it.
