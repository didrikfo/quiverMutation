# No derived-class invariant bounds d_i; the one published formula that sees J_i is a Cartan-matrix defect, and it explains why gate-admitted J_i != 0 with d_i >= 3 contradicts nothing

author: scholar · round: 038 · kind: negative (with one small derivation)
thread: T5 · bears on: E-128, E-130, E-131, E-133, E-134, H-015

## Claim

(1) d_i = dim e_iAe_v = (C_A)_{iv} is an entry of the Cartan matrix in a chosen vertex order. It is not a derived-class invariant: derived equivalence only fixes C_A up to Z-congruence (Ladkani 1001.4765 Prop 3.6 / Lemma 3.5 context; Happel-Seidel and math/0610685 Cor 3.13 in the local notes), and a congruence does not preserve single entries. Likewise "number of arrows out of v" is a property of the pair (A, v) (it is the row r^+_v of Ladkani Prop 3.6). I found no published statement, in the 25 local notes or in AI 1009.3370 as summarised, that bounds d_i or the out-degree by a derived invariant. So the answer to "which class invariant governs d_i" is: none known; the observed d_i = 2 at J_i != 0 (E-131) is a statement about the walk, not about the derived class. UNVERIFIED: the printed AI text (arXiv blocked; I read only `research/literature/1009.3370-silting-mutation.md`).

(2) AI Thm 2.32(b) is an iff with an injectivity condition and no bound on dimensions; Thm 2.31 (local note's numbering; the assignment says "Prop 2.31", UNVERIFIED which is printed) says the mutation is silting regardless. So J_i != 0 with d_i >= 3 is not contradicted by anything in the paper: it just says the mutation is silting, not tilting. E-134's hand example (d, dim J) = (3, 1) and E-128's layered (3, 2) are exactly such cases; E-128's L1 (dim J_i <= d_i - 1) is the only constraint, and it is the gate's, not the paper's.

(3) New and checkable: for acyclic A, loopless v, with B = End(mutated T) and r = r^+_v of Ladkani Prop 3.6,
C_B = r C_A r^T + H, with H_{v i} = dim J_i and all other entries 0.
Reason (derivation, stated fully): the class of each summand in K_0 gives the Euler form of T as r C_A r^T (the Lemma-3.5 computation needs only [T_i] = sum r_ij [P_j], so it holds for silting T). T is 2-term silting, so Hom(T,T[k]) = 0 for k >= 1 (AI 2.31) and k <= -2; hence Cartan(B)_{ab} = Euler_{ab} + dim Hom(T_a, T_b[-1]). By E-130 (a vertex-acyclic case) Hom(T_v, P_i[-1]) = J_i and Hom(T, T[-1]) has no other part. So the defect H is the dimension of Hom(T,T[-1]) placed in row v. Not a new invariant: it is E-130 read in K_0. Checked on three algebras (below); the identity is verified there, not proved beyond the derivation above.

Consequence: a gate-admitted step with J != 0 yields B with Cartan matrix r C r^T + H. B is derived equivalent to A only if there happens to be an unrelated equivalence; the Coxeter guard (key) then compares the Coxeter polynomial of rCr^T + H with that of A. The invariant that decides acceptance of such a row on a walk is therefore the Coxeter polynomial / Z-congruence class of C_B, which is where d_i enters: it appears through H = dim J_i <= d_i - 1 (E-128) and through r. It does not bound d_i.

## Evidence

Script `scholar_cartan_defect.py`: C_B - r C_A r^T for three gate-admitted J != 0 algebras (convention: r row v as in the note, `invariants.cartanMatrix` unchanged).

| algebra | d_i, dim J_i at the J != 0 vertex | defect matrix (C_B - r C_A r^T) |
|---|---|---|
| E-134 T1 (1 -> {2,3,4} -> 5 -> 6, one commutation, v = 5) | (3, 1) at i = 1 | single entry (row v, col 1) = 1 |
| T2 (same, two commutations, m = 3) | (3, 2) at i = 1 | single entry (row v, col 1) = 2 |
| E-080 square (v = d) | (2, 1) at i = a | single entry (row v, col a) = 1 |

In each, the defect equals dim J_i at (v, i) and is 0 elsewhere (the transposed reading r^T gives a non-sparse defect, so the orientation is as stated). Three cases, small: supports the derivation, does not prove the cyclic case (there H_vv can be nonzero, E-130).

What the literature says about "d_i small": nothing found. Local notes grepped: 1009.3370, 1001.4765, 2112.08129, 2509.12983, 1504.02617, 1305.5213, 0805.1018. Strong global dimension (1305.5213) is a derived invariant, finite iff piecewise hereditary, but it is a property of B only when B is a tilting End (J = 0), so it cannot constrain a J != 0 row.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/038/scholar_cartan_defect.py   # a few seconds
```

## Prior record

E-130 (Hom(N,N[-1]), silting-not-tilting iff J != 0), E-128 (L1), E-131/E-133/E-134 (d_i data). E-134 already says "d = 2 stays empirical and is a statement about derived classes, not the gate"; I sharpen it: it is not even a statement about derived classes, because d_i is a Cartan entry, not a class invariant (point 1). The defect formula in (3) is not in `research/` (grep "defect", "r C r", "Prop 3.6" in EXPERIMENTS/FINDINGS: only the 171-mutation check of C = r C r^T on tilting steps, 1001.4765 note). Not in RETRACTIONS.

## Code changed

None (new script only: `workshop/rounds/038/scholar_cartan_defect.py`; no tests needed).

## Next

- theorist: use (3) to turn "J_i != 0 child passes the key guard" into a Cartan statement: det/char-poly of C_B^{-T} C_B with H = dim J_i e_v e_i^T; a J != 0 row on a walk is then a key coincidence by construction, which would explain E-134's 167 / 16 split and why E-131 sees only tilting-compatible rows at d = 2. Check on the E-133 rows (d_i = 4, 4, 5): is their child's key equal to the parent's, and is the child derived equivalent? (If E-131 counts rows whose child fails the guard, they are not "on the walk".)
- toolsmith: on the E-131 rows with J != 0, report whether the child passes the guard.
- anyone with network: confirm AI numbering (2.31 as Prop or Thm) and that 2.32(b) as quoted is the printed statement.
