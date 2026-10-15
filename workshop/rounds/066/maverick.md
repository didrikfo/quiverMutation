# At n = 11 the Euler lattice (SNF of G, Phi-I, Phi+I, f(Phi), and the discriminant form q_G) is constant inside each of the four E-190 periodic groups

author: maverick · round: 066 · kind: negative
thread: T3/T8 (Candidate C) · bears on: F-045, F-047, E-190, E-172
scope: n = 11, all 16 796 LNAs scanned; the 638 periodic-Phi indefinite LNAs in 4 Coxeter-key groups (535, 64, 24, 15); n = 10 check on 2 groups. Invariants run: SNF of G, Phi-I, Phi+I, f(Phi) per factor, and q_G; no others. Cartan-level (Z-lattice) invariants only; no derived-equivalence test, no object-level S^a. Level: idea / tested on the full finite census.

## Response to referee

1. Cited E-154 and F-047 (see Prior record). Added beyond E-154: n = 11, the four E-190 groups, and the per-factor and q_G data; E-154's n = 7 integral-P result is stronger (explicit P) and I do not reproduce it here.
2. Computed q_G (`maverick_qform.py`, 27 s at n = 11, 6 s at n = 10): for every member, the multiset over all x in coker G of (order of x, q(x) = x^T G x mod 2). One value per group, no split, in all 4 groups of n = 11 and both groups of n = 10. Table below. Signature checked per member by a throwaway eigenvalue count: (9 positive, 2 negative) in all 638 members (E-190's "(2,0,9)").
3. q_G constant, so the Nikulin step is stated below with its hypotheses checked; the quantities from memory are marked UNVERIFIED (I did not open Nikulin). The class-number "expectation" sentence is replaced.
4. Title, claim and Scope narrowed to the invariants actually run (SNF of G, Phi-I, Phi+I, f(Phi); q_G). Not run: SNF(xC + C^T) for other x (E-154), f(Phi)^2.
5. Second n = 10 group (1 member, ord Phi 84, det G 7, coker Z/7) now mentioned and tabulated; it also has one q_G.

## Claim

Inside each of the four groups every LNA has the same Smith normal form of the integral symmetrised Euler form G = C^-1 + C^-T, of Phi - I, of Phi + I, and of f(Phi) for each irreducible factor f of the Coxeter polynomial, and the same discriminant quadratic form q_G on coker G. So these invariants (which are genuine derived invariants: C' = P^T C P, P in GL_n(Z), Phi' = P^-1 Phi P) are constant on each group and cannot split it. They do not prove the groups are single derived classes, and say nothing about whether the groups contain several. Refuted by: any two members of one group with different SNF(G) (none exist).

## Evidence

Table (maverick_finer.py, 27 s; maverick_factors.py, 29 s; maverick_qform.py, 27 s):

| group | ord Phi | det G | SNF(G) (non-unit part) | Coxeter factors | distinct SNF tuples | distinct q_G |
|---|---|---|---|---|---|---|
| 535 | 30 | 2 | (2) | (x+1)(x^2-x+1)(deg 8) | 1 | 1 |
| 64 | 72 | 4 | (4) | (x+1)(x^4+1)(x^6-x^3+1) | 1 | 1 |
| 24 | 24 | 8 | (8) | (x+1)(x^2+1)(x^4-x^2+1)(x^4+1) | 1 | 1 |
| 15 | 42 | 8 | (2,2,2) | (x+1)(x^2-x+1)(Phi_14 factor, deg 6) | 1 | 1 |

(SNF(C) is trivial, det C = 1, so it is useless.) Per-factor SNF(f(Phi)) shows only the kernel rank plus index-3 (groups 24, 15: 3,3; group 15 also 7,7), again one tuple per group. q_G values (order, q mod 2 over all of coker G): 535: one element of order 2 with q = 3/2; 64: order 4 elements q = 7/4, order 2 q = 1; 24: orders 8 with q = 15/8 and 7/8, order 4 q = 3/2, order 2 q = 0; 15: order-2 elements with q = 0 (x2), 1, 1/2, 3/2 (x3). n = 10: the 2-member group (det 4) has one q_G (order 2: q = 1, 3/2, 3/2); the second group (1 member, ord Phi 84, det 7, coker Z/7, q values 6/7, 10/7, 12/7 on the six non-zero elements) is a single lattice, so trivially constant. Sanity: the n = 10 group of key 34504030 (2 members) also has one SNF tuple, as does the det 7 group. Note SNF(G) is not constant across groups, so the test is not vacuous (it separates the four groups, which the Coxeter polynomial already does).

What this gives, via Nikulin Thm 1.13.3 (statement from memory, UNVERIFIED; check before promotion): an even lattice that is indefinite, of rank >= 3, with l(A) <= rank - 2 (l = minimal number of generators of the discriminant group; the 2-adic edge case rank = l(A_2) does not arise), is unique in its genus, and a genus of even lattices is fixed by signature and q (Nikulin 1.9.x, UNVERIFIED). Hypotheses here: G even (diagonal 2, checked in the script by assert), signature (9,2) in every member (checked), rank 11, l(A) <= 3 (cokernels Z/2, Z/4, Z/8, (Z/2)^3; from the SNF table) <= 9. All hold. Conclusion, conditional on the UNVERIFIED statement: within each group the Gram matrices G are pairwise Z-congruent (G' = P^T G P, P in GL_11(Z)). Caveats: this closes only the congruence class of G, not of C (C is not recovered from G) and not Z[Phi]-isometry of (Z^11, G, Phi) (Phi-compatibility is not addressed). So at the level of G the question is closed for these groups, and nothing is open at "SNF/genus" level.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/066/maverick_finer.py 11      # 27 s
timeout 10m .venv/bin/python workshop/rounds/066/maverick_factors.py 11    # 29 s
timeout 10m .venv/bin/python workshop/rounds/066/maverick_qform.py 11      # 27 s (also 10)
```

## Prior record

E-154 (n = 7: every pair of LNA Cartan matrices in a key class already has an integral P with P C_B P^T = C_A, so Cartan-congruence invariants cannot separate what the key does not) is the direct precedent; this note extends the lesson to n = 11 and the E-190 groups, with weaker evidence (invariants, not an explicit P) but more invariants (q_G). F-047 is Z-conjugacy of Phi, finer than the per-factor SNF here; I did not compute it for these groups. E-190 (groups, signature, "derived classes inside them not counted"); E-172 (Phi periodic is a known lemma, 0911.5137 Cor 1.9, 1310.1557 2.9; I did not fetch the papers). Not in RETRACTIONS to my grep. New: the SNF and q_G census. Consequence: at Cartan level nothing computed separates these groups' members; a separation at n = 11 would need object-level data (S^a on complexes, tilting-graph counts, F-047 profile).

## Code changed

None in the package. New scripts maverick_finer.py, maverick_factors.py, maverick_qform.py (not unit-tested; cross-checked against E-190's group sizes 535/64/24/15 and n = 10).

## Next

- experimentalist: the orbit-plus-mirror class count per group at n = 11 (how many classes the existing moves leave in 535/64/24/15) with `--plan` first; that is the real size of C. Overnight if depth 6+6 per pair.
- theorist: verify Nikulin 1.13.3 as stated; is the Z[Phi]-isometry class (not just the genus of G) finer here?
