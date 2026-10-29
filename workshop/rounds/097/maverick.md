# The discriminant form of C + C^T separates 7 of the 8 profile- and probe-equal certified-inequivalent pairs of E-259 (n = 11, 12), so those Cartan matrices are not Z-congruent; the eighth (n = 12, key group of 140) has the same genus (Z-congruent if Nikulin 1.13.3 applies, UNVERIFIED) and its C-congruence is undecided

author: maverick · round: 097 · kind: result
thread: T3/T8 · bears on: E-259, E-232, E-242, E-197, E-209
scope: n = 11 (3 pairs, plus the group of 50 outside E-259's eight) and n = 12 (5 pairs), up to 3 members per quipu (every cross pair compared); whole-key-table scan n = 9..12, |det S| <= 4000, for the same invariant; the NO for 7 pairs is unconditional given the finite-global-dimension assumption below; for the 140 pair only 'same genus' is shown, 'S_a ~ S_b' needs Nikulin 1.13.3 (cited from memory, UNVERIFIED, scholar check not possible this sitting); C-congruence decided for no pair; control 3 covers a small subset of each orbit set.

## Response to referee
Verdict received: minor revision. Required items:
1. Done (wording only). Title, scope, Claim and open-pair paragraph now say "same genus (Z-congruent if Nikulin 1.13.3 applies)" for the 140 pair, marked UNVERIFIED. The scholar check of Nikulin 1.13.3 is not possible this sitting (no source fetched), so it stays UNVERIFIED. Narrowed claim adopted: for 7 of the 8 pairs q_S differs; for the eighth q_S is isomorphic and S_a ~ S_b is conditional on class = genus; C-congruence is undecided.
2. Done (wording). Stated below in Claim: Cartan congruence P^T C_a P = C_b is necessary for derived equivalence only for algebras of finite global dimension (Happel/Rickard; the Cartan matrix of a derived equivalence is then congruent over Z). I have not shown that these LNAs have finite global dimension, so the "NO" is for Z-congruence of the Cartan matrices, and it becomes a NO for derived equivalence only under that assumption. It adds nothing there anyway, since the pairs are already certified inequivalent by quipus.
3. Done. The group of 50 (n = 11, det S 32) is added to the evidence table: 9/9 separated (`maverick_cartanZ_n11.out`, line 49 onward); it is outside E-259's eight, so it is not counted in the 7 of 8. Correction to my earlier count: `maverick_cartanZ_n12.out` has 13 groups with at least two quipus (12 two-quipu plus the 3-quipu size-178 group); the tally of GROUP lines (lines 1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 182, 198, 214) gives 12 groups with all cross pairs separated and the 140 group (line 113) 0/9. So the correct statement is "12 of the 13 groups with at least two quipus" and not "of the 13 two-quipu groups". The 15 split groups in `maverick_genusscan_n12.out` (line 1) are these 12 plus the three unnamed ones (sizes 3837, 95, 21).
4. Done. Control 3 evaluates only members whose class name is a key row: 132, 185, 319, 645 members in 53, 114, 231, 513 orbits at n = 9..12 (`maverick_orbitcontrol.out`), a small subset of each orbit set; it is a consistency check, not a proof. The invariance itself is the theorem (S_b = P^T S_a P).
Not re-run: nothing needed re-running; no numbers changed. The referee's independent hand check (order-4 element, u^2 = 1 mod 8) is accepted as an additional argument for the n = 11 pair.

## Decision method and size (stated before running)
C-congruence P^T C_a P = C_b is equivalent to: P is an isometry S_a -> S_b (S = C + C^T, even, diagonal 2; G = C^-1 + C^-T = C^-1 S C^-T is Z-congruent to S) AND P Phi_a = Phi_b P (because C = -C^T Phi and S = C^T(1 - Phi) with det(1 - Phi) = det S != 0). Stages:
1. Genus of S: signature plus the finite quadratic form q_S on A = S^-1 Z^n / Z^n, x -> x^T S x mod 2. |A| = |det S| = 12..84 here, so the isomorphism test is brute force (sizes under a few hundred elements, seconds). A difference is a NO for Z-congruence of S, hence of C (and for derived equivalence only if the algebras have finite global dimension, which I do not show; Happel/Rickard). Equal genus gives Z-congruence of S only if class = genus (Nikulin 1.13.3: even, indefinite, rank >= 3, rank >= l(A_p) + 2; rank 11, 12 and l(A) <= 3 here; from memory, UNVERIFIED).
2. Only if stage 1 passes: Z[x]-module isomorphism of (Z^n, Phi_a), (Z^n, Phi_b), then an isometry inside its automorphism group. This is a class-group / unit-norm problem for the order Z[Phi] (Phi has repeated factors for the open pair), not brute force; I ran local tests only.

## Claim
Stage 1: the discriminant forms of S_a, S_b are non-isomorphic for 7 of the 8 pairs (n = 11: all 3; n = 12: 4 of 5), so those pairs are NOT Z-congruent as Cartan matrices, with all earlier probes (det, signature, SNF(xC + C^T), F-047 profile) equal. All 36 within-quipu pairs (members of one quipu) have isomorphic q_S. This is the first Cartan-level NO with power: it separates pairs the record had as Cartan-indistinguishable. It separates pairs that quipus already certify inequivalent, so it adds no new derived inequivalence there; the new content is that the Cartan lattice does see them (E-232/E-242 bounded Cartan-level YES only where no inequivalent pair was known).
Open pair: n = 12, key group of 140, quipus ((1,0,1,1),(1,1,4)) vs ((1,0,2,1),(3,1,1)), reps (2,2,2,3,0,2,3,0,3,0) and (2,2,3,0,3,0,2,2,3,0), |det S| = 40, signature (11,1), same Coxeter polynomial (x+1)^2 (x^2+1)(deg-8 factor). q_S isomorphic (all 9 cross pairs of 3+3 members), so S_a and S_b have the same genus (Z-congruent if Nikulin 1.13.3 applies; UNVERIFIED). Phi-equivariant refinement: the four pieces ker f(Phi) (f = x+1, x^2+1, the degree-8 factor, (x+1)^2) have equal rank, signature, det (24, 196, -735, 24) and isomorphic discriminant forms. Local module test: the commutant lattice Hom(Phi_a, Phi_b) (rank 14) has an element invertible mod p for p = 2 (exhaustive, 65 elements to a hit), 3, 5, 7 (random, first or eighth try), so no obstruction at those primes; a global unimodular element was not found (random combinations with coefficients in {-1,0,1}, 20 000 tries, least |det| = 4). NOT decided: global module isomorphism and the isometry stage.
A refutation of the first part would be an error in the q_S routine; positive controls below.

## Evidence
| key group (n, size) | det S (abs) | cross pairs (3+3 members) separated by q_S |
|---|---|---|
| 11: 242 / 282 / 36 | 20 / 32 / 52 | 9/9 each |
| 12: 372 (second) / 754 / 54 (det 64) | 44 / 28 / 64 | 9/9 each |
| 12: 54 (det 60) | 60 | 9/9 |
| 12: 140 | 40 | 0/9 (open pair) |
| 11: 50 (outside E-259's eight) | 32 | 9/9 (`maverick_cartanZ_n11.out` L49) |
Which groups are profile-equal comes from re-running `maverick_quipuprofile.py` (E-259 list; n = 12 groups 372, 754, 54, 140, 54).
Mechanism (n = 11 size-242 pair, A = Z/20 in both): q on the order-4 element 5/4 vs 1/4 mod 2, on the order-5 elements {4/5, 6/5} vs {2/5, 8/5}. Milgram's formula holds for both (signature 9 = 1 mod 8) since the 2-part and 5-part differences compensate.
Controls: (1) q_S isomorphic after a random unimodular congruence and after a vertex swap (both reps); (2) calibration, n = 9: the F-010 group (size 36, known inequivalent, equal det/signature/SNF probes) splits into 18 + 18 by q_S; (3) verified-move orbits (`lnaMoves.closureUnderMoves`, same derived class with an explicit path): 0 orbits with two q_S classes (see `maverick_orbitcontrol.out`; members evaluated are those whose class name is a key row: 132, 185, 319, 645 members, a small subset of each orbit set at n = 9..12).
Whole key table (multiset of (order, q) values over A, a separator; |det S| <= 4000): key groups split by (det, signature, q_S): n = 9: 1 of 19, n = 10: 3 of 40, n = 11: 4 of 79 (exactly the four groups with >= 2 quipus), n = 12: 15 of 176 (12 of the 13 groups with at least two quipus, `maverick_cartanZ_n12.out` L1-L214 and `maverick_genusscan_n12.out` L1; the 140 group is not split; three further groups of sizes 3837, 95, 21 hold unnamed members and are split). The n = 10 split group sizes are 323, 176, 45; none matched to the E-178 periodic pair by me. Periodic groups (E-197, E-242, E-232) show no q_S split, consistent with their record (E-197 had q_G constant at n = 11).

## Reproduction
```
timeout 10m .venv/bin/python workshop/rounds/097/maverick_cartanZ.py 11 3   # about 2 min each; output maverick_cartanZ_n11.out
timeout 10m .venv/bin/python workshop/rounds/097/maverick_cartanZ.py 12 3   # about 5 min; maverick_cartanZ_n12.out
timeout 10m .venv/bin/python workshop/rounds/097/maverick_genusscan.py N    # N = 9..12, 3 s .. 3 min; maverick_genusscan_nN.out
timeout 10m .venv/bin/python workshop/rounds/097/maverick_control.py        # controls, seconds
timeout 10m .venv/bin/python workshop/rounds/097/maverick_pair140.py        # open pair: Phi, module search, local tests
timeout 10m .venv/bin/python workshop/rounds/097/maverick_pieces.py         # open pair: ker f(Phi) pieces
timeout 10m .venv/bin/python workshop/rounds/097/maverick_orbitcontrol.py N # N = 9..12
```
Needs numpy, sympy, python-flint, cypari (the venv lacked all; I ran `uv pip install` of numpy sympy networkx polars cypari python-flint).

## Prior record
E-197 computed q_G (same genus as q_S) only inside the periodic groups and found it constant; E-209 and E-259 used SNF(xC + C^T), det and signature, which are coarser. E-259 lists the 8 pairs as unseparated; this separates 7. grep of research/ for "discriminant form" / "genus" finds only E-197 and E-209; nothing recorded on q_S across quipu pairs. E-259's 8 pairs not in RETRACTIONS. The Phi-based F-047 profile and q_S are independent: the pairs have equal profile and different q_S.

## Code changed
None in quivermutation/ or tests/. New scripts in workshop/rounds/097/ (maverick_*.py, outputs under 200 KB).

## Next
- skeptic: independent check of q_S for one n = 11 pair (hand: Z/20 pairs above) and of the Nikulin reference; scholar: fetch the uniqueness-in-genus theorem (class = genus) and the Milgram check.
- Whether q_S is constant on derived classes is a theorem (Z-congruence of C) but only the verified-move orbit control tests it; run on larger orbits and on E-117/F-010 merged classes: a q_S split inside a known merge would be a bug.
- Open pair: hermitian-lattice genus over Z[Phi] via PARI (bnfinit on the degree-8 factor), a proposal for OVERNIGHT if wanted; I estimate hours to build, not run.
- Use q_S as a new key component: it splits 4 of 79 groups at n = 11 and 15 of 176 at n = 12, including groups with unnamed members (three at n = 12), where it certifies inequivalence for LNAs with overlapping relations.
