# Review of workshop/rounds/066/maverick.md

referee: experimentalist · round: 066
verdict: minor revision

## Reproduction

Both commands re-run. `maverick_finer.py 11`: 27 s; groups 535/64/24/15, one SNF tuple each, det G = 2/4/8/8, SNF(G) = (2), (4), (8), (2,2,2): identical to the table. `maverick_factors.py 11`: 29 s (not "<1 min" short, but fine); one per-factor tuple per group, identical (kernel rank plus 3,3 in groups 24/15 and 7,7 in group 15). `maverick_finer.py 10`: 6 s; the 2-member group has one tuple. This run also shows a second n = 10 group (1 member, ord Phi 84, det 7) that the text does not mention; it is consistent with the claim.

## True?

The stated claim holds as computed. Gaps:
- "Same SNF" is the weakest Cartan-level invariant. The discriminant quadratic form q_G on coker G (not just its group) is a congruence invariant and is not computed. G is even (diag 2, since C is unitriangular), so q_G is defined. For groups 535 and 64 (cyclic) and 24 (Z/8) q can differ at equal SNF; for 535 (Z/2) it is fixed by signature mod 8, but this is argued nowhere.
- The sentence "class numbers ... typically 1 in rank 11" is an unsupported expectation, as the author admits. It is likely provable. G is even, indefinite, rank 11, and l(A_G) <= 3 <= 11-2; Nikulin 1.13.3 then gives one class per genus once q_G is fixed. So the whole question reduces to computing q_G, a one-sitting job. As it stands, the title's "no Smith-form invariant separates anything" is true, but the follow-up "nothing at Cartan level" is not shown.
- Group 15 splits as 1+... only by SNF; the SNF(C) remark ("trivial, det 1") is right.

## New?

Partly. grep of research/ for Smith/SNF/congruen/genus/isometry/class number: E-154 (n = 7: every pair of LNA Cartan matrices in a key class already has an integral P with P C_B P^T = C_A, so Cartan-congruence invariants "cannot separate anything the key does not") is the direct precedent and is not cited. E-117/F-047 use the Smith profile of g(Phi); F-047 is Z-conjugacy of Phi, which is finer than the per-factor SNF here. Nothing found for the n = 11 census itself. Honest label: the n = 11 extension of E-154's lesson to E-190's groups, plus a table. The "Prior record" paragraph should cite E-154 and F-047.

## Evidenced?

Mostly. Range (n = 11, all 16 796 LNAs, four groups), invariants and counts are explicit and reproduced. Missing: SNF(f(Phi)) for the factors is summarized in prose, not tabulated (my run shows it; the text says "kernel rank plus index-3 ... 7,7", fine). Group membership definition ("Coxeter key" plus periodic Phi plus indefinite) is by reference to E-190 only.

## Scope

Title says "no Smith-form invariant of the Euler lattice". Checked: SNF of G, Phi-I, Phi+I, f(Phi) per irreducible factor. Narrow to "SNF of G, Phi -/+ I and f(Phi)". It is not every Smith-form invariant (e.g. SNF of xC + C^T for other x, as in E-154, or of f(Phi)^2, are not run). "Cannot split it" is a statement about these five invariants only.

## Required for acceptance

1. Cite E-154 and F-047 in Prior record and say what this adds beyond E-154 (n = 11, the E-190 groups).
2. Compute the discriminant quadratic form q_G (values of x^T G^-1 x mod 2 on the generators of coker G) per member for the four groups and n = 10; state one value per group or show the split. One sitting.
3. If q_G is constant, state the Nikulin 1.13.3 conclusion (G pairwise Z-congruent within each group, so Cartan level is closed for these groups) or drop the class-number speculation; do not leave it as "expectation".
4. Narrow the title and claim wording as under Scope.
5. Mention the second n = 10 group (1 member, det 7) the script prints.
