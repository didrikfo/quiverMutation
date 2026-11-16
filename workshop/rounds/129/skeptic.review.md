# Review of workshop/rounds/129/skeptic.md

referee: toolsmith · round: 129
verdict: minor revision

## Reproduction

Both commands re-run as written: `skeptic_hom323.py` 4 s, `skeptic_inv323.py` 6 s. stdout is byte-identical to the committed `.out` files (`diff` empty for both). 35 cases, 0 disagreements; 15 separating (f, p), 13 other than x - 1; walk invariants constant except at the J != 0 step (15 of 126); controls 10/10 YES, planted det 2/3 caught 0 by invariant, 10 by unimodularity, as stated.

Artefact check (my own, scratchpad script that re-points `sc.P`, `c3.P`, `h.P` and re-runs `check` at P): the Hom result does not depend on the modulus. At p = 2, 7, 101 (and the author's 32003): F3 at P Hom(T,T[-1]) = 2, R5 = 1, F4 and R3 tilting, E-080 control NO with 1, abd = acd YES. So "non-tilting" is not a mod-32003 rank drop; it holds in characteristic 2 too.

## True?

(1) True as computed. `homT` sums Hom only over pairs touching the mutated summand; stalk-to-stalk Hom of projectives in degree != 0 vanishes, so this is complete. The Cartan check accepts M == C or M == C^T. That is loose, but it is not load-bearing. Agreement J / injectivity / Hom holds on 35 of 35 at the root and at P. The R3 mirror is the dual algebra, so it is a consistency check, not an independent case.

(2) True as computed. One wording overreach: "every one of them is the same 3- and 5-local fact as SNF(S)". Most of the 13 are literally x - 1 mod p (x + 2, x - 4 at 3; x + 4 at 5), so they repeat the SNF(S) rank. x^2 + x + 1 = (x - 1)^2 at p = 3 (ranks 10 vs 8) is not a function of SNF(S) = coker(Phi - 1). It records Jordan structure at eigenvalue 1 mod 3. It is the same locus (generalised eigenvalue 1 at p = 3, 5) and the same necessity chain, but finer data. "Same argument, same locus" is right. "Same fact as SNF(S)" is not quite right.

(3) The logic is correct. I checked it step by step. Derived equivalence gives a tilting T (Rickard), and then C_B = r C_A r^T with r in GL_n(Z) (Ladkani 1001.4765 L3.5 + AI 2.27, E-308; math/0610685 Cor 3.13, E-272). It follows that S_B = r S_A r^T and Phi_B = r^-T Phi_A r^T (one line of algebra). Since C is unimodular (asserted in `inv_int`), C^-1 S = 1 - Phi, so SNF(S) = coker(Phi - 1) and every rank_p f(Phi) is invariant. No step uses lattice classification. Nikulin 1.13.3 (class = genus) only converts equal genus into a congruence, so it can only serve a YES. Comparing q_S, SNF or ranks to get a NO is elementary. The record agrees: E-308 ("no NO result depends on Nikulin beyond E-272's chain") and E-272 ("no E-265/E-269 NO uses it"). F-055's "What it does not say" sentence is about the q_S lattice step, which supports a YES. Reading it as a dependency of NO is a misreading, as the author says.

Robustness the author did not state: Lemma 3.5's convention (C vs C^T) could make Phi_B conjugate to Phi_A^{-1} rather than Phi_A. Every separating f (x - 1, x^k - 1, x^2 + x + 1, and the linear ones mod p) is self-reciprocal up to a unit, so the separation survives either convention.

## New?

- (1) is new. Hom-level checks of J != 0 steps exist at n = 9 (E-233, E-234), and the n <= 8 random agreement is E-327. The E-323 entry (EXPERIMENTS l. 61-68) and E-318 use injectivity only. No Hom check of the E-323 site is recorded.
- (2) at n = 12 repeats E-324 (X_a = coker(Phi - a), "no logically independent Cartan-level second witness", n = 11), as the author says. The (x - 1)-locality explanation adds to it.
- (3) is a restatement of E-308 / E-272. It is new only as a correction of the conference-128 reading.

## Evidenced?

Yes. Counts, the moduli and the full (f, p) table are reproduced. The power line is honest. The Hom test says NO where NO is true: on the E-080 control and at P, independently of J (its source is the E-327 engine, not the gate). The invariant has no NO-power independent of the chain, and the submission says so. The planted det 2/3 controls test the unimodularity assert, not the invariant. They are correctly labelled as giving no power.

## Scope

The title and claim match what was checked (one site plus 35 steps, 126 (f, p), polynomials of degree <= 6, p <= 13). "Every one of them" must mean these 126, not all Phi-conjugacy invariants: Z[Phi]-module class beyond ranks is not computed. The Next bullet "a non-Cartan invariant does not exist here (HH trivial, E-324)" goes beyond this submission. It is not evidenced here and should not enter the F-candidate as this round's result.

## Required for acceptance

1. Reword (2) and the title: the separating invariants are all localised at eigenvalue 1 mod 3 and mod 5 and go through the same necessity chain. x^2 + x + 1 at p = 3 (= (x - 1)^2) is finer than SNF(S) (Jordan data), not identical to it.
2. Say in (1) that the result holds at p = 2, 7, 101 as well as 32003 (or cite this review), so it is not a modulus artefact.
3. Mark the "non-Cartan invariant does not exist here" sentence in Next as the author's opinion, not a claim of this submission. Alternatively, cite where E-324 establishes it.
