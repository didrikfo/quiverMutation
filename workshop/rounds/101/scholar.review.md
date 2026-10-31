# Review of workshop/rounds/101/scholar.md

referee: skeptic · round: 101
verdict: minor revision

## Reproduction

Re-fetched math/0610685, math/0702060, 1404.5460 (LaTeX, own scratchpad dirs, ~1 min). Matched:
- math/0610685: the corollary at l. 723 ("If X ~ Y then 1_X, 1_Y congruent over Z"). Counting the shared `theorem` counter in section 3 gives Prop 3.11 (l. 699), Defn 3.12, Cor 3.13 (l. 723), Lemma 3.14 (l. 738), Cor 3.15 (l. 750). Numbers confirmed (not just inferred).
- math/0702060 l. 1194-1198: the sentence quoted, with P^t C_L P = C~_L. Matches.
- 1404.5460 l. 336: "Corollary 1.13.3 of nikulin ... even indefinite lattice whose rank exceeds (by at least 2) the minimal number of generators of d(Gamma) is determined by its rank, signature and discriminant quadratic form". Matches, label is Corollary.

## True?

Holds. Step 3 checked algebraically: with Phi = -C^{-T}C and C_A = P^T C_B P, Phi_A = -P^{-1}C_B^{-T}P^{-T}P^T C_B P = P^{-1}Phi_B P. The transpose convention is consistent. The form-level argument (non-degenerate form, F isometry => F Phi_A = Phi_B F) is also correct. No counterexample.

Two corrections:
1. "Equivariance not printed" is slightly understated. Ladkani math/0610685's proof of Lemma 3.14 (l. 738-748) writes M2 M2^{-t} = P M1 M1^{-t} P^{-1} with the SAME P that gives the congruence. The statement is only "conjugate", but the proof carries the intertwining P. So the equivariance is a one-line corollary of the printed computation, not an independent proof.
2. Cor 3.13 is stated for 1_X^{-1} (Euler form on simples), not the Cartan matrix; the author's "up to P <-> P^{-T}" remark covers this. State it in the H-024 rewording.
Caveat: the finite-global-dimension hypothesis is required (the author says so); E-265 already notes it is not shown for these LNAs.

## New?

Mostly already recorded: research/literature/math-0610685-sheaves-over-finite-posets.md has "Proposition 3.11 and Corollary 3.13" and "Lemma 3.14 and Corollary 3.15" sections (grep hit). Author acknowledges this. New content: math/0702060 (general algebras, not posets), the Cor-vs-Prop label for Nikulin via 1404.5460, and the H-024 "no source matched" correction. H-024, E-265, E-267 (names the label as "Prop. 1.13.3"), E-269 are the entries to change.

## Evidenced?

Yes for the literature items; each has a file, line, and quoted text. `power: n/a` is fine for a literature round. The Nikulin item is honestly marked secondary-only: 1404.5460 is a citing paper, so Nikulin's own text is unread, and E-265's hypotheses (indefinite, rank >= l+2) were not rechecked here. The author says this.

## Scope

Title says "published for Cartan congruence (Ladkani)". Fair for general finite-dimensional algebras of finite global dimension via math/0702060 and for posets via 0610685. Not Nakayama-specific, as stated. "Phi-equivariance ... two-line corollary" is fine, but add that Lemma 3.14's proof already contains it.

## Required for acceptance

1. Amend the "not printed" wording for equivariance to point at the proof of Lemma 3.14 (l. 738-748, same P).
2. In the proposed H-024 wording, note Cor 3.13 is for 1_X^{-1} and give the P <-> P^{-T} translation to C.
3. Keep E-265's Nikulin citation UNVERIFIED at source (Cor. 1.13.3, via 1404.5460 only); do not upgrade.
