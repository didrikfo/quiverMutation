# Theorist notebook (rewritten round 083)

## What I believe now
- 083 (T5): Prop. 7.6 cancellation (kM -> N quasi-iso) is PROVED by hand under: non-negative grading + no loop at i of ANY degree. Not needed: D(i) acyclic, finite dimension, J. Proof: omega = ||w|| + #stars filtration makes gr d = d_1 (d_1 psi = sum a* (a psi), d_1(phi a^-1) = +-phi a*); second filtration by W-count reduces to Anick/Morse Lemma A (monomials x a*, U and Y disjoint). Weakest step: acyclicity of the Morse matching (cited). Machine: Part A shape of d on 4 249 letters, Part B toy exact ranks to weight 12 (workshop/rounds/083/theorist_starfilt.py).
- Every `noCycles` parent in search.py satisfies the hypothesis: Q^(inf) of an acyclic parent has all arrows s->t with s<t, so no loops. E-211 also satisfies it: cancellation fine, step 7 of the repo is the failure.
- OPEN: loops of degree >= 1 at i (Step 4): d_1 is not monomial (6968 of 10206 loop letters differ). Not reached by the search.
- Older: G2/step 7 (065, 067): repo kernel = Oppermann's iff tails vertex-simple, e.g. D(i) acyclic (E-199). 073/078: guard does not force D(i) acyclic; E-211 counterexample; search only mutates acyclic parents. 075: Gamma_red = N by hand (E-217). T7: single-mutation LNA rule proved for right (E-187) and left (E-205); H-010 supported, unproved; 081 two-bystander data (16/16 admissible cells stay at overlap 2, controls reach 0).

## What I tried
- 083: read Constr/Def 3.2 and Sect. 7 LaTeX, found d = d_0 + d_1 by star count, wrote theorist_starfilt.py reusing 069's executable Def. 3.2.
- 081: theorist_twoby.py; 078: cycD; 075: Gamma_red; 073; 070; 067; 065; 062.

## Next
- Step 4 loops at i: try loop-count filtration / matrix form (psi, theta1, c, theta2).
- A real-parent machine check of kM -> N (needs products with signs); scholar to vet Lemma A and the journal version.
- Remaining links to "repo child = End(T)": d^2 = 0 beyond numerics, W1/W2, Thm 7.5.
- T7: two pair-sharing bystanders on opposite sides, m = 4; k = 4 needs OVERNIGHT job; mechanism for a run of three.
- Blind spots: proofs are relative to the printed presentation; a toy check is not a check of the real kM; capped enumerations are not laws.
