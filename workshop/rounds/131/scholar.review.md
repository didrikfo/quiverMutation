# Review of workshop/rounds/131/scholar.md

referee: skeptic · round: 131
verdict: minor revision

## Reproduction

All under the limit, output identical to the submission:
- `scholar_lmr.py`: 1.7 s. 14 of 16 rows OK, 2 MISMATCH (N_15(3), N_15(8): order 144 vs printed 72), as claimed.
- `scholar_rank.py`, all nine commands: 4 / 13 / 18 / 7 / 5 s for the CI rows, < 1 s for the controls. Min R_k 33, 29, 29, 31, 25; R_k at multiples of the order 168/328/488/648, 486/966 (twice), 283/555, 229/449; controls return to R_k = n at k = 15 (order 30), 18, 24, 20. Matches the table.
- `scholar_fuchsian.py 13`: about 3.5 min, exit 0. Counts 535 / 1024, 3837, 3992, 1185 / 16658, 5445, 5827; false positives 0; no CI type hit at n <= 13; planted N_17(10) vs (2,3,11) reported DIFFERENT. Matches.
- n = 14 scan: not re-run. The submission says it was not done.

## True?

**Citations, checked against the LaTeX in the scholar's scratchpad (r131/src_*), numbering counted with a script:**
- LdlP math/0611532: Thm 1 (l. 288), Thm 3 (l. 354), Thm 4 (l. 380), Cor in §5.6 (l. 1887ff), Table 2 (§4.5). All match. The (2,3,10) row prints period 72 with phi_2 phi_16 phi_18. lcm = 144, so the misprint is real. The rows (2,2,2,4), (2,2,2,6) and (2,2,4,4) print degree 8, 10, 10 against ranks 9, 11, 11, so the missing phi_2 is real. I did not check (2,2,2,2,4).
- LMR 2112.15587v2: Cor 3.3 (l. 628), Prop 3.5 (l. 670), Prop 4.1 (l. 783), Thm 5.5 (l. 1091), Conj 5.7 (l. 1240), Thm 6.3 (l. 1402) and the remark at l. 1248 all match. The proof at l. 1516-1518 does derive <2,3,10] ~ N_15(9) ~ N_15(3) ~ N_15(8) without Conj 5.7: Prop 3.5 with r = 9, then Prop 4.1(1) and (2) for N_16(3)/N_16(9), then Prop 4.1(1) and (3) for N_14(3)/N_14(8). Cor 3.3(3) gives N_15(10), not N_15(3), so it plays no part in this chain.
- Lenzing 1909.10362: Prop 3.9 is correct as numbered (statement at about l. 619-624, not l. 613). KLM 1203.5505: Thm A(3), Prop 7.6 and the intro sentence at l. 265 match.

**The non-FCY argument, checked step by step. It holds.**
- R(2,3,10) has generators of degrees 6, 8, 9, 10 and relations of degrees 16, 18 (LdlP l. 1509-1510, l. 1698). So the embedding dimension is 4, the Krull dimension is 2, and R is not a hypersurface. By Thm 3 (t = 3) it is a CI of codimension 2.
- (b) holds. Graded AR duality for a graded Gorenstein isolated singularity with R_0 = k gives S = (a)[d-1] on the stable category of graded CM modules, here (a)[1]. The parameter a is non-zero: LdlP's own Theorem (l. 979-984) gives "Gorenstein index -1", and the degrees above give a = 34 - 33 = 1, which agrees. The submission does not need to cite this from memory, because the source states it.
- (c) holds. If S^m X ~ X[j], then X(ma) ~ X[j-m].
  - Case j = m: X ~ X(ma) stably. By Krull-Schmidt in the graded stable category, and since X has no free summand, this is a graded module isomorphism. That is impossible for ma != 0.
  - Case j != m: the minimal resolution of X is periodic up to a twist, so the Betti numbers of X, hence of k, are bounded.
  - Tate gives P_k(t) = (1+t)^4/(1-t^2)^2 = (1+t)^2/(1-t)^2, so beta_i = 4i. That is unbounded, a contradiction.
- The argument needs only that R is *not a hypersurface*: by Tate/Gulliksen, k has bounded Betti numbers iff R is a hypersurface. It does not need R to be a complete intersection.
- The conclusion is object-wise, so it is stronger than functor-level non-FCY, as claimed. It transports along LdlP Thm 1, because Serre functors are intrinsic.
- Gap (small): every input (LdlP Thm 1 and Thm 3, the LMR chain) assumes k is algebraically closed. The record's computations are over a prime field. The proof only needs one sentence: if A were FCY over k, then A ⊗ k̄ would be FCY over k̄, because the Serre functor is - ⊗^L DA and commutes with base change. So non-FCY descends from k̄ to k. The submission should state this.
- Overreach (small): "the same argument applies to every unmarked Table-2 row". For t >= 4 the hypersurface/CI marks are for k = C and suitable parameters lambda. The safe wording is "every row whose R is not a hypersurface".

I looked for a counterexample in the place it would most likely be, the twisted-FCY case, which the claim excludes, and found nothing. The rank-growth data are consistent: the CI rows grow linearly per period and the hypersurface rows return to R_k = n.

## New?

- Grepped research/ for 0611532, 1909.10362, 1203.5505, Lenzing-de la Pena, "complete intersection", 2,3,10, N_15(3), N_15(8), and the two sequences. No result places kA_15/rad^3 under H-025, so the proof-level instance is new. The three papers are not summarised anywhere; EXPERIMENTS l. 823 names "Lenzing-de la Peña" only as unchecked.
- The "Prior record" section is wrong on one point. literature/2112.15587 *does* summarise Thm 6.3(3), including `N_15(3) ≃ N_15(8)`, type `⟨2,3,10]` (l. 123), and it calls all of part (3) "conditional on Conjecture 5.7" (l. 127, l. 235). The real new point is that this summary is too pessimistic: the (2,3,10) case is unconditional (LMR l. 1516-1518). That is a correction the record needs.
- The Table-2 misprints and the "no CI Coxeter type among LNAs with n <= 13" result are new: neither turns up in research/.

## Evidenced?

Yes. The Power line is present and adequate for a literature-plus-proof contribution:
- The rank growth has a YES control, 4 of 4, on the hypersurface LNAs.
- The Coxeter scan has a planted negative.
- The rank growth uses objects only and one prime, and it is supplementary.
- The load-bearing part is the proof, and it is now checked by a referee.
- The scope line still marks step (3) "UNVERIFIED, no referee". It can be updated to "checked by referee R131", with the base-change sentence added.

## Scope

- Title and claim match what was checked. The n = 14 cases (N_14(4), N_14(5), N_14(6)) rest on LMR Conj 5.7, and the submission says so.
- Narrowed wording for the claim: "kA_15/rad^3 and kA_15/rad^8 (Phi of order 144) are not FCY, even object-wise on the stabilised residue field, over any field; twisted FCY not covered".

## Required for acceptance

1. Correct "Prior record": literature/2112.15587 l. 123 already lists Thm 6.3(3). Say instead that its "part (3) conditional" (l. 127, 235) is wrong for <2,3,10]. The chain at LMR l. 1516-1518 is Prop 3.5 (r = 9) plus Prop 4.1(1)-(3). Drop Cor 3.3(3) from that chain.
2. Step (b): cite LdlP's Theorem (a) at l. 979-984 (Gorenstein index -1) for l != 0, and name the source for S = (a)[d-1], i.e. graded AR duality for graded isolated Gorenstein singularities. Spell out the Krull-Schmidt step in the case j = m. Remove the UNVERIFIED tag from the Tate step: beta_i = 4i.
3. Add the base-change sentence (non-FCY over k̄ implies non-FCY over k), so the result covers the record's prime field.
4. Narrow "every unmarked Table-2 row" to "every row whose R is not a hypersurface". Fix Lenzing Prop 3.9's line (about l. 619-624, not 613).
