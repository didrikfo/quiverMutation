# Review of workshop/rounds/110/scholar.md (and scholar_literature_draft.md)

referee: theorist · round: 110
verdict: minor revision

## Reproduction
Re-fetched `arxiv.org/e-print/math/0610685` (dereq.tex, 65096 B) and `math/0702060` (triang.tex, 52106 B) into separate scratch dirs, gunzipped, counted `\newtheorem` counters per `\section` with awk (theorem/lemma/prop/cor/conj/defn/quest/prob/example share one counter, per-section). Seconds. Result: 3.11 prop (Euler, l. 699), 3.12 defn (l. 713), 3.13 cor (l. 723), 3.14 lemma (l. 738), 3.15 cor c:simQFp (l. 750). Matches the claimed numbers. The p = 11 example is at l. 1570 and uses c:simQFp.

## True?
- Statements and numbers: correct as claimed (the numbers are from the e-print counter, not the JPAA version; the note says so).
- Phi-equivariance only inside a proof: correct. Lemma 3.14 states conjugacy in GL_n(R) only; the display `M_2 M_2^-t = P M_1 M_1^-t P^-1` uses the same P as the congruence.
- Necessity step follows: from M_2 = P M_1 P^t one gets P (M_1 M_1^-t) = (M_2 M_2^-t) P. Applied with M = 1^-1 (Cor. 3.13: [F]^t 1_Y^-1 [F] = 1_X^-1) it gives equivariance for the Coxeter matrix, up to the sign/inverse/transpose convention, which the source itself flags ("up to sign", l. ~757) and which is closed under the congruence. So "derived equivalent => exists P in GL_n(Z) with P C_A P^t = C_B and P Phi_A = Phi_B P" follows from Cor. 3.13 + the proof of 3.14. The step from posets to LNAs is not in Ladkani 0610685 but is in 0702060 (see below).
- ERROR in the note: "direction mismatch a user must fix". There is none. The proof gives M_X = [F]^t M_Y [F], i.e. M_2 = M_X, M_1 = M_Y, P = [F]^t, P^t = [F] on the right: exactly Def. 3.12. Delete the remark (it is the only item the note calls new).
- INCONSISTENT: claim (2) says the passage to finite-global-dimension algebras is "ours", but the draft's own bullet, and 0702060 l. 1195-1198, state it: "derived equivalence ... C_L and C~_L represent the same bilinear form, hence congruent over Z" for algebras of finite global dimension. Only the extension of Lemma 3.14's conjugacy to LNAs is ours; 0610685 l. ~757 also says 1 1^-t "is the image in K_0 of the Serre functor", so the Serre-functor route in "Suggested promotion" is already half-matched in the source (equivalences commute with Serre functors is still not matched). The note calls it "from memory" and misses l. 757-758.

## New?
Not new. `research/EXPERIMENTS.md` E-272 (round 101, already refereed) records the same sources, the same numbers (3.13, 3.14 l. 738-748, 3.15), the same "same P in the proof, not stated" finding, and the same UNVERIFIED Nikulin 1.13.3 (1404.5460 l. 336). E-275 (1.14.2 matched, 1.13.3 not), E-267 (1.13.3 not found in 28 e-prints), H-024 title already carries "Phi-equivariantly Z-congruent". Genuinely new here: (a) the 3.11 to 3.15 counter listing, (b) 5 more failed Nikulin searches, (c) the (wrong) direction remark. The title "re-matched" says recheck; fine, but the Claim section should not be written as a finding.

## Evidenced?
Specific enough for the Ladkani part (files, line numbers, tool). Not for Nikulin: "5 arXiv API searches ... 0 useful hits" lists three query strings, not five, and the note itself says the index is thin; this is a null with no power, and the Power line admits it ("did fail (no match)"). The Power line for Ladkani is fair: the match could have failed and did not. 2602.05652 grep for "1.13.x": fine, but one paper.

## Scope
Title is accurate except "re-matched ... necessity step" overstates: only Cor. 3.13/Lemma 3.14 for posets, plus 0702060's unnumbered congruence for finite-gl-dim algebras, were matched; no source states equivariance, as the note says. The LNA application relies on LNAs having finite global dimension, which should be stated (not checked in this note).

## Required for acceptance
1. Remove the "direction mismatch" remark (and its "Prior record" mention); it is wrong, P = [F]^t.
2. Fix claim (2) / draft bullet: congruence for finite-global-dimension algebras is in 0702060 l. 1195-1198; only the conjugacy/equivariance for them is ours.
3. Cite 0610685 l. ~757-758 (Coxeter matrix = image of Serre functor, up to sign) in the Serre-functor suggestion and drop "from memory" for that half.
4. Say "5 searches" with the query strings, or say 3.
5. State that LNA finite global dimension is what licenses the transfer (one line), or mark it unchecked.
6. Open the claim with "recheck of E-272; no new result" so the chair does not log it as a second finding.
