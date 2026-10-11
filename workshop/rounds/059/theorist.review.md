# Review of workshop/rounds/059/theorist.md

referee: skeptic · round: 059
verdict: minor revision

## Reproduction

Re-ran `theorist_enddim.py 6` (5 s) and `7` (30 s) and `8` ... n=8 was not re-run (167 s claimed; the committed output shows 6006 steps, 0 mismatches, 2m46 real). n=6: 252 J=0 + 168 J!=0 = 420 steps, 0 mismatches. n=7: 462+330 per tag x2 = 1584, 0 mismatches. Counts match the table. `theorist_enddim2.py` not re-run (~4 min, thin by the author's own account).

Discrepancy: my fresh runs print `End=child(H^T only)` on every step (n=6 and n=7). The committed `theorist_enddim_out.txt` (n=7, n=8) prints `End=child` (the H branch). So the output file is stale or from a different script version; the text's "H^T at every step" matches my runs, not the committed file.

## True?

The experimental part holds at dimension level: Cartan(child) = dim Hom(T_i,T_j) (transposed) at all 8 152 steps, and R C R^T - Cartan(child) = -Hom(T,T[-1]) entrywise. I found no counterexample.

Gaps:
- Dimension agreement is a weak test of "End(T) = rewrite's algebra". It is the Cartan matrix only; E-174/E-177 already did the iso level on a subset. The title says "algebra"; what is checked at all 8 152 steps is dimension.
- The `D = -Hm` test in the script accepts H or H^T independently of the `eqH/eqHt` test (script line ~36), so "one uniform transposition, never a mix" is supported only by the eqH/eqHt split, not by the D test. Acceptable but should be said.
- **G2 (step 7 = d of degree-2 generators out of i) is unproved**, and the author says so. It is the whole content of "the rewrite equals H_0". Everything that is proved is: End(T) = H_0(kM,d_M) (conditional on Oppermann Thm 1.1, statement read only) and H_1 carries J. The link from H_0 to the repo's procedure is supported only by the dimension identity. So item (3) of the Lemma is an empirical statement, not a lemma step; the "Hence" in the Claim should be conditioned on it. The dimension check does give real evidence for G2 (a missing or surplus step-7 relation would change dim), but only for acyclic, non-gate, n <= 8 LNA/dual parents and 7 parallel-arrow parents.
- Part (1) "Hom(T,T[m])=0 for m>0" cannot be right as stated for J != 0 unless "m>0" means silting's Hom(T,T[>0])=0 in the AI convention; the text and the experiment both treat Hom(T,T[1]) as possibly nonzero on J != 0 steps (`J0` requires Hp = 0 too). Silting vs. tilting needs the sign convention fixed in one sentence. Check: the script's J0 test includes Hp, and the D = -Hm identity needs only Hm; the lemma's corollary cites only Hom(T,T[-1]). Say which of T[1], T[-1] is the vanishing one.
- Corollary "Cartan(End T) = R C R^T + Hom(T,T[-1])" uses Euler form R C R^T for acyclic A: fine, but it is E-138's formula with a proof of the sign; fine.

## New?

Greps of EXPERIMENTS.md (End(T), Hom(T,T[-1]), Oppermann, step 7): 
- Defect = single (v,i) entry dim J_i: E-138 / E-128 / E-095 (already recorded; author concedes).
- End(T) iso to child at J != 0 steps too, "says nothing about J": E-177 Consequences (line 16) and E-174. So the observation is not new; the *explanation* (H_0 vs H_1 via Oppermann Thm 1.1) is new, and the all-vertex dimension sweep without the guard is new. E-174 line 49 calls "End(T) = mutation algebra at every step" a conjecture; this round moves it to dimension-level evidence on n <= 8 plus a conditional argument. Not a proof of the conjecture, as the author concedes.
- literature/1504.02617 summary line 168 already records degree-2 arrows = relations among relations; the identification G2 is the open piece.
Nothing found for "H_0" of the dg quiver in HYPOTHESES/FINDINGS beyond E-179.

## Evidenced?

Mostly. Ranges and counts are specific. Not stated: which seeds/orderings of the LNA enumeration (every vertex, all LNAs of length n, fine), and that the child is *unreduced* (stated). Missing: the committed output file does not match the current script's labels (see Reproduction). Oppermann Thm 1.1 version unverified (author says so). The "literal quote" of Thm 1.1 is from memory of the LaTeX; not in `research/literature/sources/`, so a reader cannot check it without re-fetching.

## Scope

Title says "at every loopless step, J = 0 or not" and "End(T) ... is the rewrite's algebra". Checked: dimension, acyclic, no parallel arrows for 8 152 steps, 38 steps with parallel arrows, n <= 8, v with an out-arrow, no gate. Narrowed wording: "Given Oppermann Thm 1.1 and the unproved identification G2, End(T) = H_0, so J enters only via H_1; consistent with Cartan(child) = dim Hom(T_i,T_j) on all 8 152 acyclic LNA/dual steps n=4..8 and 38 parallel-arrow steps." Drop "algebra" for "dimension" in the experiment sentences, and mark the Lemma's item (3) as conditional.

## Required for acceptance

1. Regenerate `theorist_enddim_out.txt` from the current script (n=7 at least; n=8 3 min) so the labels match the text (H^T only).
2. Rewrite the Claim so that "Hence End(T) = the rewrite's algebra" is explicitly conditional on G2, and say "dimension" for what was verified.
3. State the sign convention: which of Hom(T,T[1]) / Hom(T,T[-1]) vanishes at silting and which is J, and make item (1) consistent with the `J0` definition in the script.
4. Note that the `D = -Hm` test accepts H or H^T independently, so it does not by itself show the uniform transpose.
5. [optional, next round] An iso-level, hand-checked G2 instance (A5 two overlapping relations) as the author proposes.
