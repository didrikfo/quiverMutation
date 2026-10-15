# Oppermann's Sect. 7-8 state and visibly use no hypothesis on cycles, J or finiteness beyond Thm 1.1's (not re-proved); the cancellation is a sketch; arXiv has only v1

author: scholar · round: 066 · kind: result
thread: T5 · bears on: H-015 (open), E-189, E-193
scope: reading (not re-proof) of `research/literature/sources/1504.02617v1-oppermann-quivers-silting-mutation.tex` (sha256 prefix 7cdf284fbc3dcd97), lines 252-400 and 895-1173 (Sect. 1, 2, 3, 7, 8). Sect. 8 (l. 1048-1173) read in full, proofs of Lemma 8.2 and Prop. 8.3 followed for structure, signs not re-derived line by line (see Response). Claim is about hypotheses stated or visibly used, not about truth of Sect. 8 or the cancellation. Not covered: Sect. 4-6 (examples, Ginzburg), a re-proof of the cancellation, the journal version, the left/right convention at module level. No script, no code.

## Response to referee

Referee verdict: minor revision. Items 1-4 done; edits below.

1. Scope narrowed (title, scope line, Claim (3)): "no hypothesis on cycles, J or finiteness beyond degree-wise finiteness is stated or visibly used in Sect. 7-8; Sect. 8 and the cancellation were not re-proved".
2. Sect. 8 lines. Read in full this round: l. 1048-1144 (Thm 8.1, Lemma 8.2 and its proof l. 1061-1107, Prop. 8.3 and its proof l. 1113-1144); l. 1146-1173 is the bibliography, skimmed. In the Lemma 8.2 proof I followed the induction and the case split (st p_2 = ta p_1 = i or not) and the inductive substitutions; I did not re-derive every sign. Observed: l. 1105 ends with `(-1)^{||p_2||} dec(dp/alpha)` (unbalanced bracket) while l. 1104 has already absorbed that sign, and the statement (l. 1058) has `-(-1)^{chi(t p)}`; the signs of the last step do not visibly agree (typographical or real, undecided). Earlier skim claim "needs only d^2 = 0 and A non-loops" is therefore a reading of the hypotheses used, not a check of the identity. One further observation: l. 1124 and l. 1140-1141 use that d(red d phi) is a combination of paths starting with arrows of A; since d(x alpha) = dx alpha +- x d(alpha), this needs d(alpha) = 0 or d(alpha) in the span of paths starting in A. I find no such statement for alpha in A in Sect. 2-3 or 8 (Sect. 7 uses d alpha_* = -bar alpha only). Either implicit or a gap; not tested (a d^2 = 0 test on a small (kQ, d) with d(alpha) != 0 would settle it; not done this sitting, which was a reply to the review).
3. Sect. 7 (item ii of "Evidenced"): read l. 895-1042 (Sect. 7 up to and including Prop. 7.6). Thm 7.5 (l. 991-1004): the quotient Gamma/Gamma_red is argued acyclic from d alpha_* = -bar alpha (l. 1004, l. 924); no acyclicity or J hypothesis on Q appears in l. 895-1042. The AI line numbers (item iii) were not checked this round; "AI left mutation" stays as stated in (4), cone-triangle level only, source AI l. 928-940 unchecked here.
4. Ordering for triangularity in (3)(c): marked heuristic, not derived. What I did check: the term `partial(phi) alpha^{-1}` has paths psi_k ... psi_1 alpha^{-1} of total degree ||phi||, and all arrows of M have degree >= 0, so the generator psi_1 alpha^{-1} has degree ||psi_1|| + 1 <= ||phi||, strictly below ||phi alpha^{-1}|| = ||phi|| + 1. Not checked: the Delta-expanded terms of `dec(d phi / alpha)`: there a generator psi beta^{-1} can have degree up to ||phi|| + 1 (if psi carries all the degree), equal to ||phi alpha^{-1}||, so the strict decrease is not shown by degrees alone and a finer order (on paths, e.g. by length or by leading generator) is needed. So "triangular by degree" is the scholar's heuristic and the Prop. is still a sketch.
5. E-entry: E-189 is extended, not repeated. Old (E-189): Thm 1.1 hypotheses, cancellation as sketch, Sect. 8 unread. New here: arXiv version check (v1 only, sha match), Sect. 8 read for hypotheses, Lemma 8.2 typo and sign mismatch, the d(alpha) observation, the limit of the degree argument in item 4.
Not done: re-proof of the cancellation; d^2 test; journal version.

## Claim

(1) Only one version exists on arXiv: v1 (10 Apr 2015). `abs/...v2` and `...v3` return 404, and the API entry lists v1 only (updated = published = 2015-04-10). The e-print of v1 downloaded today has the same sha256 prefix (7cdf284fbc3dcd97) as our file, so our source is the arXiv text. No later arXiv version can change Def. 3.2 or the numbering. Whether the journal version differs is not checkable here (UNVERIFIED: I believe it appeared in Adv. Math. 2017; not fetched).

(2) Numbering matches the source (`\newtheorem` with shared counter, `\numberwithin{theorem}{section}`, l. 161-173, 238; items counted by hand within each section): Thm 1.1 (l. 263, main result); Constr. 3.1 (l. 371, red/dec); Def. 3.2 (l. 383, the differential `∂`); Thm 3.3 (l. 397, "∂^2 = 0", points to Sect. 8). Sect. 7 (counted by hand, one counter per section, subsections do not reset it): Remark 7.1 (l. 935), Def. 7.2 (l. 941), Constr. 7.3 (hat, l. 958), Lemma 7.4 (l. 963), Thm 7.5 (l. 991, Gamma_red is a dg subalgebra and quasi-isomorphic), Prop. 7.6 (l. 1027, cancelling anti-compositions). Sect. 8: Thm 8.1 (`thm.is_diff`, l. 1048), Lemma 8.2 (l. 1055), Prop. 8.3 (l. 1109, "∂ is a differential"). The Sect. 7 and 8 numbers rest on my count; cite them with the line.

(3) What each needs of the input:
- Thm 1.1 statement (l. 263-269): a non-negatively graded dg quiver algebra (kQ, d), i a vertex with no loop of degree 0. Nothing about cycles, parallel arrows, J, or finiteness of Q beyond "degree-wise finite" (l. 283). Longer degree-0 cycles through i are allowed. (Same as E-189.)
- Sect. 7 up to the cancellation: the no-degree-0-loop hypothesis is used so that the arrows alpha in A (degree 0, out of i) have tail(alpha) != i, giving the split maps alpha^*, alpha_* (l. 905-910); no use of acyclicity. Sect. 7's Gamma_red theorem (7.5) needs d(alpha_*) = -bar(alpha) (l. 924) only.
- Cancellation (Prop. 7.6, l. 1027, proof l. 1032-1042, ten lines): uses (a) the formula of Def. 3.2 for `∂(phi alpha^{-1})` and `∂(alpha phi beta^{-1})`; (b) the sign identity ||phi|| = |phi| - 1 + chi_i(t(phi)) for phi out of i, which I re-derived: it turns `(-1)^{||phi||} phi alpha^* = (-1)^{chi(t phi)} dec(d phi/alpha)` into `phi alpha^* = -(-1)^{|phi|} (d phi)/alpha`, as printed, and it agrees with the `hat`-relations at l. 1015-1020; (c) "up to terms psi beta^{-1} with psi of strictly smaller degree", which holds for the term `partial(phi) alpha^{-1}` only (the Delta-expanded terms are not covered, see Response 4; the ordering is a heuristic, not derived) because every arrow of M has ||.|| >= 0 (arrows out of i that survive have |psi| >= 1) so ||psi|| <= ||phi|| - 1 and ||psi beta^{-1}|| = ||psi|| + 1 <= ||phi|| < ||phi alpha^{-1}||. So the argument as written is a triangularity by degree (heuristic reading): it appears to use non-negative grading and degree-wise finiteness, and I see no use of J = 0, acyclicity or looplessness elsewhere; this is a reading, not derived. Step 4 (l. 354-360) needs the loop phi at i to have degree >= 1 (guaranteed by the hypothesis).
- Sect. 8 (∂^2 = 0): the hypotheses visibly used are d^2 = 0 on kQ (l. 1140) and A = degree-0 arrows out of i, none a loop; hand computation ("brute force", l. 1053). No cycle, J or finiteness hypothesis appears; but see Response 2 on d(alpha).

(4) Gaps in the text, none of which I could turn into a counterexample:
- The Proposition's proof states only that one "cancels the phi alpha^{-1} for their differential". The generator phi alpha^{-1} is cancelled against `∂(phi alpha^{-1})`, which is a combination of paths of length >= 2 (it contains phi alpha^*), not against a generator; so this is a quotient of a free algebra by the dg ideal (phi alpha^{-1}, ∂(phi alpha^{-1})), and "the quotient map is a quasi-isomorphism" is not argued (no filtration or acyclicity statement for the ideal). I regard it as a sketch. I did not prove it.
- The Proposition does not state the differential on the quotient (it is `∂` with `dec` replaced by 1, since the arrows containing alpha^{-1} are killed); the reader supplies it.
- Sect. 8 Lemma 8.2 (l. 1055-1058): the second displayed formula has `chi(t(phi))` with phi not in scope (should be p), and its proof ends with the sign `(-1)^{||p_2||} dec((dp)/alpha)` where the statement has `-(-1)^{chi(t p)}`; I did not check whether the two agree. The `∂^2` proof for `phi alpha^{-1}` states `∂(∂(phi) alpha^{-1}) = ...` with "similar to" (l. 1132). Typographical status unknown.
- AI left mutation: Oppermann's P_i^* = Cone(P_i -> sum P_tail(alpha)) (l. 899) is the AI left mutation mu^+(M; D) = add(D, cones of left D-approximations M -> D), source of AI l. 928-940, at the level of the cone triangle. Which of e_i Lambda / Lambda e_i is P_i, and so left vs right for modules, is not settled by this (the skeptic's item).

What this does not say: that the repo's rewrite equals H_0(kM) (G2, E-193, acyclic only); that J = 0 (H_1(kM) is untouched by Thm 1.1, Sect. 7, 8); that the Proposition is true (sketch).

## Evidence

Reading and matching only. Fetched from arXiv: `https://export.arxiv.org/api/query?id_list=1504.02617` (version list), `https://arxiv.org/abs/1504.02617v1|v2|v3` (HTTP 200/404/404), `https://arxiv.org/e-print/1504.02617v1` (19 150 bytes, gz of `mutation.tex`, 78 817 bytes, sha256 prefix equal to our file). Standing arXiv permission, STEERING. Purpose: settle the version question the agenda asked for.

## Reproduction

```
curl -s "https://export.arxiv.org/api/query?id_list=1504.02617"   # one <entry>, v1 only
curl -s -o /dev/null -w "%{http_code}\n" https://arxiv.org/abs/1504.02617v2   # 404
sed -n 1027,1042p research/literature/sources/1504.02617v1-oppermann-quivers-silting-mutation.tex
```
Seconds each.

## Prior record

E-189 is extended, not repeated: it read Thm 1.1 and Sect. 7 and left Sect. 8 and the cancellation open; this reads both and the version question, and adds: no hypothesis hides in them; the cancellation is the only step the theorem rests on that is not proved in the text. E-193 lists "the cancellation step of Thm 1.1 (E-189)" and "the Def. 3.2 reading" as untested: the reading is now checked on signs (b) above, not on the repo. `research/literature/1504.02617-quivers-for-silting-mutation.md` l. 142 says ∂^2 = 0 needs Sect. 8: correct. Nothing in RETRACTIONS touches this.

## Code changed

None.

## Next

- Chair: the sentence for the record is "arXiv has v1 only; the journal version is not read". If the human can supply the journal PDF, the cancellation may be written out there (UNVERIFIED).
- Theorist (optional, one round): either prove the Proposition's quasi-isomorphism (a degree filtration on kM; the ideal is generated by x and ∂x with ∂x = leading term of degree ||x||, lower terms) or state "Thm 1.1 depends on a sketched step" in the E-entry using the premise.
- Skeptic: left/right at module level on the APR example (l. 403), unchanged.
