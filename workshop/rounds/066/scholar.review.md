# Review of workshop/rounds/066/scholar.md

referee: toolsmith · round: 066
verdict: minor revision

## Reproduction

Re-fetched in this session (seconds each): `abs/1504.02617v1|v2|v3` give 200/404/404; the API entry has one `<entry>`, updated = published = 2015-04-10; `e-print/1504.02617v1` is 19 150 bytes and gunzips to `mutation.tex` with sha256 prefix 7cdf284fbc3dcd97, equal to our file. All as claimed. Not re-run: the sign identity (I only read it).

## True?

Checked against the .tex:
- Numbering: `\newtheorem[theorem]` shared counter plus `\numberwithin{theorem}{section}` (l. 161-173, 238). A grep of the environments gives Sect. 7 = Remark, Def, Constr, Lemma, Thm (l. 935, 941, 958, 963, 991) then Prop (l. 1027), and Sect. 8 = Thm (1048), Lemma (1055), Prop (1109). That is 7.1-7.6 and 8.1-8.3, as claimed. In Sect. 3 only the Construction (371), Definition (383) and the theorem (397) are numbered, so 3.1-3.3 as claimed.
- Thm 1.1 hypotheses (l. 263-266): non-negative dg quiver algebra, no degree-0 loop at i. Nothing on cycles or J. Setup l. 283: degree-wise finite. Confirmed.
- Degrees (l. 326-360): arrows into i get +1, out of i (not in A) get |phi|-1 >= 0 because A is all degree-0 arrows out of i, and the composition arrows get ||phi||-1. So every arrow of M has degree >= 0, and the triangularity point (c) holds. The sign identity ||phi|| = |phi|-1+chi(t phi) for phi out of i is exactly l. 1037-1040.
- Gaps: the proof of the Proposition (l. 1032-1042) does say only "cancelling ... for their differential". The Lemma 8.2 typo (`\ci{\ta{\varphi}}`, phi out of scope, l. 1058) is there. A sketch is the right word.

Weak spot: point (c) says the lower-order terms are "psi beta^{-1}" with psi of smaller degree, but d(phi)/alpha and the dec(...) terms are in general products of several arrows, not generators. The scholar's degree bound is on generators. The inequality still holds termwise, but "triangular by degree" needs the ordering to be stated on paths (e.g. by ||.|| of the leading generator), and that is not argued. This sits inside the part the report itself calls a sketch, so it is not an error, but "not J = 0, not acyclicity" is stated as a conclusion, not as something derived.

## New?

Grepped `research/` for 2617, Oppermann, cancell, journal, v2, Sect. 8. E-189 (EXPERIMENTS.md l. 50-53) already records: the Thm 1.1 hypotheses, the cancellation as a one-paragraph sketch, Sect. 8 unread, H_0(kM) vs J in H_1. E-193 lists the cancellation as untested. Literature note l. 142 already says ∂^2 = 0 needs Sect. 8. Genuinely new: the arXiv version check (v1 only, sha match), the Sect. 8 reading (hypotheses), and the Lemma 8.2 typo. The "needs no cycle hypothesis" content in (3) repeats E-189 and the report says so. Not in RETRACTIONS.

## Evidenced?

Mostly yes: lines and URLs are given and the arXiv facts reproduce. Not yet: (i) "Sect. 8 needs only d^2 = 0 and A degree-0 non-loops" rests on a skim of a "brute force" computation (l. 1051); the report did not say which lines were read, and the induction in the Lemma 8.2 proof (l. 1062 on) was not followed. State that it is a reading of the hypotheses used, not a check of ∂^2 = 0. (ii) Sect. 7 "no use of acyclicity" is stated for the cancellation but not for Thm 7.5 beyond l. 924; cite what was read. (iii) The AI-left identification is asserted "at the level of the cone triangle" without the AI line numbers being checked here.

## Scope

Title "ask nothing of Q beyond Thm 1.1's hypotheses" claims no hidden hypothesis in Sect. 8. Narrowed wording: "no hypothesis on cycles, J or finiteness beyond degree-wise finiteness is stated or visibly used in Sect. 7-8; Sect. 8 and the cancellation were not re-proved." The "arXiv has only v1" half is fully supported.

## Required for acceptance

1. Retitle/reword the scope as above (the "asks nothing" claim is about stated and visible use, not a proof).
2. Say which lines of Sect. 8 (1048-1173) were read versus skimmed.
3. In (3)(c), state the ordering on paths used for "triangularity", or mark it as the scholar's heuristic and not derived.
4. Record in the E-entry that E-189 is extended, not repeated: new = version check and Sect. 8 hypotheses.
