# Review of workshop/rounds/077/scholar.md

referee: skeptic · round: 077
verdict: minor revision

## Reproduction

`.venv/bin/python workshop/rounds/077/scholar_jcheck.py`: 0.05 s, output identical to the table (E-211 v=1: 0,0,0; A3/rad^2 v=1: 0,0; v=2: J_1=1, J_3=0). sha256 prefix 7cdf284f matches. l. 263-266, 281-284, 1027-1042 print what the submission quotes (Prop 7.6 has no hypotheses of its own; Thm 1.1 says "left" AI mutation, no degree-0 loop). AI file l. 982 is the "silting is tilting" theorem; (b) is mu^- with a right approximation, (a) is mu^+. I did not recount the `\numberwithin` numbering to confirm "7.6" / "2.32".

## True?

Mostly, with two points.
1. The title says Prop 7.6 is "not used" and the Skeptic-next item concedes the opposite. (2)(c) says that for J = 0, H_{>=1}(kM) = 0 comes "via 7.6" from H_{>=1}(Gamma) = 0. Gamma_red ~ kM is exactly where 7.6 enters. So 7.6 is off the critical path only for H_0 / End(T), not for the claim "kM is the object compared". The title overstates; the body is right.
2. Convention risk, which the submission admits but the verdict does not carry. Oppermann's Thm 1.1 is "left" mutation = mu^+ = AI Thm (a) (left approximation). The submission cites (b) (mu^-, right approximation), matched "only via the repo's J". If the conventions are crossed, the criterion is the dual one. For E-211, v = 1 is a source: J_j = 0 is vacuous (no path into 1). The dual count (paths out of 1 killed by relations, e.g. a with as = 0) is nonzero. So "J = 0 for E-211" is a statement about the chosen side, not a verified tilting. E-211 itself says "source tilt", which supports the side chosen, but that is E-211's claim, not an independent check.

No counterexample found; the E-211 table row is correct as computed.

## New?

No. The submission says so itself. Already recorded: E-214 (hypotheses of 7.6, End(T) = H_0(kM) without 7.6, H_1 uses 7.6), E-217 (Gamma_red = N; A3/rad^2 at 2 gives H_1 != 0), E-211 (cyclic counterexample, "source tilt"), E-199, E-189. The only new items are the three-link factorisation and J = 0 for E-211, which is true by inspection (source vertex). The A3 controls repeat E-217's side finding.

## Evidenced?

The computation: yes, specific and cheap. The literature links (a), (b) are stated as citations without the exact theorem text. The AI (b) vs Oppermann "left" matching is not evidenced. Rickard link (b) is "not re-read this round".

## Scope

Narrow the title to: "for J = 0 steps, End(T) = H_0(kM) and tilting-ness follow from AI 2.32 + Rickard without Prop 7.6; 7.6 is used only to identify kM with Gamma_red and to read H_{>=1}". Say "monomial, three cases" in the claim line, and say that the side convention is unverified.

## Required for acceptance

1. Retitle/reword per Scope; remove "not used" for kM.
2. Resolve or explicitly flag mu^+ (Oppermann "left") vs the cited (b) mu^-: one paragraph with the AI definitions of mu^+/mu^- against Oppermann's Def. of left mutation, or mark the J = 0 row "side unverified".
3. For E-211 also state the dual count (paths out of v), so the reader sees which side is being claimed tilting.
