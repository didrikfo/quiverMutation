# Review of workshop/rounds/113/skeptic.md

referee: maverick · round: 113
verdict: minor revision

## Reproduction

`skeptic_homctl.py 5` (2 s) and `skeptic_homctl.py 6 12` (3 s) re-run. Output matches the table exactly (n=5: 112/112/112/112/210/210, Rtilt 70, Rk0 112; n=6: 120s, 256/256, Rtilt 64). n=4 not re-run; the totals 262 / 506 equal 30+112+120 and 40+210+256, so the Power line is arithmetically right. Extra run (mine): I dropped the `k2 != k` exclusion in P3 so the true child is compared. Result at n=5: P3 280, P3no 210, so exactly the 70 Hom-tilting steps match their own child (the C.T convention is right). The script never reports this positive control (`Rcart` is declared and stays 0).

## True?

- Pre-registered rule is coherent. YES needs per-edge Hom vanishing, K0 det +-1 and an explicit isomorphism at the meet; NO needs a recorded-power invariant or a finished proof; a failed chain search is Inconclusive. No overlap between the three outcomes.
- "Hom failure at a J != 0 step is not NO" is correct and follows E-201 (13 of 13 have Hom(T,T[-1]) != 0) plus E-212 (8 of 9 such steps joined by Hom-verified chains). The step is not a tilting exchange; that says nothing about the class.
- P1: NO by construction holds (P_b and P_b[1] both summands, identity). But the real work is an arrow, so it is a trivial control; the dead lines `Tz`, `Tz2['k']` first assignment are leftovers. Harmless.
- P2: true, K0 det 0 in 262/262, Hom-only accepts 262/262.
- "P3 is weak" is right, with a sharper reason: the child Cartan comes from the repo's own `mutation`/`reduction`, so it is not independent, and the author did not show the positive control (70/70 above, mine). Without it 506/506 NO could be a convention error.
- Limit: none of P1-P3 touches the question T10 asks ("leaves the class"); the author says so.

## New?

Grepped skeptic_tilt, homctl, duplicated summand, K0, E-201/212/295/168/161/156. Hom test and chain verifier: E-161, E-212. 13 steps all Hom != 0: E-201. Cartan 0 failures of 89: E-295. K0-unimodularity already part of the tilting test in E-287 ("unimodular K0"; generation unverified). Generation assumption: E-168. So the rule and the "Hom failure is not NO" point are a restatement of E-201/E-212; the new piece is only the planted controls, and P2 as a measured false-YES rate (262/262) for a Hom-only test. "homctl" has no hits.

## Evidenced?

Yes for the counts (re-run). Missing: the positive control for P3, and a plain statement that P1-P3 test the detector, not the claim "13 of 13".

## Scope

Fits what was checked (LNAs n=4..6, one two-term T). The title's "no power for that [leaves the class]" is supported only by the absence of a key-equal inequivalent pair at n<=7 (E-212), not by this script. Narrowed wording: "tested on small cases; controls n<=6 LNAs."

## Required for acceptance

1. Print the P3 positive control (k2 == k matches 70/70 at n=5; give n=6) so 506/506 is not a convention artefact.
2. Say P3's child Cartan is the repo's mutation output, not independent.
3. Credit E-287 for K0-unimodularity in the test, so P2 is not presented as a gap newly found.
4. Delete the dead `Tz` lines.
