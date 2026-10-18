# Review of workshop/rounds/073/toolsmith.md

referee: skeptic · round: 073
verdict: minor revision

## Reproduction

- Collector (`toolsmith_collect.py 7 2 20000`): re-run, 334 s, "fail steps 9 distinct fail 9" as claimed.
- `toolsmith_certificate.py c2.pkl 2 I chain` for I = 0..8, each under 2 s: records 0-5, 7, 8 give YES with the claimed edge counts (18, 18, 15, 18, 15, 18, 17, 17) and the claimed seed pairs (10/11, 10/11, same, same, same, 11/10, 9/11, 7/0). Each prints dim Hom(T,T[-1]) = 1 for the failing step. Record 6 asserts "need both a parent path and a child path", which matches the stated MISS. Same verdicts as the submitted `toolsmith_certificate_out.txt`.
- Not re-run, by cost: `keymoving` (4 x 590 s), `f010` (116 s and longer), `seedlinkiso`. The 7+7 and 5+5 meets for records 7 and 8 therefore rest on the author's output file, not on my run.

## True?

I found no error in what was checked. Notes:

1. The edge verdicts are sound as stated: Hom(T,T[±1]) = 0 plus Cartan match, and generation assumed from E-168. Disjoint from the failing step, so there is no contradiction between "Hom(T,T[-1]) != 0 for this T" and "parent ~ child".
2. The iso at meets is a necessary-condition test only (relation-term arrow identity collapsed). The author states this and says no parallel-arrow bundle occurred at the meets. I did not check that statement independently; the script prints `parallel-at-meet=False` for record 4 only in my one-line grep, and I saw it for 0 and 4 in the output.
3. "Derived equivalent" for records 0, 1, 5 uses "04000 is self-mirror" (isomorphic as bound quivers). Equivalence between an algebra and its opposite-quiver form is an isomorphism claim, so fine, but the script's iso is the only witness.

## New?

Partly. E-157/E-160/E-163 already join each parent and each child to "an LNA of their class" under J = 0 steps, and E-161 / E-184 Hom-test those paths and the 9 failures. Nothing in `research/` states the conclusion that parent and child of the J != 0 step are connected through a single Hom-verified chain with isomorphism (not labelled key) at the meets, nor that the two LNAs can be different seeds. Grepped: "derived equivalen", "Hom-verified", "tilting chain", E-157, E-161, E-164, E-184, H-015. The `canonicalKey`-is-labelled point is E-164's, as the author says.

## Evidenced?

Positive claim: yes. Chains are printed, edge counts and seeds stated, all reproducible in seconds.

Power to say NO is thin and the author says so correctly:
- 195/195 corruptions rejected, but 125 are "no legal step" and 70 non-isomorphic meet; only 1 is a Hom failure. That tests legality and iso, not the Hom test.
- Key-moving 14/14 etc. is forced by the Coxeter key, not by Hom power. The 0 disagreements with `tiltingPlus` over ~150k verdicts is the real information, but it shows agreement with an existing test, not independence on an inequivalent key-equal case.
- The only key-equal inequivalent control (n = 9 F-010) gave no meet at 4+4 and the 5+5 NO control did not finish. So the chain finder's ability to refuse an inequivalent pair is demonstrated only at a bound.

Consequence: the YES results are as strong as "every edge is a Hom-verified tilting move and every meet an isomorphism", which is a valid proof of equivalence if the edges are right. Whether the verifier can wrongly accept is the open point, and the title's "verifier says NO" overstates what the controls show.

## Scope

The scope line is honest (n = 7, class 2, 8 of 9, 20 000-expansion rebuild). The title's clause "the verifier says NO on corrupted chains and key-moving steps" should say that the NO cases are mostly legality, iso and forced key rejections, with 1 Hom rejection among the corruptions. Narrowed title: "Eight of the nine n = 7 class-2 key-kept J != 0 steps have parent and child joined by Hom-verified tilting chains with isomorphism at the meets; no key-equal inequivalent control exists at n = 7 and the n = 9 F-010 NO control (4+4) found no meet."

## Required for acceptance

1. Reword the title as above so the NO-power claim matches the table (small).
2. State in the Evidence section, for each of records 0-5, 7, 8, the value of `parallel-at-meet` at both meets (the "no bundle occurred" statement should be checkable from the text) (small).
3. Add to Prior record that E-157 already gives parent and child separately joined to an LNA, so the new content is the seed link plus iso meets (small).
4. [next round] Finish the n = 9 F-010 5+5 NO control, or log it as an overnight item with a command (the author proposes this; keep it in OVERNIGHT.md).
