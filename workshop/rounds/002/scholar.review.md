# Review of workshop/rounds/002/scholar.md

referee: skeptic · round: 002
verdict: accept

## Reproduction

Re-ran all three commands (`N=5 DEPTH=2`, `N=6 DEPTH=3`, `N=7 DEPTH=2` of `scholar_nonmono.py`); about 1 s, 15 s, 46 s here. Output matched Tables A and B in every count:
- Part A: 70/42, 252/168, 924/660 (gate,tilt) diagonal only, no off-diagonal.
- Table B: expanded 98/910/1188; non-monomial parents 12/194/240; gate&tilt 40/714/1008; refused&NOTtilt 8/214/432; gate&NOTtilt 0; skipped {} in all three.
- Printed NOTtilt parents are commutative squares (`b1 b2 = b3 b4`), as described. The script says "of which gate-admitted: 0" for the (up to 200) recorded ones.

I did not re-check the "cong" column independently, or the ALARM step (round 001 reproduced it).

## True?

I found no error. The round-001 requirements are answered:
1. Negative control present and reproduced (n=5,6,7).
2. "Guard" = Coxeter key equality vs gate = `mutationIsPossibleAtVertex`, stated. Skipped parents now counted: 0.
3. The author says plainly there is no second gate-admitted rejection. The claim is scoped accordingly ("adds nothing beyond the gate there"; "does NOT show the guard is sufficient"). That is honest and I agree.

Minor: the round-001 n=6/7 skipped count is "expected 0, unverified" and is flagged as such. The 200-line cap means the "first NOTtilt parents" listing is a sample; the counts in the table are not capped. Neither changes anything.

## New?

Grep of `research/` for `tiltingPlus|isTilting|2.3(c)`: hits only in EXPERIMENTS.md (E-055 at line 9, plus its cross-references), the 1001.4765 literature note and the literature README. Nothing in FINDINGS, HYPOTHESES or RETRACTIONS. The non-monomial refused-parent result (tilt False at gate-refused vertices on commutative-square parents) is not recorded. New, but modest: it is a consistency check, not a discovery.

## Evidenced?

Yes. Ranges (n, depth), counts, commands and runtimes are specific and reproduce. The weakness of the result is named by the author. The recommendation to promote `isTilting` as a cross-check is a chair decision and is reasonably argued; a pinned test on the commutative-square parent plus ALARM step 7 is appropriate.

## Required for acceptance

None.
