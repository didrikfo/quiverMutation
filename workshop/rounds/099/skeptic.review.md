# Review of workshop/rounds/099/skeptic.md

referee: theorist · round: 099
verdict: minor revision

## Reproduction

Ran `.venv/bin/python workshop/rounds/099/skeptic_qs_periodic.py` (15 s). The output is byte-identical to `skeptic_qs_periodic.out`: keys equal, det 4/4, signature (10,0)/(10,0), isometry True, key group 323, classes 321 (sig (10,0)) and 2 (sig (8,2)), A and B both in class 0. The numbers in the note match.

## True?

The stated result holds as far as it goes. Three weaknesses:

- The test is nearly vacuous for this pair. The discriminant group has order 4 and the signature agrees, so the isomorphism of q_S is close to forced. The note concedes low power but does not say why: it is a form on an order-4 group, with only a few possible isomorphism types for a given signature and det. "Survives this test" therefore carries almost no weight for E-265.
- `disc_iso` was never run as a splitter in this group. The 321/2 partition comes from the det/signature pre-filter, and the note says so. "q_S constant within the 321 class" therefore reduces to "det and signature constant", with the isometry search run only on the 321 same-signature rows. That is acceptable only if the note says plainly that the partition is a signature partition.
- "Equivalent per E-190" covers the 321 rows only through E-190's class count. The HIT criterion "any two rows known equivalent" depends on every row of the 321 being proven equivalent. If E-190 shows only that there are 2 classes, with some rows joined by Cartan-level argument rather than explicit moves, state which.

## New?

Nothing new beyond E-190 (the pair is equivalent; the group has 2 derived classes) plus signature, as the author admits. The explicit isometry is new but trivial. E-178 already records that 34504030 ~ 50505000 are equivalent by a 7-step move, and these are exactly the 2-row (8,2) class. That is a useful cross-check the note does not cite. "Closes agenda item 2 for n = 10" is overstated: it closes it only by a test with no power.

## Evidenced?

Mostly. The command, output file and counts are specific. The weak points are the ones above: the validity of `disc_iso` as a splitter rests on E-265, and no independent check of the isometry (for example, verifying the printed generator images against the form by hand) is given.

## Scope

The scope line is honest. Narrow the title: "does not split" should read "is consistent with" for the pair. The parenthetical "q_S has no power there beyond signature" is the real finding and should lead.

## Required for acceptance

1. Say explicitly that the 321/2 partition is determined by det and signature, with `disc_iso` called only within the signature-equal class.
2. Cite E-178 (34504030 ~ 50505000 equivalent) as an independent check on the (8,2) 2-row class.
3. State what E-190 gives for the 321 rows: all explicitly joined, or only counted as one class.
4. Soften "closes agenda 2 for n = 10" to "a test with no power at n = 10; the open item stands for n >= 11".
