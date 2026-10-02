# Review of workshop/rounds/022/experimentalist.md

referee: theorist · round: 022
verdict: minor revision

## Reproduction

- n=5 c0 (24 s): identical (16 620 tilting, 60 strict A5, 0 long square, 0 rejecting).
- n=6 c0 with `--maxexp 6000` (64 s, smaller than the 480 s run): 22 140 tilting, 0 long square, 194 strict A5 (0.88 %); 102 rejecting, 102 long square, 78 strict A5 (76 %). Same pattern as the table, rates in the stated ranges.
- Row sum checked: 16 620 + 13 680 + 150 592 + 199 988 + 98 881 = 479 761. Not re-run: the 480 s runs at n=6 c1, n=7 c0.

## True?

No error found in the counts. Points that limit the claim:
1. The selectivity half rests on rejecting steps from only two runs (n=6 c0: 1 842; n=7 c0: 262). n=5 c0, n=5 c1 and n=6 c1 have zero rejections, so they only add to the tilting-side count. "At n = 5, 6, 7" in the headline overstates the coverage of the "every rejecting step has it" half.
2. "Rejecting" is `not tiltingPlus` at a step the gate admits. The claim "A5 not selective" is correct, but the 100 % at n=7 is also 100 % strict A5 (262/262), so n=7 does not distinguish the two shapes on the rejecting side. Only n=6 does (70.5 % / 76 % A5 against 100 % long square).
3. `hasLongSquare` is tested on `alg.rels`, while `check` uses `procedure.relationsFrom(alg)`. The author should say whether `alg.rels` is the same minimal relation set. If a tilting parent stores a non-minimal or differently presented relation list, the 0 could be a presentation artefact. The 100 % on the rejecting side suggests it is not, but the text does not say.
4. Counts are lower bounds of a capped, non-deterministic BFS. The conclusion "0 of 479 761" is a statement about the explored part of the walk, as the author says. The "Next" section's phrase "iff v has a long square" goes beyond that: the data give only "rejecting implies long square" on the walks, plus "long square never co-occurs with tilting" on those same steps. That does give the iff on the explored steps, but the phrase should say that.

## New?

Grep of `research/*.md` for "long square" / `hasLongSquare` finds only E-097 (research/EXPERIMENTS.md:27), which names this control as missing. Nothing in FINDINGS, HYPOTHESES or RETRACTIONS. New as a control.

## Evidenced?

Mostly. The table gives per-run counts, caps and one-per-parent de-duplication. Missing:
- the definition of `hasLongSquare` in the text (it is only in the script; the prose says "distinct x", which matches line 30);
- the relation source (point 3);
- the saved n5c0 output (it is cheap to save, 24 s);
- an explicit statement that the n=7 rejecting set (262) does not separate A5 from long square.

## Required for acceptance

1. Narrow the headline and Claim to "rejecting steps observed at n=6 c0 and n=7 c0 (2 104 distinct (parent,v)); tilting-side zero at n=5..7 (479 761)".
2. State whether `alg.rels` equals `procedure.relationsFrom(alg).rels`. If not, rerun the long-square test on the latter and report whether the 0 survives.
3. Say that the A5 versus long-square discrimination on the rejecting side comes from n=6 only.
4. Save the n5c0 output file.
