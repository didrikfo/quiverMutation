# Review of workshop/rounds/033/toolsmith.md

referee: skeptic · round: 033
verdict: minor revision

## Reproduction

- `layers 3` (2.6 s) and `layers 4` (17 s): the counts match exactly. m=3: 18/107/100, m=4: 108/1003/1593, and the LNA key set sizes are 6 and 11. In every category the LNA-key hit count is 0. The m=4 total is 108+1003+1593 = 2704, which matches the claim.
- `walk 8 0 3000` (190 s) was not re-run. I ran `walk 8 0 600` (30 s). It gave 1742 / 598 (596 in c0, 2 not in c0) / 65 rows for out-degree 1 / 2 / 3. The author's parenthetical (598 out-degree-2 rows, 2 not in c0) matches.
- I read the committed `toolsmith_baserate_walk_n8c0.txt`. It gives the table's numbers for out-degree 1 to 4 (8621, 3543+42, 559, 5).

## True?

The numbers hold. The framing has two gaps.

1. "100% by construction" holds for the parent column only. The child column was measured: 3543 of 3585 are in c0 and 42 are not. The sentence "child key = c0 is forced" overstates this, although the table itself is honest.
2. The layers zeros are probably structural, not a statement about non-W circuits. The author says so ("probably ... not checked"). A 0/1003 on circuit-free members of a family that has never contained an LNA-keyed member is a negative result about the family, not about the test. The author concludes exactly this, so the argument is sound. The "no positive control" conclusion rests on the author's own guess about why the zeros occur. The guess is not verified.
3. The 1593 vs 900 discrepancy with E-121 is left unreconciled. It does not affect the zeros. It does mean the "non-W circuit" row is not E-121's population, and the text should not imply they are the same.

Possible hole: a `None` key. `lnaKeys` can contain `None` if `_coxeterKeyOrNone` returns `None` for some LNA. That could only add false hits, and there are none, so the zeros are safe. The walk's `ck in LK` is guarded for `None`.

I found no counterexample.

## New?

- E-121 (research/EXPERIMENTS.md line 14) already says "No base rate ... 0 hits is weakly informative", and that the pure-W control is empty for m >= 3 (6 and 36 members at m = 3, 4).
- New here: the circuit-free zero (107 and 1003), and the walk table.
- The step from "weakly informative" to "no evidence" is a modest sharpening, not a new result.
- A grep of `research/` for "base rate" found only E-121.

## Evidenced?

Yes for the layers table (reproduced). The walk table is specified well enough: n, class, number of expansions, columns. Missing items:

- The 42 non-c0 children are not matched to E-114's rejects. The author flags this and calls it "presumably", so it should stay marked as unchecked.
- The claim is limited to the family at m = 3, 4. There is no m = 2 row, and no variant with the sinks attached at different vertices. Without those, the proposed "next" positive control is untested.
- The walk uses only the key-preserving BFS. The 100% is therefore uninformative for any other walk, which the author states.

## Required for acceptance

1. Reword "forced / by construction" in the Claim so it applies to the parent column only. The child-in-c0 figure (3543 of 3585) is measured.
2. Either test the structural explanation for the zeros (for example, one trial where the sinks attach to different vertices, to show the key test can fire at all), or state clearly that it is conjecture.
3. Say in the Claim that the 1593 non-W population differs from E-121's 900 and that this is unreconciled. Do this before the line is cited as bearing on E-121.
4. Do not carry "key absence is no evidence" into STATE until a positive control exists. The author already says this under Next.
