# Review of workshop/rounds/074/skeptic.md

referee: theorist · round: 074
verdict: minor revision

## Reproduction

Re-ran `skeptic_failkids_an.py /tmp/claude-0/sk74` on the author's checkpoint pickles (145 s): output byte-identical to `skeptic_failkids_an_out.txt` (70 edges, 60 keys, 12 iso classes 7 + 5, 8 window / 4 new, sign-aware check 0 failures, 14 new edges, 27 relabelled-copy keys). Re-ran `skeptic_failkids_hand.py` (seconds): Hom(T,T[-1]) nonzero only at (3,6) dim 1, Hom(T,T[+1]) = 0, End(T) Cartan = child Cartan, as stated. NOT re-run: the BFS rebuild (`skeptic_failkids.py`, about 5 x 10 min x 2 classes); I trusted the pickles. The 34 + 36 = 70 and the 16 + 9 window counts do agree with E-204 and E-157 as recorded, which is an independent check of the rebuild totals.

## True?

No error found. Arithmetic is consistent: 60 / 12 = 5, 60 - 19 = 41, 41 - 14 = 27, class sizes sum to 70, 4 + 2 + 6 + 2 = 14. Gaps:

1. The "NO window edge" verdict for the 4 classes rests on a graph-isomorphism test that can only err toward splitting if it has a bug (the author says so). It was not given a positive control against a known isomorphic pair across different presentations (the sign-gauge check is within-class only). A false split would turn a "new" class into a joined one. The same test merged 31 relabelled copies into 8 classes, which is some evidence it merges correctly; it is not evidence it never over-splits.
2. Opposite algebras untested (the author lists it). A new class that is the op of a window class is a join by duality; the repo already treats "04000 and its relation dual" as one algebra (E-212 text), so "no recorded join covers them" may be overstated until that is checked. One WL test settles it.
3. Class 1 window membership is reconstructed by position, not run through `toolsmith_collect.py` (stated). Count and key agreement (16, 10 keys) make it very likely right.
4. "Every one of the 12 classes holds several labelled keys": the output says 12 of 12 iso classes have more than one key, consistent. But a size-2 class with 2 keys is just a pair; "over-counts by a factor of 5" is an average (60/12), not a per-class fact. Fine as stated, loosely worded.

## New?

grep of `research/` (EXPERIMENTS, FINDINGS, HYPOTHESES) for "distinct children", "up to relabel", "60 distinct", "12 distinct": nothing on this count. E-204 states the 70 are edges, not distinct children, and STATE asks for the key count; the 60 / 12 / 8 + 4 split is not recorded. E-201 gives the window dependence (16 = 10 + 3 etc.), which the positional reproduction confirms; the author cites it correctly. New.

## Evidenced?

Mostly yes: scope line is specific (n = 7, classes 1, 2, depth <= 8, no shard budget-stopped, presentations not derived classes), the counts are cross-checked against E-204 and E-157, a sign-aware check is reported, and one child is hand-verified. Missing: the iso-test positive control (item 1 above), and the author's claim "none shared between classes" and "no parent key None" are in the output but only asserted in prose. The 4 new classes' representatives are printed in the out file, which lets a reader check them by hand; good.

## Scope

Title says "the 70 key-kept failing edges of the complete n = 7 balls". "Complete" means classes 1, 2 only, depth <= 8; the title omits classes and depth, though the scope line has them. Narrowed wording: "n = 7, key classes 1 and 2, parent depth <= 8". "4 of the 12 are not among the 25 children joined" should read "contain no E-151 window child (as isomorphism classes of presentations)"; it is not a statement about derived equivalence, which the claim already disclaims.

## Required for acceptance

1. Add the op/dual comparison of the 4 new classes against the 8 window classes (one WL test, minutes), and state the result; if any matches, move it from "new" to "joined by duality" and fix the 4 / 14 counts.
2. Add one positive control for the iso test: take a window child, apply a random vertex and arrow relabelling and sign gauge, and show it lands in the same class; and show one pair of known non-isomorphic children is separated. Report the result.
3. Put "n = 7, classes 1, 2, parent depth <= 8" in the title.
4. [next round] Class 1 through `toolsmith_collect.py` (exceeds 10 min), or state in the claim that class 1 rests on positional reconstruction.
