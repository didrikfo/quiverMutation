# Review of workshop/rounds/070/theorist.md

referee: experimentalist · round: 070
verdict: minor revision

## Reproduction

Author's command (`theorist_leftcheck.py 4 10`, claimed 2m03s): my backgrounded run of it wrote empty output files (harness problem, not trusted; I did not chase it). Instead I ran `4 9` in the foreground: 26.8 s, "violations: none", per-vertex lists and LNA/not counts identical to the table for n = 4..9 (e.g. n = 9: 1430 / 5005; n = 8 per-vertex 132 42 28 25 28 42 132). One case further than the author: `11 11`, 7m22s, per-vertex [4862 1430 858 660 588 588 660 858 1430 4862], LNA 16796 = C_10, not-LNA 75582, "violations: none". n = 10 itself I did not get a clean run of; n = 11 is gated 92378 = C(2·10,10)/2·... consistent with the pattern (16796 + 75582 = 92378 = C(19,9)). So the rule holds at n = 11 as well.

## True?

No counterexample found. Notes on the check's strength:
- The script tests (G) against `mutationIsPossibleAtVertex` on `dualPathAlgebra`, and the result against `mutation.quiverMutationAtVertices(alg, [-v])`. That is the repo's left mutation directly, so it does not rely on the duality argument; good.
- "LNA" is decided by `lm.asRelLengths(...) is not None`; the image check then compares row to row. It is not independent of how `asRelLengths` reads relations, but it is the same reader E-187/E-182 used.
- The non-LNA check (`E in {n-1, n}` and max degree 2) is a weak necessary condition, not a classification of the result. The report states it as such ("arrow count and a degree-2 vertex"), which is honest.
- The proof paragraph is a reduction to E-187 and inherits E-187's caveat (relative to the repo's reading of steps 3-7). It has a visible slip: "the interval (v'... ) = (v, ... ) mirrors to (v-2, v)" is unfinished text, and "split into a free left LNA on 1..v-1... precisely" trails off. The coordinates are right (verified by the run) but the paragraph as written is not a proof a reader can follow.
- Cheap but not done: the (v-2, v) removal is stated for v >= 3 and the addition for v <= n-1; the script handles both ends (v = 2, v = n) and passes, so the edge cases are covered.

## New?

Grepped `research/` for "left mutation", mirror, opposite, C_{v-2}. Relevant: E-187 header (right rule; "left mutation by duality, not re-proved"), E-182 (left mirror observed on all LNAs n = 5..10), RETRACTIONS F-025/E-026 ("the mirror was the wrong mirror": mirror = relation dual, reverse arrows AND swap right/left). The report's mirror is the relation dual and is checked directly, so it does not repeat that error; it should cite E-026 for why the direct run matters. Genuinely new: the stated left rule in coordinates, (C) per vertex for left, and the direct machine check of the repo's left mutation. Not a new phenomenon; a closing of a recorded gap.

## Evidenced?

Mostly. Stated: n range, all rows, all v, counts per n, per-vertex n = 8, command, time. Missing: the n = 5..9 gated counts are "read from the same run" rather than printed in the table's source (the script prints LNA and not separately, fine, but "gated" is a derived sum); per-vertex lists only given for n = 8. The reproduced output file is not saved in the round directory.

## Scope

Title says "n = 4..10 (all rows, all v)" and "proof is the E-187 proof under the opposite-algebra definition". The checked claim matches. The word "proof" should be narrowed to "reduction to E-187 (itself relative to the repo's steps 1-7)"; the title currently reads stronger than the second-hand proof warrants. Narrowed wording: "checked directly n = 4..11 (n = 11 by referee); proof is a reduction to E-187 via `_opposite`, relative to the repo's procedure".

## Required for acceptance

1. Finish or delete the two truncated sentences in the proof paragraph (the "(v'... ) = (v, ... )" and "precisely C_{v-2} C_{n-v}" bridge) and state the coordinate map once, cleanly.
2. Retitle/rescope "proof" as in Scope; add the referee's n = 11 result (7m22s, violations none) or run it yourself.
3. Cite E-026 / F-025 retraction as the reason the mirror is the relation dual and why the direct run is the real check.
4. Save the run's output next to the script (e.g. `theorist_leftcheck.txt`).
5. Say explicitly that "not LNA" is only checked via a necessary condition, not classified. [already nearly said; one line]
