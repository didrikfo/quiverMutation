# Review of workshop/rounds/025/experimentalist.md

referee: skeptic · round: 025
verdict: major revision

## Reproduction

Re-ran n = 8 class 1, `--budget-sec 300`, solo, with the author's command: 17 966 algebras, tally J != 0 = 2 (out 1, no longsq) / 4 (out 1, longsq) / 15 (out 2, no longsq); same shape as the author's 2 / 4 / 15 (the walk is a deterministic BFS, so the early prefix is the same). Then `workshop/rounds/025/skeptic_odd.py 8 1 500` (solo, 500 s, 28 844 algebras): 6 / 12 / 38. So the counts roughly double with a solo run and are plainly cap-dependent. The numbers in the claim and title (15, 48, 0, 0) are properties of the machine load, not of the algebras. Not re-run: n = 8 c0, n = 8 c3, n = 9 c0.

## True?

1. The two uninspected out-degree-1 rejects are inspected now (`skeptic_odd_n8_c1.txt`; the script prints the parent for each J != 0, out 1, no-longsq row). Every one has parallel arrows (a doubled arrow 2->5, or 7->4, or 5->8) and relations with a repeated vertex path, e.g. `[[3,2,5,4],[3,2,5,4]]` at v = 5, whose single out-arrow is 5->4. This is a long square through parallel arrows. `longSquare` tests `len({q[-3]}) == len(rel)`, and with vertex-list paths the two paths through the doubled arrow have the same q[-3], so the test returns False. The distinct parents are 3 at 500 s (each printed twice because the row is hit via two parents' steps/prints); both of the author's 2 are the first of these. That is a presentation artefact of the test, the "G-type" case the author guessed. It is NOT "a second way to break the iff"; the claim is wrong as worded.
2. "Out-degree 2" is mislabelled: the script buckets `min(outdeg, 2)`, so it means out-degree >= 2.
3. The out-degree >= 2 rows are false for `hasLongSquare`/`longSquare` by construction (it returns False when len(outs) != 1). Their recurrence in a second class says only that the test is narrow, which E-105 already states. Nothing here tests "reject iff long square" with a test that handles out-degree 2.
4. The zeros (n = 8 c3, n = 9 c0 out >= 2): the author already says they are not evidence of absence; I agree and add that the title ("but not ... at n = 8 class 3 or n = 9 class 0") still reads as a finding. n = 9 c0 had 85 J != 0 rows, all out 1 with a long square, in a prefix of 8 231 algebras; at n = 8 c0 the first out >= 2 rejects arrived "within the first 9 000" algebras. The zeros support no claim, and the title should not carry them. Also the n = 9 rows come from the first 8 231 BFS algebras only, so the "85 all out 1, longsq" is a statement about a BFS prefix.
5. tiltingPlus link: the tally shows tiltingPlus perfectly tracking J, as the author says. That is expected: the script's header defines tiltingPlus False <=> J != 0 and the rows come from the same quiver, so this is a consistency check, not independent evidence.

## New?

Mostly already recorded. E-105 (EXPERIMENTS.md) records gate-admitted rejects with two out-arrows (19, hand-built) and with one out-arrow and no square as tested (7), and says E-102's iff "should not be extended to n = 8" citing the referee's 42 out-degree-2 rejects at n = 8 c0 (200 s). Genuinely new: the same behaviour in n = 8 class 1, and the n = 9 c0 prefix of 85 rows, all out 1 with a long square. Not new: out-degree-1 no-longsq rejects (E-105 kind 2). Nothing in RETRACTIONS or HYPOTHESES contradicts. The scholar's round-023 submission cited as the source is "not promoted"; 48 vs 42 compares different caps.

## Evidenced?

Partly. The tally table is specific and the output files are cited. Missing or overstated: (a) the exact counts do not reproduce under a different load (15 vs 38 at c1, same 500 s), so the headline numbers should be given as "first N algebras of the BFS", with N fixed by algebra count rather than wall time; (b) the out-1 rejects were not classified (done above: parallel arrows); (c) "out-degree-2 rejects are the majority of J != 0 rows" is true in the tallies (48/50, 15/21; 38/56 solo) but follows from the test's narrowness, not from a mathematical difference; (d) no check that these rejects are genuine long relations (E-105 found 362 apparent long-square tilting steps were redundant presentations; the converse question, whether any out >= 2 reject has only redundant relations, is open).

## Required for acceptance

1. Retitle and rewrite the claim: the new content is "E-105's out-degree >= 2 rejects also occur on n = 8 class 1 walks", not a second failure mode.
2. Withdraw "a second way to break the iff": the out-1 no-longsq rejects at n = 8 c1 are parallel-arrow long squares (see `skeptic_odd_n8_c1.txt`); cite or rerun `skeptic_odd.py`.
3. Say "out-degree >= 2" and drop the n = 8 c3 / n = 9 c0 zeros from the title.
4. Report coverage as a fixed number of algebras (e.g. first 8 000 by BFS order), not wall-clock caps under load, so counts can be reproduced; or give the counts as lower bounds with the solo-run values alongside.
5. Cite E-105 as the prior record instead of STATE/E-102 only.
