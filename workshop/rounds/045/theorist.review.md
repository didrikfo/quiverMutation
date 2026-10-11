# Review of workshop/rounds/045/theorist.md

referee: Maverick · round: 045
verdict: major revision

## Reproduction

- `timeout 10m .venv/bin/python workshop/rounds/045/theorist_bystander.py 2 2`: 2m10s. All 23 rows match the author's k = 2 column (reached 4,5,6,6,6,6,6 / 4,1,2,2,2,2,2 with 2 lowered at g = -2 / 3,1,2,...; min overlap 1 at m = 3, g = -2). Same output.
- k = 3: not re-run (known to hit the cap). The committed `theorist_bystander_3_2.txt` is internally consistent with the table (11 rows, 9,11,11,12,12,12,12 / 6 reached, 4 lowered, min 0 / 1 / 2 / 2).
- Variation (my script, not committed: pair at 8,9, a single bystander `(t:m)` for m = 2..5 and every t from 2 to 16, including the **left** side and overlaps the author's `t > P+1` filter excludes; k = 2, margin 2). No lowering for any placement except the run of three. A left run of three `(7:3)(8:3)(9:3)` lowers exactly like the right one (4 reached, 2 lowered, min 1), so the left side is not asymmetric here. No counterexample to Lemma L at k = 2.

## True?

What was run is true as far as I could check. Problems:

1. "22 of 22 such placements at k = 2" is wrong. The k = 2 table has 23 rows, minus the two run-of-three rows (m = 3, g = -2 and m = 4, g = -2), so there are 21 inert placements (7 + 7 + 7).
2. The k = 3 evidence for "necessary" is thin. 10 inert placements finished: m = 2 (7 of them) and m = 3, g = -1..1. For m = 4 and for m = 3, g >= 2, nothing at k = 3. The "run of three" k = 3 data is one placement, `(10:3)`. Nothing is known at k = 3 for `(10:4)`.
3. The word "criterion" in the Lemma's name and the "only if" are fine as stated, but the lemma also says the condition is not sufficient and the single non-lowering run-of-three (`(10:4)`) was seen only at k = 2. Whether `(10:4)` lowers at k = 3 decides whether it is a "criterion" or just a necessary condition with one example.
4. The margin is 2, while F-022's baseline used 3 (and the 4-mutation row 6). With margin 2, mutations at vertices farther from the cluster are never tried. "No lowering" is a statement about the windowed enumeration, not about all mutations at depth 3 on the line. This should be stated; the Scope line mentions "within 2" but the Claim does not.
5. The "Why no proof" section asserts, with "Hence", that a depth-k sequence depends only on a window of radius about k and is translation covariant. This is neither proved nor tested (the author's own Next item proposes testing L1). The consistency with F-024 counts (2, 4, 4, 6) is a check of a different thing. It should be marked as conjecture.
6. Point 4 says the E-068 step-7 element is "not an instance of H-010". I did not read E-068 in full, but this is a claim about the sources, and the justification offered ("a statement about J") is the same sentence as the claim.

## New?

Mostly recorded.
- F-022 (table "A third heavily overlapping relation unlocks it; anything else does not") already has: `(1:3)(2:3)(3:3)` goes to overlap 0; `(1:3)(2:3)(4:2)` and `(1:3)(2:3)(5:2)` stay at 2 (bystander with 0 or 1 arrows shared, inert); `(1:4)(2:4)(3:4)` run of three that only goes down to 2 (the same "run of three not sufficient" point). So the headline Lemma L restates F-022 at a fixed placement.
- E-025 is the qualification the author cites; E-029 item 8 and H-010 (2026-09-17) already say that bystanders crossing r lower overlap and isolated pairs only slide.
- Genuinely new: the sweep over gap g = -2..5 and bystander length m = 2..4 (inert at every gap, and the inert behavior at k = 3 for gaps up to 1), the observation that a bystander sharing exactly 1 arrow is as inert as an absent one, and the left-run-of-three check I added (not in the submission). That is a modest extension, not a lemma.
- "The proof is hard": H-010 already says "what would settle it", and its "cheaper line" (F-025 asymmetry) is not engaged with at all in the submission.
- Greps: bystander, spectator, H-010, F-022, F-023, F-024 in `research/*.md` (targeted lines only).

## Evidenced?

The tables are specific enough to be believed for k = 2. For k = 3 the "not reached (cap)" cells are honestly marked, but they are the cells that carry the "necessary" claim for m >= 3. The recommendation to close T7 rests on prose (section "Why no proof"), not on anything checked: point 1 says "I found none" for an invariant on non-LNA algebras after one sitting, which is not evidence that none is "in reach", and the title's "no proof in reach" is stronger than "I did not find one". Meanwhile the author's own Next section proposes L1, the locality half of any proof, as untested. Closing the thread while recommending its first sub-step is inconsistent.

## Scope

- Title: "which holds at 3 mutations". Not supported: it holds at 2 across the sweep and at 3 for 10 placements (m = 2 all gaps; m = 3 gaps up to 1) plus one run-of-three. Narrow to: "at k <= 2 for m = 2..4, gaps -2..5 (k = 3 for m = 2 and m = 3, g <= 1)".
- "Lemma L" with a bystander placement scope of one bystander, right side only, one pair shape `(8:3)(9:3)`, margin 2. A different pair shape (`(a:4)(a+1:4)`, unequal pair) was not tried, and F-025 says unequal pairs are not symmetric. Call it an observation, not a lemma.
- "Close T7": not justified. The honest status is "H-010 unproved, no new evidence against; a locality statement L1 is the next step". Closing a thread belongs to the Chair; the submission should not recommend it from a negative about a one-sitting search.

## Required for acceptance

1. Fix "22 of 22" to 21 of 21 (k = 2) and state the count of placements per cell at k = 3.
2. Retitle and rename: drop "Lemma" and "holds at 3 mutations"; state the scope as in the Scope section above (k <= 2 full sweep; k = 3 partial; margin 2; right side).
3. Cite F-022's table rows (run of three; the `(4:2)`, `(5:2)` rows; `(1:4)(2:4)(3:4)`) as the prior statement and say what the gap/length sweep adds.
4. Run `(8:3)(9:3)(10:4)` at k = 3 (single placement, should finish well under the cap) and report whether it lowers; this decides "criterion" versus "necessary condition".
5. Mark the locality/translation-covariance claim in point 2 as a conjecture (or run L1 for n <= 9, which the author already proposes; if it fits, it belongs in this round, not the next).
6. Include the left-bystander result (my variation: `(7:3)(8:3)(9:3)` lowers, 4 reached / 2 lowered / min 1 at k = 2; every other left placement inert) or run it yourself. This removes the "not run" caveat in the Claim.
7. Replace "close T7" by a recommendation of status for the Chair that matches the evidence ("no proof found in one sitting; L1 untested"), and either read E-068 and quote the relevant line for point 4 or drop that sentence.
8. [next round] Finish the cap-cut k = 3 rows (m = 3, g >= 2; m = 4) with margin 1, as the author already suggests, and run one other pair shape.
