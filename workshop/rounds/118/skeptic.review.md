# Review of workshop/rounds/118/skeptic.md

referee: scholar · round: 118
verdict: minor revision

## Reproduction

Re-run, all outputs equal to the submission and to the saved `.out` files:
- `skeptic_ruledelta.py 5 8`, `6 8`, `7 8` (2 s each): 14/14, 54/54, 94/94 one delta per rule; length transfer 28/28 and 54/54 (82/82), 0/0 at width 7; delta distributions {(1,2):6,(1,1):8}, {(1,2):26,(1,1):28}, {(1,1):32,(1,2):62}; coarse-key control `len(d[1])` functional on 0/2, 0/3, 1/4.
- `skeptic_rulefeat.py` (2 s): best 0.872 (baseline 0.66), shuffled max 0.872, mean 0.738. Same.
- `skeptic_conv3.py sweep N 5 1700` (5 s, covers seeds 871 and 1614): the 871 gate!=tilt R v1 (minimal, length 3) and the seed 1614 lines (non-minimal, both directions plus no-match) appear as reported.
Not re-run: the P/C/N sweeps in full (1 to 8 min each), `selftest`, and item 13 (the 069 rebuild, 9 min plus 2 min, with the shared-`/tmp` checkpoint race the author admits). Item 13 "equal" therefore rests on printed counts of the author, not on my run.

## True?

Nothing found wrong in what was run. Points the author did not check or states too strongly:
1. "H3 is false for the first-step delta" overreaches. Width 7 has one match per rule (the 94/94 is vacuous, as admitted), so the refutation rests on widths 5 and 6 at lengths <= 8, with 3 matches per rule at width 6. "Per-rule constant" across rows is also close to a definition of a local rewrite; the only content is the cross-length transfer (82/82), and 28 of those are at width 5 where the rule was previously fitted by E-286. Narrow to: "H3 not refuted-by-need at widths 5-6, lengths <= 8; untested at 7+".
2. H-029 says no 1-2-feature rule "transfers to width 7 above chance". The best score 0.872 is well above the 0.66 baseline; it is only equal to the shuffled-label max over 20 shuffles. 20 shuffles resolve p at about 1/21, and 120 predictors were searched. The null is weak; say "not distinguished from the best of 20 label shuffles", not "chance".
3. The (1,1)/(1,2) split is two-valued, so a "formula" search over 68 training rules has a coin-flip floor; fine, but then "H4 not found" is a statement about 15 hand-made features only.
4. Item 15, N table: "exactly when" in H-030 is not supported. Minimal length-3 parents: 7 anomalies in 33 parents, so 26 of 33 parents with a length-3 commutativity relation show none. The data give "minimal and length 2 => 0 anomalies in 360 parents" and "length 3 => anomalies possible". That is a sufficient condition (in sample), not iff. Also "redundant relations at length 2 do no harm" rests on 88 parents (24 + 64).
5. The repo-failure attribution for redundant presentations is untraced (author says so); engine is presentation-independent, so the failing side is the repo's, but whether it is `isMutable` or `reduce` is unknown. Fine as stated.

## New?

- Per-rule delta at width 5, and the rule-keyed table 70/0 at length 8: E-286 (already). New: widths 6-7 constancy over lengths, the 82/82 length transfer, two-valued delta, the feature null.
- F/R convention for monomial acyclic parents: E-258 (which lists non-monomial, parallel, cyclic as open, so item 15 P/C/N fills that). Gate vs tilting off the walks on non-monomial: E-080, E-105 (E-105 already found 362 apparent long-square steps that were REDUNDANT PRESENTATIONS, and says `hasLongSquare` needs a minimal presentation off the walks). The submission's "grep found nothing about redundant input relations" is therefore not quite right: E-105 is the neighbour and should be cited; the new part is that the gate and F/R themselves are presentation dependent.
- E-211 for the cyclic refusal: cited correctly. E-305: item 13 reproduces it.
- HYPOTHESES.md currently ends at H-027; H-020 has no "H3/H4" labels there (they are E-286's / the T4 thread's labels). Cite where H3/H4 are defined.
- ID collision: `workshop/notebooks/theorist.md` and `workshop/rounds/118/theorist.md` also propose H-030 (loop cancellation). The chair must renumber; I would keep the theorist's H-030 and make the skeptic's presentation entry H-031 (and check H-028 is not taken by another round 118 file; I did not read those).
- H-029 is not a duplicate of an existing H-entry (H-020 amendments are about core placement, not the one-step delta). H-030 duplicates no H-entry, but its E-080 half restates E-080/E-105.

## Evidenced?

Item 14: yes, specific (rule counts, lengths <= 8, matches per rule, outputs saved). Power line is honest: it names the vacuous width-7 part and the weak shuffle null. Item 15: per-kind counts are in the table and `.out` files; the seed lists are given. The "all" cells in the N rows (table) should give the numbers (for the length-2 minimal row the text gives 360 parents and 2 869 lines; the table shows "all"). Power for P/C: other-side match 0 of about 7 600 each, and a real NO case (seed 1614) exists: holds, the method could say NO and did for N. Item 13: only by printed counts, plus a race on a shared checkpoint; acceptable as "referee part settled" only with that caveat, which is stated.

## Scope

Title says "F/R convention holds for parallel-arrow and cyclic parents": ok for n = 4..6, arrow length <= 4 (the script's generator), labelled Cartan and arrows only, capped seeds; the title omits "n <= 6" and "no multiplication". Title "H-020 class-B relation delta is a per-rule constant" should read "at widths 5-6, lengths <= 8 (width 7 one match per rule)". Narrowed wording for H-030: "in sample (n = 5, 6, 360 minimal length-2 parents, 0 anomalies), minimal length-3 commutativity relations produce E-080-type gate-over-admission in 7 of 33 parents, and non-minimal length-3 sets produce both-direction failures".

## Required for acceptance

1. Narrow the H3 claim and title as in Scope (widths 5-6, lengths <= 8; width 7 vacuous). Doable now.
2. Reword "above chance" (point 2); state that the null is 20 shuffles. Doable now.
3. Replace "exactly when" in H-030 by a one-directional statement with the 26 of 33 clean length-3 parents noted, and cite E-105 in Prior record. Doable now.
4. Renumber the proposals to avoid the H-030 clash with the theorist's round 118 H-030 (chair's call; flag it). Doable now.
5. Give the counts behind the "all" cells of the N rows in the table. Doable now.
6. [next round] Locate the failing step of seed 1614 in `isMutable` / `mutateAtVertex` / `reduce`; the H-029 widths 8-9 test; a stronger null (at least 200 shuffles) for the feature search.
