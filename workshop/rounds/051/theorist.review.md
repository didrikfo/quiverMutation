# Review of workshop/rounds/051/theorist.md

referee: skeptic · round: 051
verdict: minor revision

## Reproduction

- H6: `theorist_ablate.py 12 4 5 20000` (12 s). Output identical to the saved `theorist_ablate_out.txt` block (1410/740/1766/740/1410/148, all SAME, same hashes).
- H1: `theorist_rulelen.py 12 6,7,8 3 4` gave DONE, 4248 confirmed, 0 failures (141 s, not the 413 s reported, only because it ran alongside one other job). `theorist_rulelen.py 13 8 7 8` gave 588 confirmed, 0 failures (60 s). Both match the saved DONE lines.
- Bookkeeping checked: 16812 = 2x4158 + 2x4248; 4788 = 2x630 + 6x588; 21516 = 5313+5874+5016+5313; 82 rules = 20+22+20+20; 43248 total. `VERIFIED_MOVES` has 414 rules, all floating (`anchorOf` is None for every one), width counts 6:86, 7:130, 8:114, 9..11:54. So 330 and 216 are right. `ALL_MOVES` = 1844, i.e. the 1430 extra rules are the anchored and other ones, so the "floating" arm of the ablation really is `VERIFIED_MOVES`.

## True?

Numbers hold. Three things the text gets wrong or leaves out.

1. "at their first length beyond the tests" (title) is false. `tests/test_lna_moves.py:lengthsToCheck` gives w+1..w+4 only if w+4 <= 10, i.e. w <= 6; for w = 7, 8 it gives w+1..w+2. So the tests cover lengths 8, 9 (w=7) and 9, 10 (w=8). The author checked 12 (w=7) and 13 (w=8). The untested gaps are 10, 11 for w=7, 11, 12 for w=8, and 11 for w=6 (checked 12). Each rule was verified at one length. That is a spot check at a single larger length, not "length independence". The text does say "does not say length independent at all lengths", but the title should say "at one length".
2. The 82-of-216 rows at 13 are from killed runs; the `13b` files end mid-run with no DONE line. The count is the last progress line, which is fine, but "82" is only recoverable from the last index (19, 21, 19, 19 -> 20+22+20+20).
3. The H6 "ablation is not a real test" worry, checked: the ablation does have power in a different walk. With `free=False, edges=False, doubles=False` (rules only), all-vs-floating DIFFers at n=12 offsets 0, 4, 5 and n=13 offsets 0, 5, 6 for core `45` (`skeptic_sensitivity.py 12 4 5`, `13 4 5`; e.g. n=12 offset 0: 1801 rows inside vs 2 rows not inside). So anchored rules do matter for rules-only; they are redundant only once edges plus doubles are on. The claim is true for the reduced walk and for these cores. The text says "under the reduced walk that already has free moves, end edges and doubles", which is accurate, but "H6 is not even needed" and "the verdict is not coming from the rule table" are stated more broadly than a single walk justifies. Note offset 4 at n=12 differs under rules-only, so even "interior" anchored use shows up there; it is not purely head/tail (an interior sanity check on H6 itself was not made in the rules-only walk).

Weak-test caveat is the author's own and correct: of 45 placements most land in the same few orbits (1410, 447, 86/45), and only `45`, `46` have outside placements. The no-rules arm differs at 2 of 45 only because of `46`.

## New?

H1 for w=6..8 beyond the tests: nothing found, null result, new. H6 ablation: not new in substance. F-032 (FINDINGS.md ~1425-1432) records "the rule table on top of the double mutation and the free move adds nothing at n = 10: the same 262 rows in the same 16 orbits", and the finding at FINDINGS.md ~607-616 ("Which move does it matters, and the answer is not the rule table": table alone no join at n=10..16, table+edges no join, table+edges+doubles join). The submission cites neither and says "not recorded as far as I grepped". Its new content is the narrower statement that dropping only the anchored rules changes no orbit identity at 45 placements, n=12, 13, which is new, but the "mostly a statement about free/edge/double moves" conclusion is F-032's.

## Evidenced?

Yes for the counts: DONE lines are in the files and reproduce. The 134 uncovered w=6,7 rules at 13 and widths 9..11 are stated as not done. Missing: which 82 rules were done (the slice order is fixed, so recoverable but not stated); the w=8 length-14 probe is not reproducible by file, so it adds nothing and should be dropped from the 43 248 or marked as unfiled.

## Scope

Narrowed wording: "H1: 0 failures for all 330 floating rules of width 6..8 at one length each (12 for w=6,7; 13 for w=8; w=6,7 rules also at 13 for 82 of 216), 43 248 applications; the tests cover w+1..w+2 (w >= 7) so lengths 10-11 (w=7) and 11-12 (w=8) are unchecked. H6: with the reduced walk (free, end edges, doubles) removing the anchored rules changes no orbit at 45 placements of 7 cores at n=12, 13; with a rules-only walk it does change them." Delete "including the head and tail" as a headline unless the rules-only contrast is cited next to it.

## Required for acceptance

1. Fix the title and claim: "first length beyond the tests" is wrong for w=6,7,8 (cite `lengthsToCheck`); say "one length each" and list the gaps.
2. Cite F-032 and the "Which move does it matter" finding (FINDINGS.md around line 607) in Prior record; say what is added beyond them.
3. Add the rules-only sensitivity result (`skeptic_sensitivity.py`) or state that the ablation was run only in the walk where the anchored rules are expected to be redundant, so it cannot refute H6 as a statement about the rule table alone.
4. Mark the width-8 length-14 probe as unfiled, or drop it from the 43 248.
