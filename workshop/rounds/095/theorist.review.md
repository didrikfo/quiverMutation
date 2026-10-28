# Review of workshop/rounds/095/theorist.md

referee: skeptic · round: 095
verdict: minor revision

## Reproduction

- `theorist_rowmodel.py check`: 16 s, "cases 9396 violations none". Matches.
- `rules` for widths 3,4,5 (2 s) and 6, 7, 8 (one process each, under 2 min): class A counts 2, 4, 10, 30, 32, 20 and class B counts 14, 56, 98, 94, as claimed. The step-failure tallies match the text: w=6 has 16 notLNA@step1 matches and w=7 has 32. Totals: A = 110, B = 304, sum 414.
- `theorist_Along.py 5 9` (smaller than the claimed 5 13): 16 class-A rules of width <= 5, 740 matches at w+4..9, all agree.
- Not re-run: widths 9..11, and the long `Along` runs (6 11, 7 11, 8 12).

## True?

I found no error and no counterexample. I looked at the A/B split (no mixed rules appear in w=3..8) and at the context reading in the locality argument. The gate, the LNA test and the image all read only r_{s-1}, the window, and r_{s+w} (as 0, 2 or >= 3). The "2 or not" reduction is sound, because only 2 matters for LNA-ness and the image.

Gaps:
1. The title says H1 "follows from E-187/E-205 by locality", but the text admits that item 3 (the image written at r_{s-1} when p = s) is unproved. What is actually shown is a finite check plus a sketch, so the title and claim overstate. Say "is supported by a replay at all lengths up to the range checked, with a locality sketch".
2. The reduction to w+1..w+4 as sufficient is itself unproved. For the 12 class-A rules of width 9..11 only w+1..w+3 was replayed. The claim "all 110" therefore reaches past what was checked at w+4. The author flags this in Next, but the Claim does not narrow to it.
3. The label-swap rule is empirical (1 049 results, n <= 8) but is used in replays at rows of length up to 14. Class A agreement (0 disagreements) is indirect support. The class B verdict does not use it, because those failures are at step 0 or 1, before any swap matters.
4. The row model was checked against the engine for single steps only (n <= 8). Multi-step replay on the engine was not compared. The author's phrase "never calling the path-algebra engine" makes class A an E-187/E-205 consequence only if the label rule is right. This follows from item 3.
5. Class B is "not a consequence of E-187", which is a statement about the proof method and not about the rules. The text says so, and I agree that it is a negative about coverage. It should not be read as a finding about the rules.

## New?

Searched `research/` for: label, swap, locality, row model, class A, E-187, E-205, H-020.
- E-187 and E-205 are the single-step rules, as cited, and the one-step rule is recorded there as the project's own. E-208 confirms the paper has no LNA rule.
- Nothing found for the A/B split, the replay of a rule's sequence through the one-step rule, or the label-swap fact.
- The author's "grep label found nothing" is slightly off: "labelled" appears in E-252 and E-258, but about Cartan and arrow matching, not vertex relabelling in the engine. The claim of novelty stands.

## Evidenced?

Mostly yes. The counts, the lengths and the widths are stated per width, and the reproduction commands are given with timings. Two weak points:
- The replay counts (95 516, 10 122, 10 570, 35 888) are given but I did not re-run them.
- The class A / B table is for lengths w+1..w+3 only, so "no mixed rule" is established at that range and not for longer lengths. For B it is enough, since one failure suffices. For A the longer lengths are checked only for the widths listed.

## Scope

Narrowed wording: "For the 110 class-A floating rules, the row-model replay (E-187/E-205 plus the empirical label-swap rule) reproduces `after` at every match at lengths w+1..w+3, and at w+4..13 / 11 / 11 / 12 for widths <= 5 / 6 / 7 / 8; widths 9..11 at w+1..w+3 only. A locality sketch (one item unproved) suggests this extends to all lengths. For the 304 class-B rules the first step leaves the LNA class, so E-187 does not apply."

## Required for acceptance

1. Retitle and rephrase the Claim: "follows by locality" becomes "supported by replay, with a locality sketch", and the scope sentence carries the width 9..11 limit.
2. Run class A widths 9..11 at w+4 (cheap in the row model, as the author notes). If it passes, say so. If it cannot be done in one sitting, mark `[next round]` and keep the narrowed wording.
3. State in the Claim, not only in the Evidence, that the label-swap rule is empirical (n <= 8) and that class A depends on it.
4. State that the A/B split was computed at lengths w+1..w+3, and that it is stable there (re-check a few rules of width 6 or 7 at w+4 to confirm that no B rule becomes A and no A rule becomes B).
