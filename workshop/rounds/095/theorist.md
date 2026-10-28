# H1 for the 110 of 414 floating H-020 rules whose every step lands in an LNA is supported by an E-187/E-205 row-model replay (with a locality sketch); for the other 304 the first mutation leaves the LNA class, so E-187 does not apply

author: theorist · round: 095 · kind: proof (partial) + negative
thread: T4 · bears on: H-020 (H1 of E-159), E-187, E-205
scope: all 414 floating rules of `lnaMoves.VERIFIED_MOVES` (widths 3..11). Class A (110 rules): row-model replay at every match at lengths w+1..w+3 (all widths), at w+4 (widths 6, 7, 9, 10 re-run this round, all A/B kinds), and at w+4..13 (w<=5), ..11 (w<=7), ..12 (w<=8). Width 11 (2 A rules) only w+1..w+3 `[next round]` for w+4. The replay uses E-187/E-205 plus an **empirical** label-swap rule (1 049 engine results, n <= 8), so class A is an E-187/E-205 consequence only if that rule holds. Locality extension to all lengths is a sketch (item 3 unproved). Class B is a negative about E-187's coverage, not about the rules.

## Response to referee

Verdict minor revision; 4 required items.
1. **Retitle / rephrase: done.** Title now says "supported by replay, with a locality sketch"; scope carries the width limits (width 11: w+1..w+3 only; widths 9, 10 now to w+4); the Claim no longer says "follows".
2. **Class A widths 9..11 at w+4: done for 9 and 10, `[next round]` for 11.** New `theorist_w4.py WIDTH`: w=9 (6 A rules): 84 matches at w+4, all agree; w=10 (4 A rules): 56 matches, all agree. Width 11 does not finish inside the 10-minute limit (killed twice; `theorist_w4only.py 11`), so its 2 A rules stay at w+1..w+3 and the wording keeps that limit. Next round: run it in smaller pieces (per rule).
3. **Label-swap rule stated in the Claim as empirical: done** (see Claim, class A bullet, and scope). Class A depends on it; class B does not (failures at step 0/1, before any swap acts). Reviewer's points 3-4 (no multi-step comparison against the engine) accepted: the only multi-step check is that the row replay reproduces `applyAt`'s predicted `after` (0 disagreements), which is indirect support for the swap rule, not a proof.
4. **A/B split stable at w+4: done.** `theorist_w4.py` classifies each rule at w+1..w+3 and again at w+4 only. (kind at w+1..w+3, kind at w+4): w=6 (86 rules) A->A 30, B->B 56; w=7 (130 rules) A->A 32, B->B 98; w=9 (18) A->A 6, B->B 12; w=10 (18) A->A 4, B->B 14. No rule changes kind, no mixed rule. w+4 match counts: w=6 420 agree, 756 notLNA@step0, 28 notLNA@step1; w=7 448 agree, 1316 step0, 56 step1; w=9 84 agree, 168 step0; w=10 56 agree, 196 notLNA@step0.
Gaps 1, 2, 5 of "True?": 1 and 2 are met by the narrowed wording (the reduction to w+1..w+4 is a sketch; the all-lengths statement is not claimed as proved). 5 accepted: class B is a statement about proof coverage only (Claim already says so; sharpened below).

## Claim

Write a floating rule as (window width w, before, after, sequence). Replay the sequence with the one-step LNA rule alone (E-187 right, E-205 left), never calling the path-algebra engine. Then the 414 rules split cleanly, with no mixed rule:

- **Class A (110 rules: w=3: 2, 4: 4, 5: 10, 6: 30, 7: 32, 8: 20, 9: 6, 10: 4, 11: 2).** Every step of every match is gated and lands in an LNA, and the replay returns exactly the predicted `after` row (0 disagreements). For these, H1 at lengths w+1..w+3 (w+4 for widths 6, 7, 9, 10) is supported by the replay, which uses an **empirical** label-swap rule (n <= 8; see Evidence), so this class depends on that rule; extension to all lengths >= w+1 rests on the locality sketch below, whose item 3 is unproved. Legality (condition 2 of `verifyMove`) is then the gate itself, and condition 3 (Coxeter polynomial) is a consequence of every step being a tilting mutation, not an extra test.
- **Class B (304 rules: w=5: 14, 6: 56, 7: 98, 8: 94, 9: 12, 10: 14, 11: 16).** At every match the replay fails at step 1 (or step 2 for 16 + 32 matches of width 6/7 rules): the mutation passes the gate but its result is not a line (E-187: r_v = 0 gives n-1 arrows and a vertex of out-degree 2; r_v >= 3 gives n arrows). E-187 says nothing there, so H1 for these rules is **not** a consequence of it (a statement about E-187's coverage, not about the rules). The A/B split was computed at lengths w+1..w+3 and is unchanged at w+4 (widths 6, 7, 9, 10). The place the argument breaks: a one-step rule for mutation of Nakayama-like algebras with one extra arrow or a vertex of out-degree 2 (the gate and the next step's result on such an algebra). I do not claim B is false; it holds by E-159/E-165/E-188/E-207/E-225 evidence.

Not claimed: any statement about H2-H6, about anchored rules, or that class A is the same as "short" or "simple" rules (it is not characterised syntactically here).

## Evidence

Row model (`theorist_rowmodel.py`): LNA on n vertices = row r_1..r_{n-2}; right mutation at v is gated iff no j with j + r_j = v+1 (and v <= n-1), lands in an LNA iff v = n-1 or r_v = 2, then r_v := 0, r_{v-1} := 2; left = right on the reversed row. Verified against the repo engine on every row, n = 4..8, every v, both directions: 9 396 cases, 0 violations (gate, LNA-ness, image). New fact needed to replay a sequence (not in E-187): the engine keeps vertex labels, so later steps of a rule name labels, not line positions. Read off the engine (1 049 LNA results, n <= 8): a right mutation at position p swaps the labels at positions p, p+1; a left one swaps p-1, p. Without this the 392 "gate failures" of my first run were artefacts.

Locality argument for class A (weakest step: item 3). Take a match at window start s (arrows s..s+w-1; matchesAt forbids any relation meeting the window and leaving it).
1. A step at position p reads r_{p-1}, r_p and the relations ending at p+1. A relation ending at p+1 contains arrow p; for p inside the window it lies inside the window, so it is part of the pattern; for p = s+w it would contain arrow s+w-1.. s+w and straddle, hence is absent.
2. r_{s-1} = 0 and r_{s+w-1} = 0 at a match (a relation of length >= 2 starting there straddles). So the only spectator entry a replay can read is r_{s+w}, plus the boundary flags (v = n-1 test; s = 1, so r_0 absent).
3. Hence the outcome (gated, LNA, image, final row compared with the prediction) depends on the pattern, on [r_{s+w} = 2 or not], and on the distance of the window to the two ends (0, 1, >= 2 from the right; s = 1 or not). All of these contexts occur at lengths w+1..w+4 (r_{s+w} = 2 needs two more arrows). The step I did not prove is that the image written at r_{s-1} (when p = s) agrees with the prediction in all contexts; it is covered only by the replay.
So for class A, H1 at all lengths >= w+1 would follow from agreement at lengths w+1..w+4, plus items 1-3 (the sufficiency of w+4 is itself part of the sketch). Observed: 0 disagreements at w+1..w+3 (all 110) and at w+4..13 / 11 / 11 / 12 for the rules of width <= 5 / 6 / 7 / 8: 95 516, 10 122, 10 570, 35 888 replayed matches (0 failures). Widths 9, 10 class A (10 rules): w+1..w+3 and w+4 (84 and 56 matches, 0 failures); width 11 (2 rules): w+1..w+3 only.

Why B is not tiny: in a left rule such as (5, (0:2)(1:3) -> (0:3)(1:3)(2:3), (-5,-5)) the first left mutation sits at a vertex whose relation row entry is not 2, so the middle algebra is a line with an extra arrow. Two-step rules whose first step leaves the class are 304 of 414; this is consistent with the rule table being found by search over legal sequences, not over LNA-to-LNA steps.

## Reproduction

```
.venv/bin/python workshop/rounds/095/theorist_rowmodel.py check        # 9 396 cases, ~1 min
.venv/bin/python workshop/rounds/095/theorist_rowmodel.py rules 3,4,5,6 3   # class table, widths 3..6 (also 7, 8, 9, 10, 11 one at a time; w=11 ~ 6 min, others < 2 min)
timeout 9m .venv/bin/python workshop/rounds/095/theorist_Along.py 5 13  # also 6 11, 7 11, 8 12; each under 9 min
```
(`rules WIDTHS EXTRA`: comma list of widths, lengths w+1..w+EXTRA; widths 7..11 were run one per process in parallel, w=11 about 6 min.)

## Prior record

E-187 (one-step right rule, proved relative to `procedure`), E-205 (left mirror), E-159/E-165/E-188/E-207/E-225/E-243 (H1 checks). `research/` has no statement that H1 is a consequence of the one-step rule, nor the A/B split, nor the label-swap rule (grep "label" and "swap" in EXPERIMENTS.md found nothing). The sequences in `verifyMove` already ran through the engine's relabelling, but this is not written anywhere. The locality argument itself is elementary and may be known to the table's author (`discoverMoves` searches by windows).

## Code changed

None in `quivermutation/`. New: `workshop/rounds/095/theorist_rowmodel.py`, `theorist_Along.py`. No tests run (no library file touched).

## Next

- Theorist: extend E-187 to the one-extra-arrow / out-degree-2 algebras that appear as middle states of class B; if the gate and image there are local with radius 1, H1 is a theorem for all 414 rules. Smallest case: the 14 width-5 class B rules, first step only.
- Theorist: class A width 11 at w+4 (2 rules; split per rule to fit the time limit).
- Skeptic: the label-swap rule (1 049 results, n <= 8) is empirical; a proof from `procedure` step 7 would close it.
