# Correction to E-262: for 298 of its 304 class-B rules (not 304) the first mutation leaves the LNA class; the other 6 (widths 6, 7) leave at step 2, and in width 5 the first-step split is r_v = 2 (A, 80/80 matches) versus r_v >= 3 (B, 112/112)

author: theorist · round: 098 · kind: result (small) + negative
thread: T4 · bears on: H-020 (H1 of E-159), E-262, E-187, E-205
scope: the 24 floating width-5 rules of `lnaMoves.VERIFIED_MOVES` (14 class B, 10 class A) at lengths 6..8, every row, every match (112 B matches, 80 A); the all-width correction covers all 414 floating rules at lengths w+1..w+3 (the E-262 range; the first table below is w+1..w+2, the 6 exceptions and the 298/304 count were re-run at w+1..w+3). Row model of round 095 (checked there against the engine on 9 396 cases, n <= 8).

## Response to referee

1. A side of r_v (required 1): done. `theorist_exceptions.py` part 1 prints, per rule, the r_v at step 1 and the status. 10 A rules, 80 matches: r_v = 2 and status ok in 80/80. 14 B rules, 112 matches: r_v >= 3 and status notLNA in 112/112 (r_v = 3 for 13 rules, r_v = 4 for the rule (0,3),(1,4) -> (0,4),(2,3) with sequence 2,-4,-5; in each rule r_v is one constant value over all 8 matches, so it is a property of the rule's window, not of the row). So "A iff r_v = 2, B iff r_v >= 3" is now checked in both directions for width 5, lengths 6..8.
2. The 6 exceptions (required 2): listed under Evidence with step-1 r_v and the step-2 reason. In each, step 1 is an LNA mutation at r_v = 2 and step 2 sits at a vertex with r_v2 = 3 (outcome (1, notLNA) in every one of the 48 matches). Same mechanism as width 5, one step later. Extended to w+3 (checked: 304 B-by-replay rules at w+1..w+3, 298 leave at step 1, 6 at step 2, same 6 rules, no new exception, no mixed rule). This CORRECTS the headline of E-262 ("for the 304 class-B rules the first mutation leaves the LNA class"): it holds for 298, and is false for 6. I do not file a RETRACTION (STATE/EXPERIMENTS are not mine to edit); I ask the referee/curator to amend E-262's title to 298 of 304 with the 6 at step 2.
3. Retitled (required 3): the title now leads with the correction.
4. Definitional (required 4): agreed. Class B is defined by a step leaving the LNA class, so "first step leaves" is a restatement for 298/304 rules. The findings are (a) the 6 exceptions, (b) the reason r_v >= 3 versus r_v = 2 at the leaving step, which is what E-187's gate says directly (the step stays LNA iff r_v = 2 or v = n-1; gate fails iff a relation ends at v+1). "Window-pattern characterisation of the rule" is withdrawn as wording: the only shown fact is that r_v is constant over the matches of each width-5 rule (24/24), not a derivation from the window pattern.
5. Why w+3 was dropped from the first table: time, not principle (the 414-rule script is about 1 minute per extra length); the 6 exceptions and the 304 split were done at w+3 as asked. The 414-row table was NOT re-run at w+3; its 408/414 is lengths w+1..w+2 only.
6. Optional "next" item (Gabriel quiver of the 14 middle algebras): not done, [next round].

## Claim

Criterion, fixed before running: E-187 (right) / E-205 (left, mirror) "explains" the 14 B rules if, at every applicable match, its verdict on the FIRST mutation (gate; lands in an LNA or not) equals the engine's, and the rule-level class (A: all steps stay LNA, B: first step leaves) is determined by that verdict. Hit = 0 disagreements on both; miss = any.

Result: HIT on both counts, but the hit is narrow and largely definitional (class B is defined by a step leaving the LNA class); the content is the r_v reason and the 6 exceptions.
1. Step 1 verdict: E-187/E-205 equals the engine in 112 of 112 B matches (8 per rule): the first mutation is gated and not an LNA. In every one the row entry at the mutated position is a relation of length >= 3 (r_v >= 3, mirror for left), never 0 and never v = n-1. The 10 A rules have r_v = 2 in all 80 matches (printed, `theorist_exceptions.py`). So for width 5: B iff the first mutated vertex starts a relation of length >= 3; observed at lengths 6..8, constant per rule, not derived from the window pattern.
2. All 14 B rules replay NOTHING: the question "which replay" has the answer "none"; E-187 says which leave the LNA class (all 14, at step 1) and gives no information after that. Their second step (and the return to an LNA) acts on a non-LNA algebra where E-187 is silent.
3. Over all 414 rules, first step alone predicts the replay class for 408 (lengths w+1..w+2). The 6 exceptions (width 6: 2, width 7: 4) have an LNA first step (r_v = 2) and leave at step 2 (r_v2 = 3); 48 matches. Hence E-262's "for the 304 class-B rules the first mutation leaves the LNA class" is true for 298 and false for 6 (checked at w+1..w+3, the E-262 range). "B = first step leaves" holds for width 5 and widths 3, 4, 8..11, not for widths 6, 7.

Not claimed: any proof that B rules satisfy H1; any rule for the middle algebras; the r_v >= 3 characterisation beyond width 5 and lengths 6..8. Refuted by: a width-5 rule with a match whose first step is an LNA but which the replay classes as B (or the reverse).

Note on a pitfall I fell into first: `applyAt` returns the sequence SHIFTED to the window start; using the rule's own `d[3]` gave false "gate" failures at every non-aligned match. The first run of `theorist_b5.py` was wrong for that reason and was fixed before the numbers above.

## Evidence

Per rule (all 14 B): matches 8, E-187 == engine at step 1: 8/8, status `notLNA`, entry `>=3` in 8/8. Failing step index 0 in all. Over 414 rules (table, lengths w+1..w+2; replay class vs first-step prediction):

| width | A,A | B,B | B but first-step A |
|---|---|---|---|
| 3 | 2 | 0 | 0 |
| 4 | 4 | 0 | 0 |
| 5 | 10 | 14 | 0 |
| 6 | 30 | 54 | 2 |
| 7 | 32 | 94 | 4 |
| 8 | 20 | 94 | 0 |
| 9 | 6 | 12 | 0 |
| 10 | 4 | 14 | 0 |
| 11 | 2 | 16 | 0 |

(no mixed rule, no rule without a match at these lengths. Widths 6 and 7 B column: 54+2 and 94+4 rules; 298 = 14+54+94+94+12+14+16.)

The 6 exceptions (sequence in window labels; step 1 = LNA at r_v = 2, step 2 = notLNA at r_v2 = 3; 8 matches each, lengths w+1..w+3):

| width | rule (relations -> image) | sequence | step 1 | step 2 |
|---|---|---|---|---|
| 6 | (0,2),(2,3) -> (1,3),(2,3),(3,3) | -3,-6,-6 | -3, r_v 2 ok | -6, r_v2 3 |
| 6 | (1,3),(4,2) -> (0,3),(1,3),(2,3) | 5,2,2 | 5, r_v 2 ok | 2, r_v2 3 |
| 7 | (1,2),(2,2),(3,3) -> (0,2),(2,3),(3,3),(4,3) | 2,-7,-7 | 2, r_v 2 ok | -7, r_v2 3 |
| 7 | (1,2),(2,3),(3,3) -> (0,2),(3,3),(4,3) | 2,-7,-7 | 2, r_v 2 ok | -7, r_v2 3 |
| 7 | (1,3),(3,2),(4,2) -> (0,3),(1,3),(2,3),(5,2) | -7,2,2 | -7, r_v 2 ok | 2, r_v2 3 |
| 7 | (1,3),(2,3),(4,2) -> (0,3),(1,3),(5,2) | -7,2,2 | -7, r_v 2 ok | 2, r_v2 3 |

(pairs are (start offset, length) in the window; r_v2 is the row entry at the second mutated vertex, on the row after step 1.)

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/098/theorist_b5.py 8           # 14 B rules, step-1 verdict vs engine; under 1 min
timeout 10m .venv/bin/python workshop/rounds/098/theorist_b5_classify.py    # 24 rules, A/B vs first step; under 1 min
timeout 10m .venv/bin/python workshop/rounds/098/theorist_b_all.py          # 414 rules, lengths w+1..w+2
timeout 10m .venv/bin/python workshop/rounds/098/theorist_exceptions.py     # r_v of A and B (width 5), the 6 exceptions, 298/304 at w+1..w+3; over 2 min, under 10
```

## Prior record

E-262 / round 095 recorded the A/B split and that B's first step leaves the LNA class (grep: E-262, E-187, E-205). New here: the correction of E-262 (298 of 304, the 6 named, same mechanism one step later) and the width-5 reason r_v = 2 versus r_v >= 3. Not new: that E-187 is silent on B (E-262 says so).

## Code changed

None under `quivermutation/`. Three scripts in `workshop/rounds/098/`, importing `workshop/rounds/095/theorist_rowmodel.py`.

## Next

Theorist: a one-step rule for the middle state of a B rule (LNA plus one arrow, a vertex of out-degree 2 when r_v = 0; r_v >= 3 gives a different shape: check the Gabriel quiver of the 14 middle algebras by hand). Toolsmith: a `lnaMoves` helper returning the sequence shifted to a window so scripts stop tripping on the unshifted `d[3]`.
