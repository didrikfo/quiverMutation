# H-029 per rule at widths 8-9 (lengths w+1..w+4): the class-B first-step delta is one value per rule (106 of 106 rules, 2 332 matches), equal to the rule's isolated (no-context) value, including all 300 matches with outside context within 1 vertex of s0; the registered cross-width falsifier passes but has no power; propose closing T4 with H-029 SUPPORTED at narrowed scope

author: skeptic · round: 130 · kind: result (+ proposal to close T4)
thread: T4 · bears on: H-029, H-020 (H3/H4), E-286, E-310, E-262
scope: floating class-B rules (notLNA at step 1; not E-266's six step-2 rules) of `lnaMoves.VERIFIED_MOVES`, widths 5-7 at lengths w+1..w+3, widths 8-9 at lengths w+1..w+4 (94 rules at width 8, 12 at width 9); every admissible row and window position; delta = (arrows, relations) count change of the engine's first mutation; relation count from the engine only (not re-checked against dim I - dim(JI+IJ) here); widths >= 10 only in the family check (lengths w+1, w+2), never per rule
power: per-rule test: shuffled labels (100 shuffles) keep at most 1 of 94 (width 8) and 1 of 12 (width 9) rules constant, observed 94 and 12; (a planted wrong delta, caught 59 of 59 and 7 of 7, checks the counting code only, not power). Cross-width (family / padding) prediction: NO power -- shuffled-template control reaches 7 of 7 at width 9 and 25 of 28 padded rules at width 8, so its 7 of 7 and 28 of 28 say nothing.

## Response to referee (scholar, minor revision)

All four required items done; the claim, title, *Scope* and proposal are changed to match.

1. **Falsifier wording.** Done. "Survives its pre-registered falsifier" is withdrawn. H-029's registered falsifier (cross-width:
   build the table at widths <= 7, predict 8-9) passes, 7 of 7, but has no power (a shuffled table also gets 7 of 7). H-029 is
   supported on a different test, per rule with the isolated value as the prediction, pre-registered this round in the
   `skeptic_h029.py` docstring before widths 8-9 were run. That is a replacement falsifier, and the claim now says so.
2. **Scope.** Done. Proposed SUPPORTED scope narrowed to widths 5-7 at lengths <= w+3 and widths 8-9 at lengths <= w+4;
   class-B means notLNA at step 1 (not E-266's step-2 rules); counts are the engine's; widths >= 10 untested per rule.
3. **Near-s0 audit.** Done, and it replaces the end-vertex audit as the locality evidence. New script `skeptic_near.py`
   (outside relation = a path j -> j+a not inside the window; near = min(|v-j|, |v-j-a|) <= 1, v = |s0|). It reproduces the
   referee's numbers exactly: width 8 (lengths 9-12) has 240 near matches on 48 of 94 rules, values (1,1) 90 / (1,2) 150, and
   0 differ from the isolated value; width 9 (lengths 10-13) has 60 on all 12 rules, (1,1) 10 / (1,2) 50, 0 differ. The other
   46 width-8 rules never see an outside relation closer than 2 to |s0| at lengths <= 12 (nearest 2: 14 rules, 3: 18, 4: 14).
   For them the result shows only that context at distance >= 2 is inert, which is closer to vacuous. At widths 5-7 (lengths
   <= w+3) the near subset has 12 / 34 / 46 matches on as many rules, 0 differ, and 2 / 20 / 48 rules have no near context.
   I concede the 940 / 120 count overstated the audit: most of those relations sat far from s0.
4. **Planted control.** Done. It is now described as a check of the counting code, not as power. The power claim rests on the
   shuffled-label control alone.

Not done: re-checking the engine's relation counts as dim I - dim(JI+IJ) on the width 8-9 matches. The referee did not ask for
it, and it is a separate run. The scope line already states the limitation.

## Claim

H-029's registered cross-width falsifier passes but has no power (below). On the per-rule test pre-registered this round, H-029 holds at widths 8 and 9: every class-B floating rule has a single first-step delta over all its 22 matches at lengths w+1..w+4 (all window positions, all admissible outside contexts), and the value at length w+1 -- the window is the whole line, no context at all -- predicts every longer match (1 974 of 1 974 at width 8, 252 of 252 at width 9). Context can act only from outside the window (matchesAt forbids relations that straddle it). The sharp case is an outside relation with an endpoint within 1 vertex of |s0|: 240 matches on 48 of 94 rules at width 8 and 60 on all 12 rules at width 9, with both delta values present, and none differs from the isolated value. The other 46 width-8 rules never meet context closer than 2 to |s0| at lengths <= 12, so for them the test is weaker. Width 7, vacuous in E-310 (one match per rule), is now non-vacuous (94 rules, 8 matches each, 658 of 658). Not claimed: a formula for the delta (no instrument found has power: the delta takes only (1,1) and (1,2)); anything beyond length w+4 or about the relation ideal being the true one; any statement about later steps of a rule. A refutation would be one class-B rule with two deltas at some length, or one match differing from its isolated value.

**Proposal (for the next ledger): close T4.** H-029 -> SUPPORTED at widths 5-7 for lengths <= w+3 and at widths 8-9 for lengths <= w+4 (class B = notLNA at step 1; engine relation counts; widths >= 10 not tested per rule). Reason: "the delta is the rule's isolated (length w+1) value; outside context, including relations within 1 vertex of s0, does not change it (per-rule test pre-registered r130; the registered cross-width falsifier has no power)". What is left is (i) a proof, which is E-262's unfinished class-B locality sketch at count level -- a theorist question, not more data -- and (ii) the formula question (H-020 H3/H4), which this instrument cannot answer: with two values, best hand predictors tie shuffled labels (E-310) and cross-width family prediction ties shuffled templates (below). More widths would add rows, not information. If the ledger keeps T4 open, it should be only as "theorist: count-level locality for class B".

## Evidence

Per-rule test (`skeptic_h029.py`):

| width | lengths | class-B rules | matches | one delta | isolated predicts | length transfer | end-vertex context | shuffle max const. | planted (code check) |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 6-8 | 14 | 112 | 14 | 98/98 | 28/28 | 28 | 2 (200) | 6/6 |
| 6 | 7-9 | 54 | 432 | 54 | 378/378 | 108/108 | 108 | 2 | 31/31 |
| 7 | 8-10 | 94 | 752 | 94 | 658/658 | 188/188 | 188 | 9 | 43/43 |
| 8 | 9-11 | 94 | 752 | 94 | 658/658 | 188/188 | 188 | 9 | 44/44 |
| 8 | 9-12 | 94 | 2068 | 94 | 1974/1974 | 282/282 | 940 | 1 (100) | 59/59 |
| 9 | 10-12 | 12 | 96 | 12 | 84/84 | 24/24 | 24 | 6 | 6/6 |
| 9 | 10-13 | 12 | 264 | 12 | 252/252 | 36/36 | 120 | 1 (100) | 7/7 |

Near-s0 audit (`skeptic_near.py`): outside relation with an endpoint within 1 of |s0|, compared with the isolated value.

| width | lengths | near matches | rules with near context | (1,1) / (1,2) | differ | rules by nearest outside endpoint (d: rules) |
|---|---|---|---|---|---|---|
| 5 | 6-8 | 12 | 12 of 14 | 6 / 6 | 0 | 1: 12, 2: 2 |
| 6 | 7-9 | 34 | 34 of 54 | 18 / 16 | 0 | 1: 34, 2: 20 |
| 7 | 8-10 | 46 | 46 of 94 | 18 / 28 | 0 | 1: 46, 2: 28, 3: 20 |
| 8 | 9-12 | 240 | 48 of 94 | 90 / 150 | 0 | 1: 48, 2: 14, 3: 18, 4: 14 |
| 9 | 10-13 | 60 | 12 of 12 | 10 / 50 | 0 | 1: 12 |

Widths 5-6 reproduce E-286/E-310's values ((1,1) 8 / (1,2) 6 at width 5; 28 / 26 at width 6). Delta values at width 8: (1,1) on 704 matches, (1,2) on 1 364; width 9: 44 / 220. The first step always adds exactly one arrow (as E-286).

Cross-width prediction (`skeptic_h029pad.py`), the literal reading of H-029's falsifier ("table at widths <= 7 predicts 8-9"):

| width | class-B | padding of a class-B (w-1) rule | right | shuffle max (mean) | family template seen earlier | right | shuffle max (mean) |
|---|---|---|---|---|---|---|---|
| 7 | 94 | 16 | 16 | 14 (8.1) | 6 | 6 | 3 (0.7) |
| 8 | 94 | 28 | 28 | 25 (16.5) | 7 | 7 | 7 (2.2) |
| 9 | 12 | 0 | -- | -- | 7 | 7 | 7 (3.7) |
| 10 | 14 | 0 | -- | -- | 9 | 9 | 8 (4.9) |

No wrong prediction, but at widths 8-10 a shuffled table does about as well, so this is not evidence. Most width-8/9 rules are not paddings and have no smaller-width family member: the <= 7 table cannot even be asked about them. The decisive test is therefore the per-rule one, with the isolated value as the prediction.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/130/skeptic_h029.py 5,6,7 2 200     # 9 s
timeout 10m .venv/bin/python workshop/rounds/130/skeptic_h029.py 8,9 2 200       # 27 s
timeout 10m .venv/bin/python workshop/rounds/130/skeptic_h029.py 8 3 100         # 97 s
timeout 10m .venv/bin/python workshop/rounds/130/skeptic_h029.py 9 3 100         # 67 s
timeout 10m .venv/bin/python workshop/rounds/130/skeptic_h029pad.py 10 1         # 13 s
timeout 10m .venv/bin/python workshop/rounds/130/skeptic_near.py 8 3            # 74 s
timeout 10m .venv/bin/python workshop/rounds/130/skeptic_near.py 9 3            # 42 s
timeout 10m .venv/bin/python workshop/rounds/130/skeptic_near.py 5,6,7 2        # 8 s
```
Outputs: `workshop/rounds/130/skeptic_h029_w89.out` (the two w+4 runs), `skeptic_h029pad.out`, `skeptic_near_w8.out`, `skeptic_near_w9.out`, `skeptic_near_w567.out`. Uses `workshop/rounds/095/theorist_rowmodel.py` for rows and the step classifier.

## Prior record

E-286 (width 5 per-rule, width-independent rule fails), E-310 (82 of 82 transfers at widths 5-6, width 7 vacuous, no formula beyond shuffles), H-029 (falsifier: widths 8-9), E-262 (class-B locality extension "a sketch"), E-153 (one-step locality L1, LNA results only), E-266 (A/B split). New: widths 7-9 non-vacuous, lengths to w+4, the isolated-value predictor, the near-s0 context audit, and the measurement that cross-width prediction has no power. Nothing in RETRACTIONS bears on it.

## Code changed

None (three new scripts in `workshop/rounds/130/`).

## Next

- Ledger: close T4; H-029 -> SUPPORTED at scope (widths 5-7 lengths <= w+3, widths 8-9 lengths <= w+4, class B at step 1, count level, engine relations).
- Theorist (optional): prove count-level locality of the first class-B step (the mutated algebra's arrow and minimal-relation counts depend only on relations meeting arrows at distance <= 1 of s0) -- this would turn H-029 into a lemma and finish E-262's sketch for step 1. Lead (referee, not checked as a citation): 2112.08129's seven-step rewrite; the one non-trivial point is an outside relation ending at the window's end vertex next to s0 -- the 300 near matches above.
- Not worth doing: wider widths or longer lengths; formula search on a two-valued label.
