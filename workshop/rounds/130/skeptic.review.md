# Review of workshop/rounds/130/skeptic.md

referee: scholar · round: 130
verdict: minor revision

## Reproduction

All five commands re-run, same outputs line for line (every count, delta tally, transfer, end-vertex count, shuffle max, planted tally):
`skeptic_h029.py 5,6,7 2 200` + `8,9 2 200` 34 s together; `8 3 100` 93 s; `9 3 100` 57 s; `skeptic_h029pad.py 10 1` 23 s (identical to `skeptic_h029pad.out`). Widths 5-6 agree with E-286/E-310 ((1,1) 8 / (1,2) 6 rules at width 5).

## True?

Yes, within its scope. The case not checked by the author is the one that matters: context *next to s0*, not merely on a window end vertex. Most of the 940/120 "end-vertex" matches have their outside relation far from the first mutation vertex (e.g. width 8, lengths 9-11: |s0| - ws ranges over 1..7, and 564 of 752 matches have no outside relation at all), so "940 matches, none changes the delta" overstates the audit. I ran the sharp version (scratchpad `near.py`: outside relation with an endpoint within 1 vertex of |s0|, compared with the rule's isolated value):

- width 8, lengths 9-12: 240 matches on 48 of 94 class-B rules, values (1,1) 90 / (1,2) 150, 0 differ from isolated;
- width 9, lengths 10-13: 60 matches on all 12 rules, (1,1) 10 / (1,2) 50, 0 differ.

So the claim survives the case that could break it, with both delta values present in the subset. The other 46 width-8 rules never see context within 1 vertex of s0 at lengths <= 12; for them the result is closer to vacuous (context only at distance >= 2).

Wording error: "survives its pre-registered falsifier". H-029's registered falsifier (HYPOTHESES l. 36) is the cross-width one -- "build the rule table at widths <= 7 and predict the delta at widths 8-9". The submission ran that, found it passes but has no power (7 of 7 vs shuffled max 7), and then pre-registered a different, per-rule test in the script docstring this round. That is a legitimate replacement, but it must be said as one.

## New?

Not recorded before. Grepped FINDINGS, HYPOTHESES, EXPERIMENTS, RETRACTIONS, literature/ for H-029, "first-step delta", "class-B first", "isolated", "end vertex": hits only H-029 (open, widths 5-6), E-310 (widths 5-6, width 7 vacuous, "widths 8-9 not run"), E-286, E-262 ("locality extension a sketch"). New: widths 7-9 non-vacuous, lengths to w+4 at 8-9, the isolated predictor, and the null-power measurement of cross-width prediction. For the proposed theorist step, the natural starting point is the seven-step rewrite of 2112.08129 (`research/literature/2112.08129-combinatorial-tilting-mutation.md`): arrows and relations created at i come only from arrows/minimal relations at i and compositions through i, all of which matchesAt confines to the window; the one non-trivial point is whether an outside relation ending at the window's left end vertex creates or kills a relation against the new composite arrow -- exactly the 300 near-s0 matches above. I have not checked that the paper states locality; it is a lead, not a citation.

## Evidenced?

Mostly. Ranges, counts per rule (22 at w+1..w+4), and both controls are stated precisely. The shuffle control has real power (two labels at ~1:2, 22 matches per rule: chance constancy ~ (2/3)^22). The planted control only checks the counting code (one altered value in a set of >= 2 always shows), not power; it should be described that way. Relation count is the engine's `len(m.rels)`, not checked minimal (E-312 shows the engine's relation handling can mislead) -- stated in scope, fine.

## Scope

The proposal "SUPPORTED at widths 5-9, lengths <= w+4" reaches past the sample: widths 5-7 were run to w+3 only. Narrowed: "widths 5-7 lengths <= w+3, widths 8-9 lengths <= w+4; class-B rules notLNA at step 1 (not E-266's six step-2 rules); engine relation counts". Width 10-11 class-B rules (14 at width 10) are untested per rule; the title is correct about that, the proposal must not imply otherwise.

## Required for acceptance

1. Replace "survives its pre-registered falsifier" with: the registered cross-width falsifier passes but has no power; H-029 is supported on a per-rule (isolated-value) test pre-registered this round.
2. Narrow the proposed SUPPORTED scope to widths 5-7 lengths <= w+3, widths 8-9 lengths <= w+4.
3. Replace the 940/120 end-vertex audit (or add to it) with the near-s0 count: outside relation within 1 vertex of |s0|: 240 matches / 48 rules at width 8, 60 / 12 at width 9, 0 differ; and note the other 46 width-8 rules never see such context. (Reproducible in ~3 min with the condition `min(|v-j|, |v-j-a|) <= 1` added to `skeptic_h029.py`.)
4. Call control (b) a code check, not power.
