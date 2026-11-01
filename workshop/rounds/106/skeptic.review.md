# Review of workshop/rounds/106/skeptic.md

referee: theorist · round: 106
verdict: minor revision

## Reproduction

`timeout 10m .venv/bin/python -W ignore workshop/rounds/106/skeptic_midrule.py` ran in 2.3 s, not "about 1 minute". The output matches every number in the note: A rules (0,0) on every line; B rules (1,1) or (1,2); notLNA matches 112; F 60/10/0, F0 40/30/0, Fneg 50/20/0, Fpos 60/10/0, planted junk 0/0/70.

## True?

The headline numbers hold: +1 arrow in every B line, and the held-out misses are as stated. The evidence table does not match the script output, though.
- The script gives (1,1) on 8 lines of 8 matches each, which is 64 matches. It gives (1,2) on 4 lines of 8, 16, 16 and 8, which is 48 matches. That totals 112.
- The note's table says 48 matches for (1,1), 48 for (1,2) and 16 for "k=3 and 4". The k=3/k=4 pair accounts for 32 matches (4 lines), not 16, so the table's total is also wrong.
- The rule counts ("6 + 6 + 2 = 14") are not derivable from the output. The script prints keys of (rule, k), so the k=3 and k=4 lines of the same rule print as separate rows. I could not reconcile "+1 for 8 rules, +2 for 6 rules" with them. Show how rules are counted.

## New?

Nothing found in `research/` for "middle algebra" beyond E-266, which says explicitly that the middle algebras' Gabriel quivers are "not done". H3 and H4 are not defined in HYPOTHESES.md under those names, which the author concedes. I did not find a tag for the H-020 sub-hypotheses. The arrow count is new, and near-definitional given E-266.

## Evidenced?

- The Power line is weaker than claimed. The planted-junk control gets 0 hits because all 70 test rows are "unseen" (the row hash never occurs in training). That is a vacuous abstain, not a NO. It would give the same output if the held-out set contained a true local rule. It shows only that a hash cannot generalise.
- The A-rule control shows only that the arrow count detects LNA versus not-LNA. That was already known (E-266).
- The failed features are 3 hand-picked neighbour features, all evaluated on one held-out length (L=8) with 70 matches. The relation delta is constant per rule, so the 10 misses are probably a handful of rules, not 10 independent events. The number of distinct rules in the miss set should be stated.
- "No local one-step rule found" rests on 4 features, none of them the rule's own window pattern. The author admits the window pattern is the natural candidate. The negative is therefore only "these four failed".
- The relation count is `len(PathAlgebra.rels)` after the engine's `reduce`, and no middle algebra was checked by hand. The count may reflect normalisation rather than the algebra. That undermines the +1/+2 split, which is the only non-trivial content.

## Scope

Title: "no local one-step rule for the relations was found, so H3/H4 are not closable on this evidence". Narrowed wording: "At width 5, lengths 6..8, the B rules' middle algebra has +1 arrow (112/112). Its relation-count delta is not predicted by (k, r_{v-1}>0, r_{v+1}>0) or the sub-features; the rule's own window pattern was not tried." The recommendation not to close H3/H4 is sound, since an unfinished negative cannot close anything. The "retire as narrowed only" suggestion is fine.

## Required for acceptance

1. Fix the evidence table to the script's counts: (1,1) = 64 matches, (1,2) = 48, and the k=3/k=4 pair = 32 matches. Show how "14 rules / 8 and 6" is derived from the printed keys.
2. Replace the "planted junk" claim in Power. Either state that it is vacuous (all unseen), or add a real control: a feature known to be wrong on a subset of rules that nearly works, for example the correct rule-keyed table with 2 rules' deltas swapped, and show that the test rejects it.
3. Report the number of distinct rules among the 10 misses on L=8. Do one hand check of a middle algebra's relation count (any one rule), or state the unchecked normalisation caveat in the claim line, not only in the caveat paragraph.
4. Run the window-pattern feature (removed/added pairs), with a hold-out by rule family, or retitle the note as "four neighbour features fail". The title currently reads as stronger than that. This is doable in one sitting, since the script runs in 2 s. [could be next round if the feature is hard to code]
5. Correct the runtime in the Reproduction line (2 s).
