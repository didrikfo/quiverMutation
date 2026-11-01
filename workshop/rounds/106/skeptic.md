# The 14 width-5 class-B rules: the first step adds one arrow (also at widths 6 and 7), and the relation-count delta is fitted 14/14 on width 5 by the number of added window relations, a post-hoc rule that FAILS at width 6; no width-independent one-step rule found, so H3/H4 are not closable on this evidence

## Response to referee

1. Table fixed (see Evidence). Per rule index in the 24 floating width-5 rules, from `skeptic_midrule2.py`: 14 rules have notLNA matches (indices 0,1,4,5,6,7,10,11,12,13,16,19,22,23), 8 matches each. Delta (1,1): 8 rules, 64 matches (4,5,6,10,12,13,22,23). Delta (1,2): 6 rules, 48 matches (0,1,7,11,16,19). Total 112. The delta is constant per rule (asserted in the script). The old "6 + 6 + 2" and "16 matches for k=3 and 4" were wrong: I had merged rules sharing a window pair and counted the printed (rule, k) lines. Rules 12/22 are the k=4 versions and 13/23 the k=3 versions of the two window pairs; each of the four is its own rule with 8 matches, all delta (1,1).
2. Power: the "planted junk" control is withdrawn as vacuous (all 70 L=8 rows unseen, it would say NO to anything). Replaced by controls with seen rows: the rule-keyed table gets 70/0 hit/miss on L=8; the same table fitted with 2 rules' deltas swapped (three different swaps) gets 60/10 and is rejected, as are F0 (k) at 40/30 and the constant baseline at 40/30. These are "nearly works" wrong rules; the test rejects them.
3. The 10 misses of F on L=8 come from exactly 2 distinct rules (13 and 23, 5 matches each: the two k=3 rules on the window pair (0,3),(1,4) <-> (0,4),(2,3); their k=4 twins 12 and 22 are predicted correctly), so the "10 independent events" worry is right: the F failure is 2 rules. Relation count check: (a) by hand on one case (rule 0, L=6, row [2,3,0,0], vertex -5): parent rels 123, 2345 (2 relations); engine middle algebra has arrows 1-2,2-3,2-5,3-4,5-4,4-6 and rels 123, 125, comm(2-3-4 = 2-5-4), 546 (4 relations): delta (1,2). Longer paths agree (1234 = 1254 = 0 via 123 / 125; 2346 = 2546 = 0 via 546), so none of the 4 is redundant by hand and none is missing among those I traced. (b) Independent check on all 112 matches: minimal relations computed as dim I - dim(JI+IJ) from the engine's generators, by linear algebra on the path basis, equals `len(rels)` in 112 of 112, so the count is not a normalisation artefact of redundancy. Not checked: that the engine's generators generate the true relation ideal (that needs the tilting endomorphism algebra independently). That caveat stays in the claim line.
4. Window-pattern feature run, hold-out by length and by rule (leave-one-rule-out, LORO, 14 held-out rules). Pattern features (nrem, nadd) = numbers of relation entries in the removed and added window lists: LORO 14/14 and L=8 70/0. Which part carries it: nadd alone (nadd = 3 gives (1,2), nadd = 2 gives (1,1)), i.e. relation delta = nadd - 1. Other pattern features: (nadd - nrem, k) and (len add - len rem, k) get 6/14 LORO (no better than the constant baseline 8/14) and 50/20 on L=8. BUT the feature was chosen after seeing the table (post-hoc, 2 classes, 14 rules), so I tested it out of sample on widths 6 and 7 with the prediction fixed beforehand (`skeptic_midrule3.py`): it FAILS. Width 6 (lengths 7, 8; 54 notLNA rules): 84 hits, 78 misses of 162. Width 7 (length 8 only; 94 notLNA rules): 12 hits, 82 misses of 94. So the width-5 fit is a coincidence of small width or needs a width-dependent form. The arrow claim does hold out of sample: +1 arrow in 162 of 162 (width 6) and 94 of 94 (width 7, L=8 only), 112 of 112 at width 5.
5. Runtime corrected (about 2 s per script).

Net change: title, claim and scope are rewritten. Narrower on the relation delta (a rule exists at width 5 but is not width-independent), wider on the arrow (+1 at widths 5, 6, 7 for notLNA first steps). Notebook: not changed (no notebook entry of this note was made by me this round).

author: skeptic · round: 106 · kind: negative
thread: T4 · bears on: H-020 (H3 completeness, H4 mechanism as worded in STATE.md; I could not find H3/H4 defined elsewhere)
scope: arrows: widths 5 (lengths 6..8), 6 (7..8), 7 (length 8 only); relation counts: width 5 fitted, widths 6, 7 as out-of-sample test. First mutation only; every match of every floating rule of `lnaMoves.VERIFIED_MOVES` with an engine-notLNA first step (width 5: 14 B rules, 112 matches; 10 A rules, 80 matches). Gabriel-quiver counts of the middle algebra via the engine (`len(rels)` after `reduce`; redundancy checked independently, generation of the true ideal NOT checked); not its derived class; no step 2.
power: held-out length with SEEN rows (rule-keyed table 70 hit / 0 miss; same table with 2 rules' deltas swapped 60 / 10, rejected; (k) and constant baselines 40 / 30, rejected); leave-one-rule-out over the 14 rules; out-of-width test (widths 6, 7) of the one post-hoc rule found, which it kills. Earlier "planted-junk" control withdrawn (vacuous: all rows unseen). Weakness: the width-5 rule was found after seeing the data; the width 6 / 7 samples cover lengths 7..8 and 8 only.

## Claim

For first steps that leave the LNA class (engine notLNA, r_v = k >= 3), the engine's middle algebra has exactly one more arrow than the parent: 112 of 112 at width 5 (14 B rules, k = 3 and 4), 162 of 162 at width 6 (lengths 7, 8; 54 rules), 94 of 94 at width 7 (length 8). The A rules (k = 2, LNA): 0 extra, 80 of 80 at width 5. That is near-definitional given E-266. The relation-count delta is constant per rule and takes two values at width 5: (1,1) for 8 rules (64 matches), (1,2) for 6 rules (48 matches). A local neighbour-feature (k, r_{v-1} > 0, r_{v+1} > 0) fitted on L = 6, 7 gets 60 right and 10 wrong on L = 8, and the 10 misses are 2 rules. The rule's own window pattern does fit at width 5: relation delta = (number of added window relations) - 1, 14 of 14 leave-one-rule-out and 70 of 70 on L = 8. That rule is post-hoc and FAILS out of sample: width 6 gives 84 right, 78 wrong of 162; width 7 gives 12 right, 82 wrong of 94. So there is no width-independent one-step rule for the relations that I found; I do not claim none exists. Caveat in the claim: relation counts are the engine's `len(rels)`; redundancy was checked independently on all 112 width-5 cases, but that the engine's relations generate the true ideal was not checked. Refuted by: any window-derived feature with 0 misses at widths 5, 6 and 7 together.

Recommendation for the round-108 ledger: do NOT close H3 (completeness) or H4 (mechanism) on the strength of the B rules. Record "B rules have a gated, non-LNA first step with exactly one extra arrow, at widths 5, 6, 7 (E-266 + this note); middle-algebra relation counts: no width-independent local rule found (a width-5-only fit by nadd - 1 is refuted at widths 6, 7)". If the chair wants to retire T4's last open item, retire it as "narrowed only", not "settled".

## Evidence

Per rule (index in the 24 floating width-5 rules, as printed by `skeptic_midrule2.py`; 8 matches each at lengths 6..8; delta = middle algebra minus parent as (arrows, relations)):

| (arrows, rels) delta | rules (index) | k | matches |
|---|---|---|---|
| (1, 1) | 4, 5, 6, 10, 12, 13, 22, 23 (8 rules) | 3 for 4,5,6,10,13,23; 4 for 12, 22 | 64 |
| (1, 2) | 0, 1, 7, 11, 16, 19 (6 rules) | 3 | 48 |
| total | 14 rules | | 112 |

The relation delta is not a function of k alone (k = 3 gives both), so (k) is a weak rule. Rules 12/22 (k = 4) and 13/23 (k = 3) share window pairs; each is its own rule here, so the old "k = 3 and 4, 16 matches" row was an error (it is 4 rules, 32 matches).

Held-out length test (train L = 6, 7; predict L = 8; 70 L=8 matches):

| feature | L=8 hit / miss / unseen |
|---|---|
| rule-keyed table (positive control, seen rows) | 70 / 0 / 0 |
| rule-keyed table, 3 different pairs of rules' deltas swapped (negative control) | 60 / 10 / 0 each |
| (k, r_{v-1}>0, r_{v+1}>0) | 60 / 10 / 0 (the 10 are rules 13 and 23) |
| (k, r_{v+1}>0) | 60 / 10 / 0 |
| (k, r_{v-1}>0) | 50 / 20 / 0 |
| (k), constant | 40 / 30 / 0 |
| window pattern (nrem, nadd), (nrem, nadd, k), (nrem, nadd, lenrem, lenadd) | 70 / 0 / 0 |
| window pattern (nadd - nrem, k), (lenadd - lenrem, k) | 50 / 20 / 0 |

Leave-one-rule-out (rules right of 14): (nrem, nadd) 14, (nrem, nadd, k) 14, 4-feature pattern 14, (k, nadd - nrem, a>0, b>0) 14; (nadd - nrem, k) 6, (lenadd - lenrem, k) 6; (k) 2; constant (1,1) 8. The working part is nadd alone: nadd = 3 gives (1,2), nadd = 2 gives (1,1).

Out of width (`skeptic_midrule3.py`, prediction (1, nadd - 1) fixed in advance): width 5 112 of 112 hit; width 6 84 hit / 78 miss; width 7 12 hit / 82 miss. Arrow +1: 112/112, 162/162, 94/94.

Relation-count checks: hand check of rule 0 (L = 6, row [2,3,0,0]) in the Response; independent minimal-relation count (dim I - dim(JI+IJ)) equals `len(rels)` in 112 of 112.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/106/skeptic_midrule.py    # about 2 s
```
`skeptic_midrule.py` (original, per-(rule,k) table and neighbour features, ~2 s); `skeptic_midrule2.py` (per-rule table, misses by rule, independent relation check, swap controls, pattern features LORO, ~2 s); `skeptic_midrule3.py 5,6,7 8` (out-of-width test, ~2 s). It reuses `workshop/rounds/095/theorist_rowmodel.py` (row model and replay), so the "match" list and the parent rows are shared with the audited E-262/E-266 work; the arrow and relation counts come from the engine (`quiverMutationAtVertices`), not from the row model.

## Prior record

E-266 (r_v split, 112 of 112, "no one-step rule beyond the first step", middle algebras' Gabriel quivers "not done"); E-262 (B needs a one-step rule for Nakayama-like algebras with out-degree 2); E-187/E-205 (the LNA row rule). New here: the middle-algebra arrow and relation counts for the 14 rules, and a negative for the neighbour-feature rule. I grepped EXPERIMENTS and RETRACTIONS for "middle algebra" and Gabriel-quiver counts of B rules: no entry.

## Code changed

None in the repository code. New scripts `workshop/rounds/106/skeptic_midrule.py`, `skeptic_midrule2.py`, `skeptic_midrule3.py` only; no tests touched.

## Next

Theorist: the width-5 pattern fit (relation delta = nadd - 1) is a coincidence or has a width-dependent form; the width 6 misses (e.g. a window with 4 added relations giving (1,2), where nadd - 1 predicts (1,3), and nadd = 3 giving (1,1) where it predicts (1,2)) suggest the delta is not a function of nadd alone; the removed/added relation lengths at width 6 are the next candidate. Not done: independent check that the engine's relations generate the true ideal; step 2 for B rules; Gabriel quiver shape; width 7 at lengths above 8.
