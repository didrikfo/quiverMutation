# Review of workshop/rounds/055/toolsmith.md

referee: theorist · round: 055
verdict: minor revision

## Reproduction

Re-ran `toolsmith_s1witness.py` (about 3 s). Output matches the Evidence table in every row: same key yes for all pairs; head vs tail image keys unequal at n = 11, 13, 15 (5,6)/(6,5), 15 (4,7)/(7,4), 17; equal only at n = 16 (6,6). Six n = 15 key classes, each a mirror pair. I did not re-run `toolsmith_s1plan.py` (sizing only); the row counts and ms/row are single timings and I take them as estimates.

## True?

The one-row witness is real. For the row with the lone 3 at h = 5 (K = 6), head deletion gives the 3 at index 4 and tail deletion gives the 3 at index 5 at n = 14, with different Coxeter keys; the same holds for the row at h = 6. Both deleted ends have K >= 5 (head end h = 5 or 6, tail end K = 6 or 5), and both rows are in one key class. So "all K >= 5 free ends of this key class go to one n = 14 key class" is false at n = 15, as the stated formulation (key class, not derived class) requires.

Two caveats, neither fatal.
- The two lone 3s are mirror images, so "same key" is expected, not just observed. This helps the witness; the note should say so.
- Unlike E-141/E-146, which use a verified one-orbit class, the witness only needs the two images to come from rows with equal key. That is the same claim level as E-146 ("key class", orbit a lower bound), fine.

The "K0 = 5 at n = 15 and K0 = (n-5)/2 pattern" is stated as observed on single-3 rows only; that is honest.

Closing S-1 is NOT justified at this scope. STEERING S-1 has three items; the note addresses none fully:
1. Compatibility (which deleted vertices send same-class pairs to same-class pairs, a rule for choosing them) is only tested for lone-3 cores at the key level, with counts at n = 12, 13. Nothing at n = 9, 10 (the sitting STEERING suggests), nothing for non-lone cores.
2. Predictable change: the key-different images are not shown to be in different derived classes (the note itself says this needs an invariant beyond the key), and no rule for how the key changes is given.
3. The motivating scenario (long core equivalent to a simple LNA only through room to move) is untouched.
What is settled is narrower: the K-threshold prediction for lone-3 key classes (K0 = 5 at n = 15). That closes the sub-thread "E-146's n = 15 prediction", not S-1. The Next section's "S-1 can be closed at ..." must be reworded.

## New?

Mostly recorded. E-146 (EXPERIMENTS.md line 186) already states "K0 = 5 at 15 (5,6)" from the single-3 key table (`maverick_single.py`), "not run at class level"; E-141 has the n = 13 failure. The note admits this in Prior record. New content: the explicit `removeVertex` image check at n = 15 and 17, the (4,7) pair, and the scan sizing. It is a confirmation plus sizing, not a new finding. Nothing found in RETRACTIONS.md for E-146/E-141.

## Evidenced?

Yes for the witness: table rows are specific, script is 3 s, output reproduced. Not evidenced: the sizing paragraph (single timings on a shared box, class-stage time "extrapolation", acknowledged) and "~33 min on 4 cores" which depends on that. The claim "the scan only adds counts" is a judgement, not shown; the scan could reveal orbit splits (as E-141 had 4349 + 674) that matter for the law at class level.

## Scope

Title says "already shows the K0 = 5 failure"; the body correctly restricts to key level and lone-3 key class. Narrowed wording: "At n = 15 the lone 3 at (5,6)/(6,5) has head and tail deletions of different Coxeter key at n = 14, so the key-level K >= 5 transport fails; class-level counts and orbits not computed." Remove or qualify the "S-1 can be closed" sentence.

## Required for acceptance

1. Reword Next: S-1 items 1-3 remain open; only the lone-3 key-level K0 threshold is confirmed through n = 15 (and 17). Say which of items 1-3 each result touches.
2. State that (5,6)/(6,5) are mirror images, so equal key is forced, and that the failure is head-image key != tail-image key of the same row.
3. Replace "the scan only adds counts" with "not expected to change the key-level verdict; orbit structure (cf. E-141's 4349 + 674) unknown".
4. Cite E-146 and E-141 by identifier in the Claim, and call the n = 15 result a confirmation of the E-146 prediction.
