# Scholar's position (round 020)

## Most promising question

**Why does the Cartan congruence fail exactly when dim ker g_i ≠ 0?** This is the single crux: the (k,i) entry of the rewrite equals dim coker g_i when the congruence holds (n=5..7, 807 steps, 0 disagreements), and the discrepancy is always -dim ker. E-095 shows every rejecting parent has exactly one such vertex. If we can derive this from step 7's commutativity element, we turn an observation into a theorem that explains why `tiltingPlus` is sound even when the gate alone is not.

## Weakest claim

**E-095 (dim ker histogram) is conditional on A5-shaped parents, unchecked.** The 1,050+ parents are assumed to have A5 structure, but that is not verified. If they don't, the pattern may be fragile or apply only to a subclass. This is the pivot between "observation" and "principle."

## What I need

- **Theorist:** derive "(k,i) entry = dim coker" from step 7 (requested). Does the derivation hold for off-diagonal (k,i) entries? Can you show why the rewrite's (k,i) must equal coker dimension when the congruence holds?
- **Experimentalist:** A5-shape check of the 1,050+ rejecting parents (requested). This confirms whether the dim ker pattern is universal or restricted.
- **Toolsmith:** once A5 shapes are checked, hardcode the Cartan assertion in mutateAtVertex as a cheap sanity test.
