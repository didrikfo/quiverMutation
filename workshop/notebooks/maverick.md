# Maverick's notebook

## What I now believe (after round 047)
- T1/T2 (H-021', s = n - k(c)): should be closed as description. Literal H-021 false on all readings (E-056); pairing holds for cores that pair (E-060/E-062), centre formula survives null (E-067); k(33x) = 2x lower bound derived (E-065), aax drift in E-071; k(c) in general has no rule after three fits (E-060, E-061, E-068). Round 047 was record-only, no run.
- The pairing is a property of the forward reduced walk, not of derived classes (E-065 limits); its theory waits on T10/T5. T3/T8 closes with it (same equal-size singleton pair phenomenon).

## S-1 lone 3, key level (round 042)
- Free-end K-threshold law, three lengths: n = 11 K >= 4 holds, K >= 3 fails; n = 12 same (one orbit, 2746 LNAs); n = 13 K >= 4 fails (E-139, two orbits). K0 = 3, 4, 5 first fail at n = 11, 13, 15 (15 unrun at class level).
- The failure is the lone 3 with h != K, not orbit structure; self-mirror (4,4) never fails.
- Compatible deletion rule for lone relations: delete from the shorter free side; untested for other cores.

## What I tried
- r047: reading E-059..E-071 for T1/T2, closure note. r042: maverick_n12.py. r039: n = 13 class. r037: mirror redo, single-relation keys.

## Watch for
- Image comparisons by key, not label: "differs" sound, "same" not.
- Orbit within key class = lower bound for derived class.
- One class per length is a small sample; lone 3 is the simplest core.

## Next
- S-1 n = 15 K0 = 5: build from the lone-3 orbit, `--plan` first (requested in STATE).
- Non-lone cores (pairs, lone 4 / 7); K = 2 failures at n = 11.
- Depth 7 at n = 9 for H-017 (T6 positive control, n = 12 signature).
