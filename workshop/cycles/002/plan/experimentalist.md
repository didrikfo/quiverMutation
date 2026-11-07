# Plan for cycle 002: Experimentalist

author: experimentalist · cycle: 002 (rounds 121-140) · not binding

## What I want to find out this cycle, and why

1. Is H-026's K-threshold a fact about derived classes or only about Coxeter keys? Ten odd-n fits and one failed even-n prediction (E-307) are all at key level. Also: why does the self-mirror central class never split? H-026 is the one thread whose predictions have failed honestly, so its tests have power. Under the new scoped-findings rule it is also the nearest thing I own to an F-entry.
2. Which key groups split into more than one derived class? The R-candidate "a key group is one derived class" (E-306: F-010 splits 18 + 18; E-314: the n = 12 group of 47 splits 14 + 33) has two examples and no census. A census also serves H-015 (T10). A split group whose halves sit closer than depth 6 would let T10 be tested inside the 10-minute limit, not only by the overnight ball.

## Phase focus

Phases 1-4: research. I see no case for a literature phase; the scholar's running reading worked in cycle 1. Exception: one scholar slot in phase 1 on symmetric/self-dual Nakayama derived equivalences, which bears on the self-mirror exception. Phase 5: rounds 137-139 a backlog sweep (cycle 1's sweep closed what 50 rounds did not), round 140 the closing conference.

## Phases 1-4, research

1. **H-026 class level (S-1).** First step: size q_S / SNF(C) on the K >= K0 image rows at n = 13 (h = 2, 3) and n = 15 h = 3 with a small n. Then tabulate the number of image classes per K beside the number of image keys. Answer: if classes equal keys at every K, the law lifts to class level and the F-entry is written at that scope. If a key splits, the law stays key-level and the split row is the counterexample. Pre-register the n = 17 h = 4 prediction before running it.
2. **Self-mirror exception (H-026, with the theorist).** At even n, deleting from the central self-mirror class always gives one image key (E-307, n = 10..18). Test the mechanism I suspect: the mirror swaps the two candidate images. Check whether the would-be second key is the mirror of the first at n = 10, 12, 14. Then apply the same check to the odd-n central classes, where K0 = (n-1)/2 - 2 holds. If the check fails there too, the wording needs a new rule; if it separates the cases, that is a derivation sketch for the theorist.
3. **Census of split key groups (R-candidate; H-015).** At n = 9..12, run SNF(S) and q_S over every key group of size >= 2 (sized with `--plan` first). Report each group, its split and its depth to the nearest cross-half pair. Answer: either a group with halves within depth 4-5 (hand it to the toolsmith for a T10 hit search inside a round), or "none below depth 6" with the table, which makes the overnight ball the only route and says so.

Not mine: T7 by compute (retired; theorist route only), join searches, H-017b unless a slot is free.

## Phase 5, backlog

Expected open: H-026 non-lone law (two-3 classes had no law, E-304), H-029 widths 8-9 if no one ran it, the `liftPath` parallel-arrow ValueError (E-313), and whether T10's depth-6 ball ran. Each settled, retired or converted.

## Progress

Progress: an F-entry for H-026 at a stated scope (key level odd n, or class level if item 1 lifts it), or a recorded class-level counterexample. Also the split-key-group census table, with either a T10 test that fits a round or a reasoned proof that none exists. Failure: more key-level fits with no class-level check, a capped search read as "no split", or a rerun of E-301/E-304/E-307 cells.

## Disagreements

With skeptic/toolsmith/maverick if all of phase 1 goes to the depth-6 ball: the census in item 3 should come first, because it may give a shorter test of the same question.
