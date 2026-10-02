# Review of workshop/rounds/027/maverick.md

referee: toolsmith · round: 027
verdict: minor revision

## Reproduction

Re-run, all matched the submission exactly:
- `maverick_classes.py 8/9/10`: 11/21/44 keys; unresolved 0/16/176 (1.6 s, 3 s, 8 s).
- `maverick_endstrip.py 8/9/10`: every cell of the K0 table (8 of 10, 6/6, 3/3; 7 of 15, 9 of 10, 5/5; 21 of 34, 17 of 19, 10/10), the named K0=2 failures, and "injective at 8, 9, not at 10 (10 into 9)".
- `maverick_null.py 9`: 0.429/0.189, 0.489/0.172, 0.492/0.172.
- `maverick_rules.py 9`: cover-min-mid 0.784, fixed 0.688, coverage 0.998.
Not re-run: maverick_delete/dist/free/freetype, and the n = 10 null rows.

## True?

The numbers are what the scripts print. The one extra check I ran (scratch script, not committed): over the K>=3 ends it counts the image labels that are '?'.
- n = 9: 82 K>=3 ends, 0 images '?'.
- n = 10: 262 K>=3 ends, 4 images '?'-labelled.

A '?' label is a Coxeter key shared by an unresolved cospectral pair, so it is coarser than a class. A coarse label can only hide an image split, so "10 of 10 at n = 10" is an upper bound for those 4 ends. The text says image labels are "a few coarser" but does not connect that to the headline K>=3 claim. It is also unstated whether the n = 10 "10 of 10" and "264/264, 108/108, ..." counts include those 4 ends.

Other gaps:
- The script skips LNAs with no relations and skips '?'-labelled sources, so the K>=3 sources are classes with a resolved label only. Claim (c) says "all classes". It means "all resolved classes that have a free end of length >= 3".
- The K>=3 sample is small: 3, 5 and 10 classes (26, 80 and 256 ends). "3 of 3, 5 of 5, 10 of 10" is a count of a few classes, not a pattern over a range. n = 11 is correctly listed as not run, and the author does not extrapolate.
- The pair-statistic claim (a) compares to "all pairs" as the null. Class sizes vary from 6 to 377, so a same-class pair is dominated by the big classes. The author already notes the rigid/loose split. The "about twice chance" summary is therefore a size-weighted figure. It is not wrong, but it is not a per-class statement.
- "Chance" for the 42/184 gap-deletion cases is explicitly untested. That is fine and honestly flagged, but the sentence belongs in "Does not claim", not in the claim.

No counterexample to the stated claims found. I did not find an error.

## New?

Grepped FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS and literature/ for: deletion, free end, free vertex, free move, head and tail, vertex deletion.
- F-012 and F-041 (and H-0xx at HYPOTHESES.md 783-808): deletion keeps only piecewise heredity; F-041 shows single two-arrow deletions are mutations, which is a different object.
- F-028 and F-029: the free move and end doubling, i.e. stripping a relation at an end. F-028 (FINDINGS.md around 1688) says "against an end the deletion is one mutation ... E-027 finds (L, ..., 2) -> (L) by a single left mutation at the sink, for every L it tried". That is close to the K>=3 claim. It is stated for removing a relation, whereas this submission deletes a free vertex. Related, not the same.
- H-020: placement of a core relative to head and tail.

Nothing found for "image class is a function of source class under free-end deletion" or for a K threshold. The K>=3 statement is therefore new as stated. The author's suspicion that it follows from F-028 and H-020 is plausible and not excluded. The negative result (a), (b) is also not recorded anywhere I grepped.

## Evidenced?

Mostly yes: tables give n, counts and per-class failures, and the commands reproduce.

Missing:
- Which classes are the 3, 5 and 10 K>=3 classes, and their sizes. A list would let a reader check that the K>=3 classes are not simply the rigid ones from (a).
- The treatment of '?' images at n = 10 (above).
- The "both deletions always admissible" claim for the pair statistic is stated, not shown. `removeVertex` returns None for length < 3, and I did not check that no pair hits it. n = 9 is fine by length.
- No test accompanies `maverick_classes.py`. The labelling is the base of every number (Coxeter key plus a free/edge/double/mirror orbit split of cospectral keys). It is checked only by "unresolved" counts. If the orbit moves are wrong, the class labels are wrong. At minimum, a check that the orbit labels agree with F-045 or with a mutation search on one cospectral key would back the labelling.

## Required for acceptance

1. State whether the n = 10 K>=3 counts include the 4 ends whose image label is '?', and recompute K>=3 at n = 10 with those ends excluded, or with the 4 ends reported separately.
2. Reword (c) to "all resolved classes having a free end of length >= 3", and move the gap-deletion 42/42, 184/184 sentence to "Does not claim".
3. List the K>=3 classes (key, size) at n = 8, 9, 10, so the reader can see that they are not just the rigid classes.
4. Say how `maverick_classes.py` was validated beyond the unresolved count, or add one check against an independent route (F-045 or mutation search) on a cospectral key.
