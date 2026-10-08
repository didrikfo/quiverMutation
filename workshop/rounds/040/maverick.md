# Maverick, round 040 (conference)

**Most promising question:** Does the deletion threshold for lone relations depend only on the Coxeter key of the lone relation's core, so that the length at which "K >= K0 holds" first fails (K0 = 3, 4, 5 at n = 11, 13, 15) can be read off without computing classes?

**Why:** If yes, S-1 becomes a checkable law with a prediction per length, and H-020's {h, K} dependence is the mechanism to look for. If no, S-1 is a list of two data points. Speculation level: `tested on small cases`.

**Weakest claim the workshop relies on:** the S-1 K-threshold law (K >= 3 fails first at n = 11, K >= 4 at n = 13, K >= 5 at n = 15). It rests on two lengths, n = 11 (E-118, E-133) and n = 13 (E-139). The n = 12 check and the n = 15 prediction are both unrun. The n = 13 "class" failure is also partly trivial: the witness is the lone 3 whose two ends have different keys, so it fails by construction under a class-level test. The n = 11 numbers have already moved once (E-125's 82 ends became 94 after `mirrorRow`). I would not write K >= 3 into `research/FINDINGS.md` on the present data.

**What I need:**
- experimentalist: the n = 12 run with `--plan` first. Report whether K = 3 fails at (3,5) and K >= 4 holds, with the key of the lone 3 at each end.
- theorist: whether H-020's {h, K} dependence gives a formula for the threshold, or a reason it cannot, before n = 15 is run.
- toolsmith: the `maverick_single.py` label-block fix (crashes per E-133), and a positive control for the "orbit = key class" test at n = 13.
- skeptic: a null. Does a simple length statistic of the core (arrow count, free-stretch length) predict the threshold as well as the key? If it does, the key adds nothing and the law should be stated without it.

**Watch:** image comparisons are by key, not label. "Differs" is sound; "same" is not. An orbit within a key class is only a lower bound for the derived class.

**Not for this round:** lone a-relations for a = 4, 7 (would test whether the law depends on a); the n = 15 run (overnight, after sizing).
