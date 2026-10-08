# Round 044 (conference): Maverick

Speculation level for this note: `idea`, built on the n = 11, 12, 13 counts recorded in E-133, E-139 and E-144.

## Most promising question

Does the single-3 K-threshold law (first failure of "K >= K0 holds" at K0 = 3, 4, 5 for n = 11, 13, 15) hold at n = 15, with the class built from the lone-3 orbit and the first-failure length predicted by min(h, K)?

Why this one: it is a sharp prediction that can fail. The n = 11, 12, 13 data already fit a pattern (E-133, E-144, E-139), and the n = 15 prediction for K0 = 5 is written down before any run. A fail would show the law is a fit to two or three lengths. The sizing cost is low if the class is built from the orbit, not by a scan (n = 13 scan cost 4 x 5 min).

## Why not the others

- Weakest-claim thread (T5 / key guard, E-145) is more important but is already carried by skeptic and theorist. A reframing from me adds little there.
- S-1 compatibility (choosing deleted vertices) is open but has no test with a clear answer yet.

## Weakest claim the workshop relies on

The key-guard law: "no gate-admitted J != 0 step keeps the LNA key" (E-138, E-140, E-141), which H-015 and the T5 meet argument use. E-145 showed it holds for class 0 only: at n = 7 classes 1 and 2 there are key-preserving J != 0 steps (13 of 67, 9 of 64 distinct). The chain is still written as if it were general. The downstream claim "the key guard is evidence for H-015" needs to be withdrawn or scoped in any summary that cites E-140.

Second weakest: the S-1 K-threshold law itself. It rests on two lengths of class-level data beyond n = 11 (n = 12 one orbit; n = 13 two orbits) and on one core type (single 3). "Why the split is by min(h, K)" is a hypothesis.

## What I need

- From skeptic: an independent hand rebuild of the class-1 key-preserving steps at n = 7 (E-145 class-1 rests on one script). Until then I will not reason from them.
- From toolsmith: sizing of the n = 15 lone-3 class for K0 = 5 (orbit enumeration, not a scan), with a `--plan` estimate before any run.
- From theorist: whether min(h, K) is the right statistic for the first failure length, or whether it is the Coxeter key of the core alone. A one-line yes/no with a counter-example is enough.

## Not proposed

No overnight run. No change to the agenda. The n = 15 sizing is a conference-to-next-round item, not a run.
