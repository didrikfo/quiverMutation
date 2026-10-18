# Experimentalist, round 072 (conference)

**Most promising question (next few rounds):** At n = 10 -> 9, does the D >= 3 free-stretch deletion give one image class for all 42 labels when the labels are exact (`classifyLength(10)`, not key+orbit)? Why: E-206's "10 of 10" uses key+orbit labels, with 176 LNAs unresolved and dropped, so it is not yet a class-level result. Exact labels either confirm it or find the first failure. Size with `--plan` first; the labelling is an overnight job (resume=True, about 3 runs).

**Question outside the thread that got most of the last 4 rounds (T10/T5, the J = 0 premise):** Does the n = 10 equivalent pair 90000000 ~ 00000030 (equivalent, not joined by the orbit-plus-mirror moves, E-190 / T3/T8) keep its equivalence under deletion? Test: apply the D >= 3 deletion to both members at n = 10 -> 9 and compare images. If the images are in one class, deletion transports this equivalence, a case the walk does not reach. If not, the pair is a control for the S-1 question. Cheap (seconds on cached data if the classes are cached; check `research/EXPERIMENTS.md` first). It is a dormant thread (T3/T8) crossed with S-1, not a new line.

**Weakest claim the workshop relies on:** that the n = 10 labels are classes. E-206 and the S-1 table rest on key+orbit labels; the mirror and the 176 dropped LNAs mean "one image class" is a statement about proxy labels. The D = 2 failures inside one exact class (E-206) are the only case where proxy and exact would disagree, so they are the first thing to recheck.

**What I need:**
- toolsmith: `classifyLength(10)` with `resume=True` sized with `--plan`, and a `--plan` estimate for the exact n = 10 labelling of the 42 keys, for the chair to put in `OVERNIGHT.md` if approved.
- maverick: confirm the 90000000 / 00000030 pair and its n = 9 members in the cached CSVs, so I do not rebuild them.
- theorist: whether D >= 3 has a reason (a vertex far from every relation support is a free stretch of the walk) that predicts where the rule fails at n = 11, so the n = 10 result can be a test and not only a count.

Not run: no commands in this round; nothing edited outside this file.
