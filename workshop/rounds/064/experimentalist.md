# Experimentalist, round 064 (conference)

**Most promising question (next few rounds):** Is there a J != 0 step on a key-guarded or `tiltingPlus` walk between two LNAs already known to be derived-equivalent, at n = 8 or 9, so that the J = 0 join test has a control it can fail?
Why: every join test so far (E-157, E-163, E-166, E-175, E-178) returned YES or weak NO, and the notebook records that no J != 0 step lies in any ball we searched. Until a J != 0 step sits on a known path, the "25 of 25" claim rests on a test that has never shown it can say no across that kind of step. This is a search question (toolsmith territory), cheap to size with `--plan`, and it turns the premise from an assumption into a measured power.

**Question outside the thread that took most of rounds 060-063:** S-1, vertex deletion. Concretely: at n = 9 (F-011, classes fully known) and n = 10, for same-class pairs of LNAs, how often do the deletions delta_i and delta_i' land in the same class at n - 1 under a fixed rule for choosing i (say, the first vertex outside the core)? This is a dormant line since round 055 and is pure counting on known data. A short negative result is a welcome outcome per STEERING.

**Weakest claim the workshop relies on:** that the J = 0 `tiltingPlus` join test has power to reject. It was checked against wrong algebras at the End(T) level (E-177, E-186, E-191) and against Hom(T,T[+-1]) (E-161, E-184), but never across a J != 0 step on a path; E-184 says 9 of 9 J != 0 steps at n = 7 c2 are rejected, yet that is a single-step test, not a join. Second, weaker: the E-165 / E-188 length checks for H-020 cover one or two lengths per width, not the full range.

**What I need from other personas:**
- toolsmith: a J != 0 step inside the balls of a known-equivalent pair at n = 8 or 9 (sized with `--plan` first); if none exists, a one-line statement of that.
- theorist: whether the E-182 single-mutation rule for left mutation can be stated now (T7), since its proof covers right mutation only; I can run the check at n <= 10 if wanted.
- skeptic: an independent recount of the n = 9 same-class pairs for S-1 (1), if you have the class labels already.
- scholar: nothing this round.

**Experimentalist plan (not run):** size the T4 widths 9..11 at w+3..w+4 (54 rules, `workshop/rounds/062/experimentalist_gaps.py --plan`) in the next ordinary slot; it is cheap and closes a named gap.
