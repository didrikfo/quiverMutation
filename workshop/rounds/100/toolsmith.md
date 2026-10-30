# Round 100 (conference) · Toolsmith

Nothing run. Read: persona, notebook, STATE, DIGEST rounds 096-099, STEERING suggested questions.

**Most promising question (next few rounds):** Can the J = 0 join test return NO with demonstrated power on a key-equal pair? Concretely: for 3060000 vs 3304000 (n = 9, F-010 / E-057 certified inequivalent), does a 7 + 7 J = 0 search, with the recorded key-kept J != 0 edges added as moves, fail to meet while a YES control on the 3304000 side meets at the same total? Why: it is the only test that can show the premise failing where we know the answer. Every NO so far is one-sided (E-219, E-226) or a move-set MISS (E-255, E-268), and the workshop has named this gap in every conference since round 080. Cost is about 2 h by my estimate, not sized with `--plan`, not approved.

**Promising question outside the thread that got most of the last 4 rounds:** Dormant H-017 (Euler signature). Scholar found no source for "relations minus cords" (E-261), and Tits form = Euler form as matrices holds iff gl.dim <= 2 (derivation). Question: does the defect statement hold on the gl.dim <= 2 subclass alone? That could be checked on the 19 defect-1 members already computed (E-221, E-203), with no new search. This is a line the last four rounds did not touch (T6 was closed in round 080 and not reopened).

**Weakest claim the workshop relies on:** "25 of 25 failing n = 7 children are joined to an LNA", and so "in the class", under the J = 0 premise. The join chains are checked at tiltingPlus and Cartan level only (E-222, E-228). Hom+Cartan has no NO power: 0 of 26 key-kept corruptions rejected (E-228), and End(T) cannot see J (E-174, E-181). Generation of K^b(proj) is proved per loopless step (E-168), but the step-7 cancellation for loops of degree >= 1 at i is still open (E-235). The chair's own checklist already says this; I would add that the verifier's power is the thin part, not only the premise.

**What I need:**
- theorist: the E-235 step 4 check against the Lemma A hypotheses (E-240), so the cancellation claim has a status that is not "conditional".
- skeptic: an inequivalence certificate for one key-equal pair with a stated power number, to set the 7 + 7 NO against.
- experimentalist: a `--plan` line for 3060000 7 + 7 with checkpointing and its wall time, so the human can decide on approval.
- maverick: the one-line result of a power table for the `q_S` test on the n = 12 group of 140 (E-269 says not Z-congruent; say whether that is a NO with power or a separation at the Cartan level only).
- human: a yes or no on the 7 + 7 run once it is sized (STEERING, overnight menu).

**Suggested question (optional):** none. The S-1 line in STEERING is already taken by the experimentalist (E-263) and is not a toolsmith item.

Most promising question: can the J = 0 join test return NO with demonstrated power on the key-equal pair 3060000 vs 3304000 at n = 9?
