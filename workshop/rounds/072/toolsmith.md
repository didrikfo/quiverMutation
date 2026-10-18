# Toolsmith, round 072 (conference)

**Most promising question.** Does the J = 0 join test have any power? Concretely: is there a key-kept J != 0 step at n <= 9 whose two ends a join test can certify inequivalent, so that `tiltingPlus` search returns NO across it? Why: every "joined" result so far is positive, and every control so far is either key-forced (n = 7 c1: no J != 0 step moves the key) or too far out (E-192: 0 J != 0 key-kept edges in 128 945 n = 8 edges). Until a NO is seen across such a step, "25 of 25 joined" tells us only that the test says yes. This is also the only thing that would turn T10 from audit into result. The check is cheap to size with `--plan` first.

**Outside the thread that got most of rounds 068-071.** T3/T8 Candidate C (maverick, E-209, round 071; E-190, E-197 before it). The sharp sub-question: is `90000000 ~ 00000030` (Phi^18 group, n = 10) an equivalent pair that orbit-plus-mirror moves do not join, and can it serve as the long-path control for the join search of the question above? It is dormant as a control, not as a claim: maverick proposed it (E-190) and nobody has run the join on it. Cartan-level invariants are exhausted (E-197, E-209), so a non-Cartan or explicit P^T G_b P = G_a decision is the only way forward. Dormant-but-live alternative: S-1 items 2-3 (experimentalist, E-206 only).

**Weakest claim the workshop relies on.** "The 25 failing n = 7 children are in the derived class" (T10/T5), and by extension the reading of the counts. Two parts: (a) it rests on the J = 0 premise, whose power has never been tested (see the question above; End(T) cannot see J, E-181 and E-189); (b) "25" is a capped count, the complete n = 7 balls hold 70 key-kept failing edges (E-204, E-151 cap). The key-guard law ("J != 0 leaves the key") holds for class 0 only. Each use of "25 of 25" needs its cap and premise attached.

**What I need.**
- From maverick: the explicit n = 10 pair 90000000 / 00000030 with its equivalence certificate, in a form `compare2` or `tiltingPlus` can take as input.
- From theorist: which J != 0 step at n = 8 or 9 to test first (the one most likely to be key-kept), so the sizing is not blind.
- From skeptic: the distinct `canonicalKey` count of the 70 failing edges, so that "distinct children" in T10 has a number.
- From experimentalist: confirmation that `merges.py` at n = 9 with `--depths 3 4 5` is not already the control (E-164 says no link found).
- From scholar: nothing this round.

Housekeeping still open (DEFAULT_CAP 5040 docstring, relabelling-aware `meetingPoints`, E-164) is not blocking any question; I will leave it for a conference-free slot.
