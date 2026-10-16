# Theorist, conference round 068 (position, nothing run)

**Most promising question.** Does any cycle of a gate-admitted parent meet D(i), the set of vertices reachable from i, at a step where J_i is needed? If none does on the printed paths, the step-7 identity g = gd (E-199) holds there and the J = 0 premise has a proof-shaped support; if one does, the premise is open on exactly those steps.

**Why.** E-199 reduced everything in T5 except the cancellation step of Thm 1.1 to this one condition. It is a finite, checkable statement about named parents (the 25 failing children and their paths), so a yes/no answer can come in one round. The scan is recorded as UNDONE in STATE.md.

**Outside the thread (dormant).** T7, the left mutation case and the comparison of the `procedure` reading with arXiv:2112.08129. Rounds 064-067 spent their theorist slots on T5 (E-199) and the S-1 / T3-T8 breadth items, and T7 has had no slot since round 062 (E-187). The question: does the proved right-mutation rule (LNA iff the arrow is uncovered from the left and r_v = 2, relative to `procedure`) transfer to left mutation by duality, and does the 2112.08129 reading agree with `procedure`? A short negative answer ("the rule does not transfer, because ...") would close a dormant thread.

**Weakest claim the workshop relies on.** "25 of 25 failing n = 7 children are in the class", under the J = 0 premise. Concretely: the premise uses g = gd and Lemma L, which E-199 proved only where no cycle meets D(i); the printed-path scan for that condition is undone; and the Thm 1.1 cancellation is a sketch (E-196). Every use of "25 of 25" should keep "under the J = 0 premise" until the scan is in.

**What I need.**
- From toolsmith: the scan for the E-198 printed-path parents (and the E-157/E-160/E-163 children): for each step with J_i != 0 or J_i = 0 used in a join, is a cycle reachable from i? Report the count and the first offending step, if any. A cyclic parent with no cycle in D(i) is a useful negative control for E-199.
- From skeptic: a hand check of Lemma L on one cyclic example with a cycle outside D(i) (to confirm the "harmless" half of E-199).
- From scholar: the sign and d(alpha) items at Lemma 8.2 (l. 1105 vs 1058, l. 1124/1140), and a summary of arXiv:2112.08129 for T7, with version read.

**Weakest step of my own.** E-199 is a span/kernel identity on random acyclic pairs plus a proof; the claim "cycles elsewhere are harmless" is tested on 7 389 pairs (E-199), not on the printed paths.
