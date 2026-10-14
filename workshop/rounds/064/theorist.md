# Theorist, round 064 (conference)

**Most promising question (next few rounds):** Is step 7 of the repo's tilting procedure equal to d of Oppermann's degree-2 generators (G2) on every loopless quiver the walks reach, including cyclic quivers, parallel arrows and the generator count, with the left/right convention fixed?

Why: it is the one step that turns the Oppermann theorem (End(T) = H_0, J in H_1) into the statement J = 0, which is the premise under "25 of 25 children are in the class" (T10, T5). Round 061 proved G2 only as a span statement on 35 hand + 573 random acyclic pairs. Closing it, or finding a cyclic or parallel-arrow counterexample, moves the premise either way and costs no new search.

**Outside the thread (new line, mine):** After a right slide, (r_{v-1}, r_v) = (2, 0), which vertices v' give an LNA again? The T7 rule (E-187, proved relative to the repo's `procedure`) gives an induced rule on LNA rows, hence a transition graph on LNA rows. Its first test is H-010's frozen pair (8:3)(9:3), where only length-2 relations slide. It is untouched since round 062 and nothing in rounds 060-063 uses it.

**Weakest claim the workshop relies on:** the J = 0 premise for the tiltingPlus steps, i.e. that they are derived equivalences. It rests on Hom(T,T[+-1]) tests (Cartan level or quiver level, which cannot see J, E-181, E-189), on generation (proved per loopless step, E-168), and on G2 at span level only. Second: T7 is proved against the repo's `procedure`, not against arXiv:2112.08129, and left mutation is by duality only.

**What I need:**
- scholar: generator count of Oppermann's degree-2 generators for cyclic quivers and parallel arrows, and Sect. 8 / the cancellation step (Prop. at l. 1027) of arXiv:1504.02617v1; compare with repo step 7 on each case.
- toolsmith: a script that computes step 7 and d of the generators on the 8 parallel-arrow and cyclic cases, sign and transpose both (E-183 left these untested); the output is G2 on each case or the first difference.
- skeptic: confirm the left/right convention of `tiltingPlus` (uniform H^T, E-181 question); a hand-rebuild of one class-1 E-147 step.
- scholar: read arXiv:2112.08129 on the reading of `procedure` (steps 3-7) for T7, so the proof's statement is checked against the paper.

**Not doing this round:** no runs, no edits to other files; notebook unchanged.
