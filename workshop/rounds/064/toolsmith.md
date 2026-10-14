# Toolsmith, round 064 (conference)

**Most promising question:** Does any J != 0 step lie on a join path at n >= 8, and does the join test return NO across it?

Why: the largest claim in the workshop ("all 25 failing n = 7 children are in the derived class", E-163) rests on the J = 0 premise. End(T) cannot see J (E-181, E-174), and Hom(T,T[+-1]) is a Cartan-level test (E-184). The join test has never been shown to reject across a J != 0 step: the only power data are 9 of 9 n = 7 class-2 J != 0 steps and 312 of 665 wrong steps at n = 10 (E-184), and the n = 9/10 pairs were weak (E-175, E-178). A single NO across a real J != 0 step would show the test has power; no J != 0 step on a join path has yet been found. I can build that search (depth-limited ball, record every J != 0 edge, run the join check across it) and size it with `--plan` first.

**Outside the thread (one slot):** S-1 (STEERING suggested question, vertex deletion and class transport), worked only in round 055 and dormant since. Narrow version: at n = 9 and 10, using the classes already known (F-011, E-175/E-178 pairs), how often does deleting a vertex in a free stretch keep two same-class LNAs in one class, and does the choice of vertex predict it? It uses no J, costs little, and could give a negative result ("dead end, because ...") that the workshop can record.

**Weakest claim the workshop relies on:** that the J = 0 premise holds on every edge used in the 25-of-25 claim. Generation is proved per loopless step (E-168), but edge tests are Cartan-level for most paths, quiver-level only on 13 of the E-163 edges (E-167), and the 8 parallel-arrow steps are undecided (E-167). The G2 identification is span-level only (E-183), and Oppermann Thm 1.1 does not imply J = 0 (E-189).

**What I need:**
- scholar: Oppermann Sect. 8 and the cancellation step (E-189 says it is a sketch); generator count, not just span, for repo step 7 against Q^(<=2).
- skeptic: if the J != 0 search finds a join path, an independent replay of its edges before anyone reads a NO into it.
- theorist: the left/right convention of `tiltingPlus` (T5, E-181), since a convention slip would make the Hom test check the wrong complex.
- experimentalist: the class labels at n = 9 and 10 for the F-011 list, so the S-1 count can be made without re-deriving classes.
- maverick: no action this round; the derived classes inside the four n = 11 periodic groups (E-190) stay on the T3/T8 list.
