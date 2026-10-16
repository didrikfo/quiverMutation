# Toolsmith, round 068 (conference, position statement; nothing run)

**Most promising question (for the next few rounds):** Can the J = 0 join test return NO across a real J != 0 step, with an independent certificate that the two endpoints are inequivalent? Why: every join so far is positive, End(T) cannot see J (E-181, E-189), and a power control at n = 7 gives 9 of 9 J != 0 rejections only through the same Hom test. A search that can fail is the only thing that would move T10 off "assumed". Sizing (E-192) says it is multi-hour at n = 8, so the first step is a shallow n = 7 control with a certificate, not a deep search.

**Outside the thread (dormant, cheap):** T3/T8 orbit-plus-mirror. At n = 10 the two mirror-closed orbits P, Q of the even-n key-coarser cores share one key (E-079, E-082), and the "no shift by 1" question is unproved. Question: does any Z[Phi]-isometry or SNF(xC + C^T) invariant separate P from Q at n = 10 (a 2-orbit, one-group case)? It is dormant, not closed, and the n = 10 group is small enough to test in one command. Maverick's n = 11 SNF result (E-197) does not cover this case.

**Weakest claim the workshop relies on:** "All 25 failing n = 7 children are in the derived class." It rests on the J = 0 premise: the Hom(T,T[+-1]) test is Cartan-level, End(T) = H_0 and J lives in H_1 (E-181, E-189), and the cancellation step of Oppermann Thm 1.1 is a sketch (E-196). Edge tests cover 370 n = 7 path edges and 35 n = 10 edges; the 8 parallel-arrow steps are decided only by quiver-level End(T). Every "25 of 25" must keep "under the J = 0 premise".

**What I need from other personas:**
- theorist: the cancellation step of Thm 1.1 (or a statement of where it fails) so that J in H_1 has a computable test; a cycle meeting D(i), since cyclic children are untested.
- scholar: Lemma 8.2 sign (l. 1105 vs 1058) and the d(alpha) condition (E-196); whether Sect. 8 needs J = 0.
- skeptic: a hand check of one J != 0 step's Hom(T,T[-1]) != 0 (silting but not tilting) on a case the code did not generate.
- experimentalist: a control pair at n = 7 with a J != 0 step on the path and a known inequivalent endpoint (the "can say NO" power control); the 9 class-2 steps are the candidate set.
- maverick: the n = 10 Phi^18 pair 90000000 ~ 00000030 as the long-path control (E-190), if the power control above is sized.

**Toolsmith backlog, not blocking:** cyclic-parent count in the E-192 balls (only the printed paths are done, E-198); `canonicalKey` DEFAULT_CAP 5040 and the guard/`canonicalKey` docstring rewords (waiting for a toolsmith round, as the chair decided in round 052); `merges.py --witness` on a real link; `compare2` dimension assertion is in place (E-180, E-191).

**Suggested question:** none new. S-1 is taken by experimentalist (E-200); the T3/T8 item above is the one I would take up in a breadth slot.
