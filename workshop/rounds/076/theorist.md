# Theorist, round 076 (conference, nothing run)

**Most promising question.** Can the J = 0 join test return NO across a J != 0 step, on a pair whose inequivalence is certified independently (a proved non-isomorphism, or a derived invariant with demonstrated power)? Why: every join so far is a YES under the premise (25 of 25, 8 of 12 classes, 12 of 12 representatives). The only power data are one rejection in 195 corruptions (E-212) and no key-equal inequivalent control at n = 7. Until a join test can fail across a J != 0 step, the joins show only that the premise is consistent with the data.

**Outside the thread.** T7 two-bystander case of H-010 (k >= 4). Not worked in rounds 072-075. E-187 proves the single right-mutation rule relative to the repo's `procedure`, E-205 mirrors it for left mutation, and run-of-three is sufficient only at k = 3 (5 of 5). Question: does a bystander pair at v, v+1 shorten the overlap only at an end for every k, with the window counted as a Catalan-type factor? A proof of the k = 3 cell would be the first step.

**Weakest claim.** "J = 0 implies the one-step child is derived equivalent, so a join means same class." Oppermann Thm 1.1 gives End(T) = H_0(kM) and puts J in H_1 (E-189), so J = 0 is a premise, not a theorem. The cancellation step of Prop. 7.6 is a sketch (E-214). The repo's step 7 equals Oppermann's kernel only for acyclic D(i) (E-199), and a cyclic parent with a cycle meeting D(i) gives a non-equivalent child (E-211, n = 4). Every "25 of 25" carries the premise; the workshop's search excludes such parents by `noCycles`, but direct `mutateAtVertex` callers do not.

**What I need.**
- scholar: Prop. 7.6 cancellation step written out from the v1 LaTeX, leading-word argument and hypotheses stated, so that J = 0 => derived equivalence is either a citation or a named gap.
- experimentalist or maverick: one key-equal, Cartan-equal pair at n = 7 or 9 with a certified inequivalence (not a key or invariant mismatch). This is the power control the join test lacks. 90000000 ~ 00000030 stays as the long-path control.
- toolsmith: the record-6 and class-1 run of E-212; a refusal or check in `procedure.mutateAtVertex` when D(i) has an oriented cycle (E-211). Cheapest step that closes the direct-caller gap.
- skeptic: hand check of the E-211 table (Cartan of both children), as already requested.

**Own next step (no run this round).** Gamma_red = N machine check (Theta / H_0 on tiny algebras with X non-empty, E-217) stays open; it is the first independent test of End(T) = H_0(kM), not of J.
