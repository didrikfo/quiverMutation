# Experimentalist, round 056 (conference)

**Most promising question:** Can the J = 0 join test ever fail? Run `tiltingPlus` and the Hom(T,T[+-1]) test on a pair known to lie in different derived classes with equal key (the n = 10 P/Q pair, E-079/E-082, or a Phi-group pair from E-172), and on the 25 children with the control in place.

Why: every join so far (E-157, E-160, E-163, E-166, E-169) is a positive. E-169 records that the J test has not been shown able to fail (0 of 1202 gate-admitted edges). Without a demonstrated reject, "25 of 25 joined" passes by construction, and the control is cheap next to any new census.

**Outside the thread (dormant T3/T8):** Do the 4 classes of the one n = 10 Phi-group (E-172) get separated by a non-Cartan derived invariant with a power control? The group is small, so this fits in minutes. A null result, with the control, is a valid dead end.

**Weakest claim:** "25 of 25 failing children lie in the LNA class." It rests on the J = 0 premise, whose test has no demonstrated reject. The quiver-level End(T) check covers the 13 E-163 edges only; the 8 parallel-arrow steps are undecided (E-167). Hom is tested at Cartan level on 324 + 33 edges.

**What I need:**
- toolsmith: a decision on the 8 parallel-arrow failing steps, and a relabelling-aware `meetingPoints` so the E-169 group-A paths are found by the library, not only by my script.
- skeptic: an independent replay of the 5 E-169 paths (35 edges) with the Hom test.
- maverick: a non-equivalent control pair with equal key (P/Q at n = 10, E-079/E-082) for the power test.
- theorist: whether generation per loopless step (E-168) covers the parallel-arrow case, or whether those 8 steps need a separate argument.

**Notes:** The F-037 cross-check (notebook step 1) waits for the control. No run proposed this round.
