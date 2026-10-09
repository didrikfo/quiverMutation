# Scholar, round 052 (conference)

Read: notebook (after 046), STATE, DIGEST rounds 048-051, STEERING. Nothing run; arXiv still blocked by the proxy, so every literature statement below is from memory and marked UNVERIFIED.

## Most promising question
Does any gate-admitted J_i != 0 step that keeps the LNA key join an LNA by a path where every step is tilting (Hom(T,T[-1]) = 0 for the whole T, not only at the mutated vertex), with a positive control at the same depth? A hit would promote `tiltingPlus` and settle H-015 for the 25 children; a miss is only bounded, but it tells us which of the two tests carries the weight. This is the one answer that changes the library, so it is the one to spend a round on.

## One question outside the thread
T3/T8 (dormant, scholar and toolsmith). The 8 of 10 odd-n key-coarser cores and all even-n cores are two mirror-closed orbits P, Q sharing one Coxeter key (E-077, E-080). A shared key means equal Coxeter polynomial, so the key cannot separate them. Question: is there a finer derived invariant (Cartan matrix up to congruence, or the Grothendieck-group data with the Euler form) that separates P from Q, and does a 2-orbit split ever occur across a single class? One n = 10 check with a control would say whether the parity name is a class fact or only an orbit fact.

## Weakest claim the workshop relies on
That a `tiltingPlus` pass at each step certifies a tilting step, and so a derived equivalence (the J = 0 premise). The local check is AI 2.32(b) at the mutated vertex (UNVERIFIED, from memory). The E-159 Hom test rests on the assumption that the printed paths generate, so it checks Cartan-level congruence, not the derived class. Neither has been checked on a path where the whole-T condition and the vertex-level condition could disagree.

## What I need
- **skeptic:** Hom(T,T[m]) replay for the E-158 paths with generation checked (or stated as an assumption in the record), and the whole-T test on one E-155 step.
- **toolsmith:** the group-A witness path (05040330 -> 33460000) by a relabelling-aware search, and a control for the tilting-only meet at depth 6.
- **theorist:** a statement of which finer invariant could separate P from Q (T3/T8), with a power check on n = 10.
- **maverick:** the S-1 n = 15 K0 = 5 sizing, unchanged; no new request.
- **chair:** this is a conference round; the ledger item for H-015 should keep the status OPEN and name the vertex-level versus whole-T distinction as the open point.
