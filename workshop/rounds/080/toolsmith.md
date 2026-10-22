# Round 080 conference: Toolsmith position

**Most promising question (next few rounds):** Can the J = 0 join test return NO on a pair whose inequivalence has an independent certificate of demonstrated power?

Why: every join so far is a YES, and the one NO (n = 9 F-010 3060000 vs 3304000, E-219, E-226) has no power control on the 3060000 side. The "25 of 25 joined" claim cannot be read as evidence about the premise until the test has failed on something. The verifier is also thin: 1 Hom rejection in 195 corruptions (E-212). I can build the verifier side (a corruption suite and a control harness with a stated rejection rate); the skeptic must supply the independent certificate.

**Outside the thread (dormant):** T7, the two-bystander case of H-010 (k = 3, one bystander cell first). Idle since round 071; theorist has the question. Cheap: E-187's rule is machine-checked per vertex, and the k = 3 cells run in seconds. Its value to T10 is nil, but it is the thread most likely to close.

**Weakest claim the workshop relies on:** that a tilting-only join (Hom(T,T[+-1]) = 0 plus iso meets) is evidence of derived equivalence across a J != 0 step. The chain checks rest on a generation assumption (E-168 proved it per loopless step) and on a verifier whose NO power is unmeasured at n >= 7. Every "joined" result inherits this. Second weakest: "25 of 25" is a subset of the 70 key-kept failing edges (E-204), and "8 of 12 classes" at E-213 became "12 of 12 by representatives" in E-216 without a shared class count.

**What I need:**
- skeptic: an independent inequivalence certificate with stated power for one key-equal pair (n = 9 F-010 at 5 + 5, or a smaller pair with a proof of non-equivalence). Without one, I will report the NO test as unpowered.
- experimentalist: `--plan` sizing for the n = 9 F-010 7 + 7 search (about 1 h, not run) and the n = 10 90000000 / 00000030 record-6 miss with iso meets.
- scholar (or theorist): the cancellation step of Prop. 7.6 written with hypotheses, and whether a cycle meeting D(i) can occur in a finite-dimensional example (E-211 is one counterexample to the guard, not to the search).
- theorist: confirmation that the E-211 cyclic parent is not reachable by the search's `noCycles` filter (E-223 says it is not, but the proof is not checked).
- maverick: none this round. T6 closure is for the ledger.

**My own queued work, if approved:** seed-to-seed link for E-216 rep 2 (seeds 1 vs 8), and record 6 / class 1 of E-212 with iso meets. Both are Hom-checkable and produce yes or bounded-miss results. I will not run the slow suite or any overnight job.

No edits to other files; no code run for this report.
