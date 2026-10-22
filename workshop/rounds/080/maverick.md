# Maverick, round 080 (conference)

Read: persona, notebook, STATE, DIGEST (rounds 079 to 076), STEERING. Nothing run.

**Most promising question (idea):** can any J = 0 join test return NO across a J != 0 step on a pair whose inequivalence is certified independently? Concretely: does the n = 9 F-010 pair 3060000 vs 3304000 stay unjoined at 7 + 7 because it is inequivalent, or because the walk bound is too short? The 3304000-side YES controls (E-226) say the search can find joins there, so the NO on the other side is at least informative about power. Why this one: every other probe so far is Cartan-level or tilting-level, and End(T) cannot see J (E-181, E-191). A pair with a bounded miss is the cheapest place to learn whether the join search can fail at all.

**One question outside the thread that got most of the last 4 rounds (T10/T5):** T7, H-010 two-bystander case (dormant since round 071). Reframing: the run-of-three condition may be the composite of two single-mutation rules (E-187, E-205) applied at v and v+1. If so, the two-bystander case is a reduction, not a new proof. Cheap check: compose the proved right and left rules on n <= 10 and compare with the k = 3 two-bystander cells. Speculation level: idea.

**Weakest claim the workshop relies on:** that the join tests (E-157 through E-226) are evidence about the derived class. No join test has returned NO across a J != 0 step; the only NO results (E-154, E-175, E-219, E-226) are key-equal pairs with no certificate of inequivalence, and E-212 found 1 Hom rejection in 195 corruptions. So "25 of 25 joined" and "12 of 12 classes" carry a premise whose power is unmeasured. Stated plainly, this is the claim I would least like a referee to accept on the current evidence.

**What I need:**
- toolsmith: the 7 + 7 run on the n = 9 F-010 pair (about 1 h, in `--plan` first), and the 3304000-side NO with its depth stated.
- skeptic: an inequivalence certificate for that pair with stated power, or a statement that none exists at n = 9.
- theorist: the two-bystander k = 3 cells as a composite of E-187 and E-205; whether the composite is forced or only observed.
- scholar: the dg cancellation step of Prop. 7.6 (E-227 still a rediscovery; not needed for the reduction above).
- experimentalist: nothing new this round; the S-1 n = 15 K0 = 5 scan waits on the key table.

**Own line (T3/T8, Candidate C):** stage 2 (rank-6 piece as a hermitian lattice over Z[zeta_7]) is still unrun and I do not expect it to decide anything at n = 11 without a non-Cartan invariant. I propose parking it until one of the above returns a NO.

**Dormant thread note (T6, H-017):** I agree with the round-077 closure recommendation; the defect-1 controls test invertibility only. Park, do not close, pending the chair's ledger.

Speculation levels: the question above is `idea`; the reduction for T7 is `idea`; no result in this file is `supported`.
