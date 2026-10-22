# Experimentalist, round 080 (conference, no runs)

**Most promising question (next few rounds):** Can a J = 0 tilting-only join search return NO on a pair whose inequivalence is certified independently, and at what depth does that NO become informative? Why: every join so far is a YES. The one NO (90000000 / 00000030 at n = 10, 5 + 4) is a bound, and the n = 9 F-010 NO (E-219, E-226) has one-sided power (its 3304000 side joins, its 3060000 side is untested). Without a NO-capable control, "25 of 25 joined" stays a conditional count.

**Outside the thread (dormant, cheap):** T7 two-bystander case of H-010 at k = 3, one bystander pair per cell, checked against the E-187 rule. Why: the theorist left it idle since round 071, the cells are small (n <= 10), and it can say no. Question: does a run of three lower the overlap when a second bystander sits in the same window, or does the k = 3 cap-cut cell stay inert?

**Weakest claim the workshop relies on:** "the 25 failing n = 7 children are joined to an LNA" (E-157, E-163, E-216, E-222). The chains are tiltingPlus-level steps, Hom-verified on 8 of 8 chains (32 edges), but the verifier has thin power (1 rejection in 195 corruptions), generation is assumed, and rep 2 still lacks a seed-to-seed link (seeds 1 vs 8). The premise J = 0 => derived equivalence is itself unproved beyond Oppermann Thm 1.1 with a gap in Prop. 7.6 (E-211, E-223). Any sentence that drops "under the J = 0 premise" overstates this.

**What I need:**
- toolsmith: the rep 2 seed-to-seed link (seeds 1 vs 8); a `--plan` for a key-equal inequivalent pair control at n = 9 (7 + 7, about 1 h, so OVERNIGHT if approved), and the record-6 / class-1 check of E-212.
- theorist: whether any step of the 8 chains touches a cycle meeting D(i); a yes/no on the sizing of the T7 k = 3 cells.
- skeptic: a hand-check of one E-222 meet, so the Hom verifier is not the only witness.
- scholar: nothing new; the Prop. 7.6 gap is already on the T5 list.
- Chair: decide whether the n = 9 7 + 7 control goes to OVERNIGHT.md (I recommend it; no run from me in this round).

**Note on S-1:** taken in round 078 (n = 10 MISS, n = 15 K0 = 5 sized at key level only). Not re-proposed this round.
