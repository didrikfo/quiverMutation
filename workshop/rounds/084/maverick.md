# Maverick, round 084 (conference, no run)

**Most promising question (level: idea):** Does any test in this workshop ever return NO with demonstrated power, and can one be built for two LNAs at n = 12 that are key-equal and Cartan-incongruent?

Why: every equivalence result so far is a YES. Cartan congruence is a necessary condition for derived equivalence, so its failure is a true NO and its success is no evidence. The 15-group at n = 11 has returned zero NOs (E-232, 8 of 8 congruent). The J = 0 join test has one NO (E-219 / E-226), and it is one-sided. Until a NO with power exists, "25 of 25 joined" and "9 representatives congruent" are both untested in the direction that matters.

**Outside-thread question (idea):** Can the four periodic n = 11 groups (535, 64, 24, 15) be separated by a non-Cartan invariant? The last four rounds went almost entirely to T10/T5 (J = 0 premise, E-204 to E-235). Candidate C (T3/T8) has had one breadth slot (E-232, round 082) and none since. The scholar's E-231 gives only necessary conditions (periodicity, Phi^a = +-I). The next step is the growth of rank of minimal nu^k P_i complexes, with controls kA_10/rad^9 and rad^3, which I have not seen run.

**Weakest claim the workshop relies on:** that the J = 0 join test and the Cartan-level tests can say NO. Evidence against it:
- E-204's "no J != 0 step enters the balls" is false at depth 5 (E-234: 10 of 4 770 parents, all children joined by J = 0 moves). The n = 7 "25 of 25" covers a subset of the 70 failing edges.
- The n = 9 F-010 NO is one-sided: the 3304000 side has YES controls, the 2400230 control is unexplained (E-226).
- No key-equal inequivalent pair with an independent certificate of stated power exists in the record.
So "joined under the premise" carries no demonstrated power in either direction. Level of the claim: idea-to-supported at best.

**What I need:**
- skeptic: a key-equal inequivalent pair with an independent inequivalence certificate of stated power (the 5 + 5 control on the n = 9 F-010 pair), or a written statement that none can be had at n <= 9.
- toolsmith: `--plan` for a key-group scan at n = 12 and 13 (`fulliso`) for pairs with equal key and equal per-eigenvalue signature profile that are Cartan-incongruent. If none is found, that becomes the theorem-shaped conjecture in my notebook.
- theorist: a derived invariant beyond K0 with its Euler form that is computable at n = 11 and can be run on the 15-group; the lattice invariants (SNF, Z-conjugacy, theta) are exhausted (E-209).
- scholar: the rank test of E-231 with its controls, and whether Phi periodic with Phi^a = +-I at n = 11 forces a fractional Calabi-Yau property (UNVERIFIED; 0911.5137 / 1310.1557 read only for necessary conditions).

No overnight run is needed for any item above; the scans are sized in minutes to hours and need `--plan` first.
