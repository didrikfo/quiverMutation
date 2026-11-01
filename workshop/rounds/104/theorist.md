# Theorist, conference 104 (phase 1, literature)

No new work. Formed from my notebook, my 101 and 103 submissions and reviews, E-265..E-280, STATE and the cycle plan.

**Most promising question (T5/H-023).** Can the "no loop of degree 0" hypothesis of Oppermann Thm 1.1/7.5/7.6 be carried to a Step-3/4 elimination for a loop of degree >= 1 on an explicit parent (two-stage: kill Step-4 letters on gr, then Step-3 letters on E^1)? Why: the phase showed the source has no loop hypothesis beyond degree 0 and that library steps never meet a degree >= 1 loop (E-279). So what is left of H-023 is one stated lemma of ours (Lemma A for loop letters), checkable on one parent with Ext^2(S_i,S_i) != 0 and an independent dim End(T). It is the cheapest way to turn "J = 0 premise" from an assumption into a theorem or a counterexample.

**Outside the thread.** T4/H-020 width-5 class-B rules: E-266 left 14 width-5 middle algebras and 6 step-2 failures (w6/w7) without a one-step rule. A rule for the width-5 first step, with the 6 held out, is small and decidable, and H3 (completeness) is open only there.

**Weakest claim the workshop relies on.** "25 of 25 joined" and the 70 n = 7 edges resting on the J = 0 premise: the premise is proved only for one step on acyclic D(i), under Construction 2.2 minimality (assumed in my loop argument) and not through H_0 chains (E-249). Second: the Cartan-level NOs rest on E-269's float bound plus a necessity step whose Nikulin number (1.13.3) is still UNVERIFIED and secondary.

**Closest to a finding.** E-279: a degree >= 1 loop at i in a Construction-2.2 model forces i onto an oriented cycle of Q0 (341 random parents, 0 violations, 1000+ loops, non-vacuous). Needed check: it is a proof sketch by induction on degree plus a random sweep with the Ext route only; one independent check is the graded-quiver Step 3/4 rule on a parent with a cyclic i (E-248 style), or a hand proof of the m = 1 base (d phi is a relation) with parallel arrows and non-monomial relations written out, refereed by the skeptic.

**What the literature changed for my thread.** Lemma A is not in Oppermann (E-273); the cancellation is Prop. 7.6 and needs only no degree-0 loop. AI Thm 2.31/2.32 and Keller (Rickard) have no quiver hypothesis, so a cyclic or loop case is not excluded by any theorem; our gap is ours. Anick Th. 3.1/3.2 and Forre Th. 2.6 are known only second-hand (Gartner, Efrat); Rickard's own number is UNVERIFIED. Ex. 3.2 (cyclic D(1)) does not reproduce with `allowCyclic=True` (KeyError, E-277).

**Phase 2 should not.** Re-prove Oppermann's Step 7 kernel or the cancellation (published, E-273); run more join searches (E-250, E-268, plan); sweep more random monomial cyclic parents through the library (the default refuses by construction of the gate, so counts show nothing, E-193); claim a theorem for Anick or Nikulin numbers without the primary; try to source the Catalan factorisation or n 2^(n-5) in the literature (searched, none found, E-280).

**What I need.** Scholar: Anick 1982 primary (Th. 3.1/3.2), Rickard 1989 theorem number, Nikulin 1.13.3 text. Toolsmith: a Hom/End check of one `allowCyclic=True` child with i on a cycle (it returns unproved vertex-simple paths in 90-101 of 300). Skeptic: referee the m = 1 base with parallel arrows and non-monomial relations, and the E-235 step 4 hypotheses against Lemma A (E-256 a).
