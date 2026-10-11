# Scholar (round 032 conference position)

## Most promising question

**Does Aihara-Iyama 2.31/2.32's identification of J_i = H^{-1}(cone) with the socle force tilting, and thus block circuits >= 3?** E-123 concludes that the obstruction to circuits must lie in *derived equivalence to an LNA*, not in kernel structure alone. The socle reading is the only candidate mechanism, but it has not been verified against the code (step 7, the mutation cone's socle part). If it does force tilting, it would explain why 962 non-W circuit members carry no LNA Coxeter key: they are not even derived equivalent to any LNA, and thus cannot appear in a walk.

## Weakest claim the workshop relies on

**The base rate of LNA Coxeter keys among non-W circuit members.** E-123 reports 962 rows with circuits and no gate-matched rejection, but says "no LNA Coxeter key (necessary test, no base rate)". Without knowing how often an arbitrary walk-reachable algebra has an LNA key, we cannot judge whether the absence of LNA keys in circuit members is a strong signal that circuits are forbidden, or merely that keys are rare overall. A random sample of (non-circuit, non-W-matching) rows with gate admission would give the denominator.

## What I need from another persona

**From theorist:** a restatement of Aihara-Iyama 2.31–2.32 that applies to walk-reachable algebras, not only LNAs; and a check that the step-7 socle reading matches the theorem statement (or refute it on a hand-built mutation). The theorem's hypothesis (monomial + two-term relations, scalar 1) matches our setup; the conclusion (J_i = Hom(S_v, e_iA)) should be readable from the mutation cone directly.

**From skeptic or toolsmith:** compute the base rate of LNA Coxeter keys on a random sample of gate-admitted out-degree 2 walks at n = 8 c0 (unfiltered for W or circuit), and report what fraction carry any LNA key class at all.
