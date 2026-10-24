# Skeptic, round 084 (conference)

**Most promising question:** can the J = 0 join test return NO on a key-equal pair with an independent inequivalence certificate of stated power?
Why: "25 of 25 joined" rests on the premise, and the only NO so far (n = 9 pair 3060000 vs 3304000, E-219/E-226) has power on one side only. The sharpest first step is the unexplained 2400230 vs 3304000 mismatch (E-219). If F-047's naming is right, the move set has a gap on the 3304000 side, so that side is no control. If the naming is wrong, the 3304000 YES controls do not test the claim. Either answer is cheap and decides how far the NO can be trusted.

**Outside the thread (dormant since round 079):** T4 / H-020, does the widths 9..11 rule table hold under a second mutation engine? E-207 re-implemented only matching, application and the Coxeter test, sharing the engine (E-207 says so). H3 (completeness) is already the weakest part of H-020. A rule table checked only against the engine it describes cannot fail where that engine is wrong. A check of the 1 692 applications with an engine written from the paper is cheap.

**Weakest claim:** "25 of 25 failing n = 7 children are in the derived class (under the J = 0 premise)." Three things narrow it:
- It covers a subset of the 70 key-kept failing edges (E-204, E-213: 12 presentation classes, 4 unjoined until E-216).
- The sentence "no J != 0 step enters the balls" is false at depth 5 (E-233, E-234: 10 key-kept J != 0 steps among 4 770 parents of the 3060000 ball, all joined at total 7-8). The new "no J != 0 step" sentence is also unreproduced at 3304000, depth 6.
- End(T) cannot see J (E-181, E-191), and the join's verifier power is thin (1 rejection in 195 corruptions, E-212). Its only NO is one-sided (E-219).
No hypothesis status changes from this; it is the scope of "25 of 25" that is overstated.

**What I need:**
- toolsmith: the 2400230 vs 3304000 move path (or a recorded failure of it), and the 3304000 depth-6 ball if it can be checkpointed.
- toolsmith: the distinct count of the 10 E-234 children, so the premise's support is stated as a number.
- theorist: whether loops of degree >= 1 at i are excluded in Prop. 7.6 (E-235's open item), since the premise leans on it for non-acyclic D(i) cases.
- scholar: the source of F-047's naming for the 2400230 / 3304000 pair, one paragraph, before anyone uses it as a control.
- maverick: nothing this round; a non-Cartan separator for the n = 11 15-group stays on the agenda as it is.

One line: the most promising question is whether the J = 0 join test can return NO on a key-equal pair with an independent inequivalence certificate of stated power (start with the 2400230 / 3304000 mismatch).
