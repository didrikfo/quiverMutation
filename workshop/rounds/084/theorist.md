# Theorist, round 084 (conference)

**Most promising question (next few rounds):** Does any J != 0 step ever leave the derived class, and can a join test with an independent inequivalence certificate return NO on such a pair?

Why: every "joined" result in T10 (25 of 25, the 4 unjoined classes, the 3060000 ball) is a YES under the J = 0 premise, and no test has yet been able to fail across a J != 0 step. End(T) cannot see J (E-181, E-174: J lies in H_1, End(T) = H_0), so the next test needs an invariant that sees H_1. Ask: which invariant of the silting object changes exactly when J changes? That is a mechanism question, not another count.

**One promising question outside the thread that got most of the last 4 rounds (dormant line):** Why do floating H-020 rules not depend on length (H1), and does the Catalan factorisation of E-187 give the reason?

Why: H1 rests on pass counts (8 248 applications, E-188; 1 692, E-225) and nobody has asked for a mechanism. T4 has had only breadth slots since 071, and its H3 (completeness) is the weakest of its six. Data says what; a length-free argument would say always. Not run; a sketch only.

**Weakest claim the workshop relies on:** "25 of 25 failing n = 7 children (and the 4 unjoined classes) are in the derived class." It rests on J = 0 => equivalence, and my own part of that chain is weaker than the notebook's tone suggests:
- Cancellation kM -> N (E-235) is proved by hand only modulo the cited Anick Lemma A and a printed-degree fact; about 16k letters machine-checked, not a full proof.
- Loops of degree >= 1 at i are not covered; step 4's d_1 is not monomial (6 968 of 10 206 loop letters differ). Open.
- Repo steps 1-7 = End(T) is proved only for acyclic D(i) (E-199); E-211 is a cyclic counterexample. The search only reaches acyclic nodes (E-223), so the claim survives for the search, but not as a general statement.

**What I need:**
- scholar: Lemma A's theorem number and statement matched in the v1 source (or UNVERIFIED); whether the journal version of 1504.02617 changes Prop 7.6.
- skeptic: a hand check of the omega filtration in E-235 (the step the proof leans on: acyclicity of the Morse matching, cited).
- toolsmith: one parent with a loop of degree >= 1 at i (n <= 5) and a kM -> N check with signs; I cannot run it from here.
- maverick: nothing this round; the NO test (question above) needs an invariant before any run.

No run proposed this round. No change to any H-, F-, R- status is proposed by me.
