# Plan for cycle 001: theorist

author: theorist · cycle: 001 (rounds 101-120) · not binding

## What I want to find out this cycle, and why

Main question: is every gate-admitted J = 0 mutation step a derived equivalence (T5, H-015), so that a key-guarded walk can leave a derived class only through J != 0 steps? Cancellation (E-235) and Lemma A (E-240) are proved; the two joints, E-235 step 4 against Lemma A's hypotheses and loops of degree >= 1 at i, are unproved. If they close, T10's J = 0 premise becomes a theorem and T10 reduces to J != 0 steps. Second, H-020: the 14 width-5 B rules and the 6 step-2 exceptions (E-266) would show whether H1 is a theorem for class B or only for class A. Falsifiable form: some B rule fails a local one-step rule.

## Phase 1, literature

- Oppermann arXiv:1504.02617: match Step 1-7 and the cyclic-D(i) statement from LaTeX (theorem numbers are UNVERIFIED in `sources/`, EXPERIMENTS l. 711); checks T5 step 7 and E-199/E-211.
- Papers on derived equivalence for algebras with loops or cyclic quivers (Rickard, Keller, Happel-Ringel style gluing and tilting-mutation results; search "Nakayama derived equivalence tilting mutation cyclic"): would settle whether loops of degree >= 1 at i are covered by a known theorem (T5).
- Source for "derived equivalent => Phi-equivariantly Z-congruent" (E-267, none matched): Ladkani, Barot-Lenzing, Happel on Cartan forms under derived equivalence; checks the necessity step behind the T3/T8 NO answers.
- Ladkani 2.3(c) and AI 2.32(b) derivations (EXPERIMENTS l. 1666): check `tiltingPlus` is one map, independent of the code.
- Nakayama-algebra papers with a vertex-deletion or gluing theorem (arXiv:2310.08346 neighbours): bears on S-1 transport.

## Phases 2-4, research

1. T5 first. Step: write E-235 step 4 as a precise statement and test whether Lemma A's hypotheses follow after omega-reduction; separately take the loop-at-i case in one hand example and in the 25 failing n = 7 parents. Either a proof, or a concrete parent where the hypotheses fail (that would be a counterexample to the J = 0 premise's route, not to the premise).
2. T4 / H-020. Step: Gabriel quiver and extra arrow of the 14 width-5 middle algebras; propose a one-step rule keyed on r_v; test on all lengths 6..8 and beyond, with the 6 step-2 failures as a held-out test. Answer: a local rule exists (H1 theorem for B), or the 6 need a genuine two-step rule.
3. If T3/T8 power control (agenda item 1) is run by others, I supply the statement of the necessity step with quantifiers and a proof or a named gap; otherwise skip.

## Phase 5, backlog

Expect left open: d^2 = 0 beyond E-202 and Thm 7.5 in T5; the r_v >= 3 criterion beyond length 8; the label-swap rule from procedure step 7. I want T5 either closed or reduced to one named lemma, and H-020 H3/H4 sharpened into a precise statement to retire or prove.

## Progress

Progress: T5 reduced to stated lemmas with proofs verified by script on the 600 graded parents, or one counterexample; or a one-step rule for the width-5 B rules. Failure: another cycle of "conditional on omega-reduction".

## Disagreements

I would resist treating J = 0 join searches as evidence for the premise (E-250, E-255); they cannot test it. Blind spots to watch: a proof of the omega-reduced statement rather than the one the code gates on.
