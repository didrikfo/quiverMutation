# Skeptic, round 040 (conference)

No new work. Read: persona, notebook, STATE.md, DIGEST.md (rounds 036-039), STEERING.md.

## Most promising question for the next few rounds

Does any gate-admitted step with J != 0 (silting, not tilting) have a child that keeps the LNA Coxeter key, when the parent is reachable from an LNA by key-preserving walks?

Why: E-138 says no at n <= 8 (all 229 J != 0 rows refused by the key guard), and E-137 says the E-134 meeting needs exactly such a step. If the answer is yes, the Coxeter guard (H-015) is not the tilting test on walks, and the "16 fans lie in the LNA derived class" claim fails. If no, the guard becomes a checkable theorem-shaped statement. Either way it is a bounded, decidable question at n = 6, 7, 8, and it decides which of the two readings of E-134 survives.

## Weakest claim the workshop relies on

The E-134 membership claim: the 16 n = 6 fans with an LNA key are in the LNA derived class. It rests on a key-preserving meet whose hit-side first step is gate-admitted, key-preserving and fails tiltingPlus and Cartan congruence (E-137). With non-tilting steps removed, 0 of 16 hits meet the LNA side, but that search is not closed (LNA side 27.5k algebras in 150 s). So the claim is unsupported, not refuted. Downstream: E-135/E-136 and the round-036 item 1 framing lean on it.

Secondary: "d_i = 2 whenever J_i != 0" (E-129). It holds only on capped walks with out-degree <= 2; E-131 has counterexamples (d_i = 4, 4, 5) past the cap. It should not appear in any summary without the cap.

## What I need from other personas

- toolsmith: the tilting-only hit closure with a positive control and a closure flag, plus a count of how many steps the closure rejects for tiltingPlus. Without the control, a null means nothing.
- theorist: whether the Coxeter guard, as coded, admits the hit -> M step (`coxkey-pres = True` in E-137). If it does, the guard is the defect; if not, the E-134 path was never a guard-admitted path.
- scholar: Ladkani 2.3(c) and AI 2.32(b) read against the tiltingPlus code for one non-monomial step, so that tiltingPlus is not checked only against itself.
- experimentalist: the n = 8 c1, c2 J != 0 key-guard column (already requested in E-138) so the answer to the question above covers more than class 0.

## Proposed (not run)

Replay the 13 remaining E-134 hits with tiltingPlus AND Cartan congruence, not the guard that generated the path, and print the shared keys. This is the same check the notebook lists under "Habits", applied to all hits.

## Summary line

Most promising question: does any gate-admitted J != 0 step with a key-preserving parent keep the LNA Coxeter key?
