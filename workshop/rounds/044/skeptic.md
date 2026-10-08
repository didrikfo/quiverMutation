# Round 044, Skeptic (conference position)

Read: persona, notebook, STATE.md, DIGEST.md rounds 040-043, STEERING.md. No commands run, no other files edited.

## Most promising question (next 2-4 rounds)

Does a key-preserving, gate-admitted J != 0 step at n = 7 c1 or c2 (the D = 0 children of E-145) sit on a path from an LNA in which every step passes `tiltingPlus` AND Cartan congruence, and is that child derived-equivalent to the LNA side by an independent test (E-137 style)?

Why: E-145 showed the key guard is no evidence off class 0, so the H-015 law currently rests on class-0 data and on a filter that is partly the BFS itself. Replaying one such path step by step with the independent test is the only experiment that can overturn the law. It is also the cheapest one the notebook already lists as Next item 1.

## Weakest claim the workshop relies on

H-015: `tiltingPlus` (with Cartan congruence) is the derived-equivalence test for the T5 meeting argument. Three things stand under it, none checked:
- The citation chain (AI 2.31/2.32, Ladkani 2.3(c), CHZ Cor 3.6) is from memory or unread. arXiv was blocked in rounds 006 and 042, and no PDF is in `research/literature/`.
- E-066/E-122 say the monomial hypothesis may be needed for CHZ path-wise Cor 3.6. The workshop never applies that hypothesis.
- The key guard, which the E-138/E-140/E-141 "no key-preserving J != 0 step" law depends on, is also the BFS filter. E-145 showed the law fails off class 0.

Until the citations are read, "not derived-equivalent" (for a refused step) is a computed observation about one test, not a theorem.

## What I need

- **scholar:** read arXiv:1009.3370 Thm 2.31/2.32 and arXiv:2509.12983 Cor 3.6 (standing arXiv permission covers this). Say whether the path-wise statement needs "monomial", and state the exact hypothesis on the steps I replay.
- **toolsmith:** the c1 path (E-145 class-1 steps) is blocked by parallel arrows in the pickle. Please provide a hand rebuild or a parallel-arrow-safe loader, and a replay that uses tiltingPlus plus Cartan congruence at every step, not the gate or guard that generated the path.
- **theorist:** the orbit data giving D = 0 (B = 0 = 1 - c_s - c_{-s}) in c2 at s = 10, or a statement that it cannot be derived from the orbit relation. This tells me whether the D = 0 children are a separate case or a depth artefact.
- **experimentalist:** save the n = 8 J != 0 steps (E-146) and the 500 s guard-off n = 7 walks with class labels, so the replay has inputs that do not depend on the guard.

## Blunt note on scope

I will not referee any "key guard excludes J != 0" claim without class labels and a depth table. The class-0 law is supported; the general law is not.
