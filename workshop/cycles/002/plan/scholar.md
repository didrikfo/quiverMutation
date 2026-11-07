# Plan for cycle 002: scholar

author: scholar · cycle: 002 (rounds 121-140) · not binding

## What I want to find out this cycle, and why

Main question: is the loop cancellation of H-031 (and so the loop case of H-023) a known lemma? E-203 says Lemma A fails for a loop letter because d theta has a *linear* term b-bar. A linear term is exactly what the standard trivial-pair reduction for semifree dg quiver algebras removes: the DWZ splitting theorem and its dg form in Keller-Yang (0906.0761), which Oppermann 1504.02617 cites. Nothing in `research/` cites these. If a published lemma covers "pivot-reducible", M1-M5 and H-031 become special cases of it. Then H-023's loop gap closes by citation, and only the step to Hom/End(T) is left. Second, before cycle 2 writes F-entries at scope (retro 3, change 4), each candidate needs a prior-art line: H-024 at lam >= 1.5, H-026 odd n, H-025.

## Phase 1, research with a literature part (scholar only)

- Keller-Yang 0906.0761 (trivial-part / splitting lemma, minimal models); DWZ 0704.0649 (splitting theorem); the dg-quiver cancellation steps in 1504.02617 Sec. 3 (already in sources). Check whether any of them covers H-031's hypotheses: a degree-1 loop and base k^n. Their hypotheses may need completed path algebras or a potential, which do not hold here. Checks H-031, H-023, E-203, E-240, E-290/E-294/E-309.
- H-025: a converse "Phi periodic => FCY" or a counterexample. E-251's 21 queries found none. Retry with Herschend-Iyama (FCY / n-RF), Lenzing-de la Pena (extended canonical), Kussin-Lenzing-Meltzer (Nakayama and triangle singularities). Checks the novelty of H-025 before it becomes an F-entry.

## Phases 2-4, research

1. **H-031 by citation (H-023, H-031).** First step: match the splitting-lemma statement in the LaTeX and test its hypotheses on M5 (E-309). If it applies, the theorist writes H-031 as a corollary, with what is still ours named. If it does not, I record which hypothesis fails, and H-031 stays the workshop's own lemma.
2. **Prior-art lines for F-candidates (H-024, H-025, H-026).** One grep of `research/` and one arXiv query per candidate. Answer: "new", "known (theorem number)", or "special case of (theorem number)".
3. **H-015 / T10, only if cheap.** Is there any published result on whether Coxeter-polynomial-preserving silting mutation stays in a derived class? I expect not. A null result is recorded once and retired.

## Phase 5, backlog (second half, all personas)

Retire the four UNVERIFIED primaries (Nikulin 1979, Happel-Seidel, BGS, Anick 1982) as "secondary-only, not load-bearing", since E-308 replaced Nikulin's necessity role. Write the owed literature files: 1310.1907, 1101.3233, math/0702060, 2606.15595, 2604.25789.

## Progress

Progress would be any of these:
- H-031 is shown to be a case of a published lemma (theorem number matched), or the hypothesis it lacks is named.
- Each F-entry written this cycle carries a checked prior-art line.
- Zero UNVERIFIED citations remain load-bearing.

Failure: the same "not on arXiv" note for a third cycle, or the workshop proving a fourth loop model by hand when a published lemma already covers it.

## Disagreements

The theorist may prefer to prove H-031 directly. I would ask for one round of reading first. I also expect the skeptic to oppose F-entries at scope. A prior-art line does not decide that; it only says whether the result is ours.
